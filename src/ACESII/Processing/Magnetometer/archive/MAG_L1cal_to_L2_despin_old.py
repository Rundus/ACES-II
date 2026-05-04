# --- MAG_L1cal_to_L2_despin_old.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: use the insights from cal1_determine_ACS_timeOffsets.py to adjust the DCM
# timebase and apply the despin. Include the model field with the output




# --- bookkeeping ---
import time
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
#################
# --- IMPORTS ---
#################
from scipy.integrate import simpson
from src.ACESII.data_tools.my_imports import *
import datetime as dt
start_time = time.time()



######################
# --- DATA TOGGLES ---
######################
just_print_file_names_bool = False
rocket_str = 'low'
wInstr = 'EEPAA'
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['L2', [[0],[0]]],
    'attitude':['attitude',[[0],[0]]]
}
outputData = True

# --- OLD OFFSETS ---
TimeOffset = [0.12756724137931032, 0.12078947368421052] # apply this to the CHAOS data before interpolating onto RKT timebase
offsetDC = np.array([[0.10010010010010717,-1.2901290129012892,4.590459045904595], [0,0,0]])

# --- NEW OFFSETS ---
# from cal1_determine_ACS_timeOffsets
slope_init = {'high':-0.000006,
              'low':0.000089}
intercept_init = {'high':0.12756724137931032,
                  'low':0.12756724137931032}

find_DC_offset = False
replaceNANS = True
outputData = True # plots RKT XYZ and CHAOS model Spun-up XYZ


def RingCore_L1_to_L2_Despin(data_dicts):

    #######################
    # --- LOAD THE DATA ---
    #######################
    stl.prgMsg(f'Loading data')
    data_dict_mag = deepcopy(data_dicts[0])
    data_dict_attitude =deepcopy(data_dicts[1])
    stl.Done(start_time)


    stl.prgMsg('Despinning RingCore Data')

    # --- Define some useful variables ---
    T0 = dt.datetime(2022,11,20,17,20)
    B_rkt = np.array([data_dict_mag['Bx'][0], data_dict_mag['By'][0], data_dict_mag['Bz'][0]]).T
    T0_attitude = stl.EpochTo_T0_Rocket(data_dict_attitude['Epoch'][0], T0=T0)
    T0_ringCore = stl.EpochTo_T0_Rocket(data_dict_mag['Epoch'][0], T0=T0)

    #####################################################
    # --- [1] Calculate CHAOS model on Attitude Epoch ---
    #####################################################
    B_CHAOS_ENU_attitude = stl.CHAOS(
                    lat=data_dict_attitude['Lat'][0],
                    long=data_dict_attitude['Long'][0],
                    alt=data_dict_attitude['Alt'][0]/stl.m_to_km,
                    times=data_dict_attitude['Epoch'][0])

    ################################################
    # --- [2] Spin Up CHAOS model into RKT frame ---
    ################################################
    # --- get the DCM ---
    DCM = np.array(
        [
            [[data_dict_attitude['a11'][0][i], data_dict_attitude['a12'][0][i], data_dict_attitude['a13'][0][i]],
             [data_dict_attitude['a21'][0][i], data_dict_attitude['a22'][0][i], data_dict_attitude['a23'][0][i]],
             [data_dict_attitude['a31'][0][i], data_dict_attitude['a32'][0][i], data_dict_attitude['a33'][0][i]]]
            for i in range(len(data_dict_attitude['Epoch'][0]))
        ]
    )
    B_CHAOS_rkt = np.array([np.matmul(DCM[i].T, B_CHAOS_ENU_attitude[i]) for i in range(len(data_dict_attitude['Epoch'][0]))])

    ##############################################
    # --- [3] interpolate CHAOS onto mag Epoch ---
    ##############################################
    # Adjust the timebase
    T0_attitude_adjusted = T0_attitude + TimeOffset[wRocket - 4]
    B_CHAOS_rkt_magTime = np.zeros_like(B_rkt)

    for i in range(3):
        cs = CubicSpline(T0_attitude_adjusted, B_CHAOS_rkt[:, i])
        B_CHAOS_rkt_magTime[:,i] = cs(T0_ringCore)

    # adjust the magnitude
    B_CHAOS_rkt_magTime = B_CHAOS_rkt_magTime + offsetDC[wRocket-4]

    #######################################
    # --- [4] De-Spin the RingCore Data ---
    #######################################
    # Get the DCM on the RKT timebase
    DCM_mag = np.zeros(shape=(len(data_dict_mag['Epoch'][0]),3,3))
    for i in range(3):
        for j in range(3):
            cs = CubicSpline(T0_attitude, data_dict_attitude[f'a{i+1}{j+1}'][0])
            DCM_mag[:,i,j] = cs(T0_ringCore)

    # Apply DCM to RingCore Data
    B_ringCore_ENU = np.array([np.matmul(DCM_mag[i], B_rkt[i]) for i in range(len(T0_ringCore))])

    #################################################
    # --- [5] Get the Time-adjusted DCM for CHAOS ---
    #################################################
    DCM_timeAdjusted = np.zeros_like(DCM_mag)
    for i in range(3):
        for j in range(3):
            cs = CubicSpline(T0_attitude_adjusted, data_dict_attitude[f'a{i+1}{j+1}'][0])
            DCM_timeAdjusted[:, i, j] = cs(T0_ringCore)

    # apply DCM to CHAOS data
    B_CHAOS_ENU = np.array([np.matmul(DCM_timeAdjusted[i], B_CHAOS_rkt_magTime[i]) for i in range(len(T0_ringCore))])

    #######################################
    # --- [6] Handle NAN values in data ---
    #######################################
    if replaceNANS:
        def nan_helper(y):
            return np.isnan(y), lambda  z: z.nonzero()[0]

        newB = []

        # linearily interpolate the Nans
        for i in range(3):
            Bcomp =B_ringCore_ENU[:, i]
            nans, x = nan_helper(Bcomp)
            Bcomp[nans] = np.interp(x(nans),x(~nans),Bcomp[~nans])
            newB.append(Bcomp)

        B_ringCore_ENU = np.array([ [newB[0][i],newB[1][i],newB[2][i]] for i in range(len(B_ringCore_ENU))])

    ###############################
    # --- Determine the DeltaB  ---
    ###############################
    B_mag = np.array([np.linalg.norm(vec) for vec in B_ringCore_ENU])
    B_ringCore_ENU = B_ringCore_ENU - B_CHAOS_ENU
    stl.Done(start_time)

    if outputData:

        stl.prgMsg('Writing out Data')

        data_dict_output = {**{},
                            **{
                                'B_E' : [B_ringCore_ENU[:, 0], deepcopy(data_dict_mag['Bx'][1])],
                                'B_N': [B_ringCore_ENU[:, 1], deepcopy(data_dict_mag['Bx'][1])],
                                'B_U': [B_ringCore_ENU[:, 2], deepcopy(data_dict_mag['Bx'][1])],
                                'B_model_E': [B_CHAOS_ENU[:, 0], deepcopy(data_dict_mag['Bx'][1])],
                                'B_model_N': [B_CHAOS_ENU[:, 1], deepcopy(data_dict_mag['Bx'][1])],
                                'B_model_U': [B_CHAOS_ENU[:, 2], deepcopy(data_dict_mag['Bx'][1])],
                                'Bmag': [B_mag,deepcopy(data_dict_mag['Bmag'][1])],
                                'Epoch':deepcopy(data_dict_mag['Epoch'])
                            }}

        for key in ['B_E','B_N','B_U','B_model_E','B_model_N','B_model_U']:
            data_dict_output[f'{key}'][1]['LABLAXIS'] = key

        # Write out the File
        fileoutName_despin = f'ACESII_{ACESII.payload_IDs[rocket_str]}_l2_RingCore_ENU_fullCal'
        outputPath = f'{DataPaths.ACES_data_folder}/L2/{rocket_str}//{fileoutName_despin}.cdf'
        stl.outputCDFdata(outputPath, data_dict_output, instrNam='RingCore')
        stl.Done(start_time)



#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(RingCore_L1_to_L2_Despin,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)
