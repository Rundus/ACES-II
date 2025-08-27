# --- L1_to_L1_rocket_convection_cal.py ---
# Description: Calculate the rocket convection terms to the EFI dataset
# THEN rotate everything into auroral coordinates



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

from src.my_imports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
from scipy.interpolate import CubicSpline

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 4 -> ACES-II High Flier
# 5 -> ACES-II Low Flier
wRocket = 5
wFiles = [0]
outputData = True
use_RingCore_data = True # if False, just use CHAOS data
Plot_correction_term = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none

def L1_to_L1_rocket_convection_cal(wRocket, justPrintFileNames):
    inputFiles_elec = glob(f'{DataPaths.ACES_data_folder}\\L1\\{ACESII.fliers[wRocket-4]}\*E_Field_ENU_no_corrections*')
    inputFiles_mag = glob(f'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket-4]}\\*RingCore_ENU.cdf*')
    input_names = [ifile.replace(f'{DataPaths.ACES_data_folder}\\L1\\{ACESII.fliers[wRocket-4]}\\', '') for ifile in inputFiles_elec]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles_elec):
            print('[{:.0f}] {:80s}'.format(i, input_names[i], round(getsize(file) / (10 ** 6))))
        return

    # --- get the data from the B-Field file ---
    stl.prgMsg(f'Loading data')
    data_dict_mag = stl.loadDictFromFile(inputFiles_mag[0])
    data_dict_EFI, GlobalAttrs = stl.loadDictFromFile(inputFiles_elec[0], getGlobalAttrs=True)
    data_dict_traj = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\trajectories\\{ACESII.fliers[wRocket-4]}\\*_ECEF.cdf*')[0])
    stl.Done(start_time)

    ########################################
    # --- ADD IN ROGER'S TIME CORRECTION ---
    ########################################
    stl.prgMsg('Adding In Rogers timebase correction')
    timeCorrection = (0.1157 * 1E6) # in us
    data_dict_EFI['Epoch'][0] = np.array([tme + dt.timedelta(microseconds=timeCorrection) for tme in data_dict_EFI['Epoch'][0]])
    stl.Done(start_time)

    #########################################################
    # --- interpolate E-Field data onto RingCore TimeBase ---
    #########################################################
    stl.prgMsg('Interpolating E-Field Data')
    Epoch_EFI_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_EFI['Epoch'][0]])
    Epoch_mag_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_mag['Epoch'][0]])

    for key in data_dict_EFI.keys():
        if key != 'Epoch':
            cs = CubicSpline(Epoch_EFI_tt2000,data_dict_EFI[key][0])
            data_dict_EFI[key][0] = deepcopy(cs(Epoch_mag_tt2000))

    data_dict_EFI['Epoch'][0] = deepcopy(data_dict_mag['Epoch'][0])
    stl.Done(start_time)


    ################################################
    # --- Calculate the Rocket Convection Effect ---
    ################################################
    stl.prgMsg('Correcting Rocket Convection Effect')

    # Get the payload velocity vector
    rkt_VEL_ECEF = np.array([data_dict_traj['ECEFXVEL'][0], data_dict_traj['ECEFYVEL'][0], data_dict_traj['ECEFZVEL'][0]]).T

    # Convert payload velocity to ENU
    data_dict_transform_ENU = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\transforms\\{ACESII.fliers[wRocket-4]}\\*ECEF_to_ENU.cdf*')[0])

    rkt_VEL_ENU = np.array([np.matmul(data_dict_transform_ENU['ECEF_to_ENU'][0][i], vec) for i,vec in enumerate(rkt_VEL_ECEF)])

    # Get the CHAOS geomagnetic field
    B_model = 1E-9 * stl.CHAOS(lat=data_dict_traj['Lat'][0],
                               long=data_dict_traj['Long'][0],
                               alt=data_dict_traj['Alt'][0],
                               times=data_dict_traj['Epoch'][0])  # CHAOS in ENU coordinates

    if use_RingCore_data:
        # Load the RingCore - assumes the RingCore data does NOT have the CHAOS model included
        data_dict_RingCore = stl.loadDictFromFile(glob(rf'C:\Data\ACESII\L2\{ACESII.fliers[wRocket-4]}\\*RingCore_ENU.cdf*')[0])

        # Interpolate RingCore Data onto Trajectory timebase
        T0_RingCore = stl.EpochTo_T0_Rocket(data_dict_RingCore['Epoch'][0],T0=dt.datetime(2022,11,20,17,20))
        T0_Traj = stl.EpochTo_T0_Rocket(data_dict_traj['Epoch'][0], T0=dt.datetime(2022, 11, 20, 17, 20))
        for key in ['B_East', 'B_North', 'B_Up']:
            cs = CubicSpline(T0_RingCore, data_dict_RingCore[f'{key}'][0])
            data_dict_RingCore[f'{key}'][0] = cs(T0_Traj)

        B_vec_Ringcore = 1E-9*np.array([data_dict_RingCore['B_East'][0], data_dict_RingCore['B_North'][0], data_dict_RingCore['B_Up'][0]]).T
        B_model = B_model + B_vec_Ringcore

    # Calculate the vxB electric field in ENU
    vxB_term = np.array([np.cross(vec, B_model[i]) for i, vec in enumerate(rkt_VEL_ENU)])

    # plot the calibration term
    if Plot_correction_term:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3)

        ax[0].plot(data_dict_EFI['Epoch'][0], data_dict_EFI['E_East'][0],color='blue')
        ax[0].plot(data_dict_transform_ENU['Epoch'][0], vxB_term[:,0], color='red')

        ax[1].plot(data_dict_EFI['Epoch'][0], data_dict_EFI['E_North'][0], color='blue')
        ax[1].plot(data_dict_transform_ENU['Epoch'][0], vxB_term[:,1], color='red')

        ax[2].plot(data_dict_EFI['Epoch'][0], data_dict_EFI['E_Up'][0], color='blue')
        ax[2].plot(data_dict_transform_ENU['Epoch'][0], vxB_term[:,2], color='red')

        for i in range(3):
            ax[i].set_ylim(-0.12,0.12)

        plt.show()
    stl.Done(start_time)

    # interpolate the -vxB term onto the magnetometer timebase
    Epoch_traj_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_traj['Epoch'][0]])

    # X
    cs = CubicSpline(Epoch_traj_tt2000,vxB_term[:,0])
    vxB_E = cs(Epoch_mag_tt2000)

    # Y
    cs = CubicSpline(Epoch_traj_tt2000, vxB_term[:, 1])
    vxB_N = cs(Epoch_mag_tt2000)

    # Z
    cs = CubicSpline(Epoch_traj_tt2000, vxB_term[:, 2])
    vxB_U = cs(Epoch_mag_tt2000)

    # form the new, downsampled vectors
    vxB_ENU = np.array([vxB_E,vxB_N,vxB_U]).T

    # --- prepare the E-Field ---
    E_Field_ENU = np.array([deepcopy(data_dict_EFI['E_East'][0]), data_dict_EFI['E_North'][0], data_dict_EFI['E_Up'][0]]).T

    ########################################
    # --- Convert to auroral coordiantes ---
    ########################################
    stl.prgMsg('Transforming Coordinates')
    data_dict_transform_auroral = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\transforms\\{ACESII.fliers[wRocket - 4]}\\*ECEF_to_auroral.cdf*')[0])

    # Interpolate the transformation matrices onto the magnetometer timebase
    ENU_to_ECEF_transform = np.zeros(shape=(len(vxB_ENU), 3, 3))
    Epoch_tt2000_ENU_to_ECEF_trans = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_transform_ENU['Epoch'][0]])
    for i in range(3):
        for j in range(3):
            cs = CubicSpline(Epoch_tt2000_ENU_to_ECEF_trans,data_dict_transform_ENU[f'a{i+1}{j+1}'][0])
            ENU_to_ECEF_transform[:,i,j] = cs(Epoch_mag_tt2000)

    ENU_to_ECEF_transform = np.array([val.T for val in ENU_to_ECEF_transform])

    # Interpolate the transformation matrices onto the magnetometer timebase
    ECEF_to_auroral_transform = np.zeros(shape=(len(vxB_ENU), 3, 3))
    Epoch_tt2000_ECEF_to_auroral_trans = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_transform_auroral['Epoch'][0]])
    for i in range(3):
        for j in range(3):
            cs = CubicSpline(Epoch_tt2000_ECEF_to_auroral_trans, data_dict_transform_auroral[f'a{i+1}{j+1}'][0])
            ECEF_to_auroral_transform[:, i, j] = cs(Epoch_mag_tt2000)

    vxB_ECEF = np.array([np.matmul(ENU_to_ECEF_transform[i], vec) for i, vec in enumerate(vxB_ENU)])
    vxB_auroral = np.array([np.matmul(ECEF_to_auroral_transform[i], vec) for i, vec in enumerate(vxB_ECEF)])

    E_Field_ECEF = np.array([np.matmul(ENU_to_ECEF_transform[i], vec) for i, vec in enumerate(E_Field_ENU)])
    E_Field_auroral = np.array([np.matmul(ECEF_to_auroral_transform[i], vec) for i, vec in enumerate(E_Field_ECEF)])

    stl.Done(start_time)

    data_dict_output = {'vxB_N':[vxB_auroral[:,0], {'LABLAXIS': 'vxB_N', 'DEPEND_0': 'Epoch',
                                                                'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': ACESII.epoch_fillVal,
                                                                'FORMAT': 'E12.2', 'UNITS': 'V/m',
                                                                'VAR_TYPE': 'data',
                                                                'SCALETYP': 'linear'}],
                           'vxB_T': [vxB_auroral[:,1],{'LABLAXIS': 'vxB_T', 'DEPEND_0': 'Epoch',
                                                                'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': ACESII.epoch_fillVal,
                                                                'FORMAT': 'E12.2', 'UNITS': 'V/m',
                                                                'VAR_TYPE': 'data',
                                                                'SCALETYP': 'linear'}],
                           'vxB_p': [vxB_auroral[:,2],{'LABLAXIS': 'vxB_p', 'DEPEND_0': 'Epoch',
                                                                'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': ACESII.epoch_fillVal,
                                                                'FORMAT': 'E12.2', 'UNITS': 'V/m',
                                                                'VAR_TYPE': 'data',
                                                                'SCALETYP': 'linear'}],
                            'Epoch': deepcopy(data_dict_EFI['Epoch']),

                            'E_N_raw': [E_Field_auroral[:,0],deepcopy(data_dict_EFI['E_East'][1])],
                            'E_T_raw': [E_Field_auroral[:,1],deepcopy(data_dict_EFI['E_North'][1])],
                            'E_p_raw': [E_Field_auroral[:,2],deepcopy(data_dict_EFI['E_Up'][1])],
                           }

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')

        fileoutName = rf'ACESII_{ACESII.payload_IDs[wRocket-4]}_l1_E_Field_ENU_rkt_convection.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\\calibration\\EFI_rkt_convection_calibration\\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output, globalAttrsMod=GlobalAttrs, instrNam='EFI')
        stl.Done(start_time)

# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
L1_to_L1_rocket_convection_cal(wRocket, justPrintFileNames)