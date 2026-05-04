# --- cal1_determine_ACS_timeOFfsets.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: adjust the timebase of the ACS DCM to see if a cleaner MAG despin can be derived. Also remove any NaNs from the data

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2026-03-30"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np

# --- --- --- --- ---


######################
# --- DATA TOGGLES ---
######################
just_print_file_names_bool = False
rocket_str = 'low'
wInstr = 'MAG'
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['L1', [[0],[0]]],
    'attitude':['',[[0],[0]]],
}
outputData = True

slope_cal = {'high':-0.000006,
              'low':0.000089}
intercept_cal = {'high':0.12756724137931032,
                  'low':0.12756724137931032}

#################
# --- IMPORTS ---
#################
from src.ACESII.data_tools.my_imports import *
import datetime as dt
import ppigrf
from matplotlib.widgets import Slider
from geopy import distance
start_time = time.time()

def MAG_L1_to_L2_despin(data_dicts):

    #######################
    # --- LOAD THE DATA ---
    #######################
    stl.prgMsg(f'Loading data')
    data_dict_MAG = deepcopy(data_dicts[0])
    data_dict_attitude = deepcopy(data_dicts[1])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
                       }
    data_dict_temp = {
        f'a{i}{j}':[[],{}] for i in range(1,4) for j in range(1,4)
    }
    data_dict_temp = {**data_dict_temp,
                     **{'Lat': [[], {}],
                        'Alt': [[], {}],
                        'Long': [[], {}], }
                     }

    # --- --- --- --- --- --- ---
    # --- ADJUST THE ACS TIME ---
    # --- --- --- --- --- --- ---

    # --- [0] Calculate the time lag due to speed of light ---
    # The IGNORES any curvature due to the Earth, which might skew results
    Alat = 69.294167
    Alon = 16.020833
    Aalt = 12.7  # in meters
    pt1 = np.array([[Alat, Alon, Aalt] for i in range(len(data_dict_attitude['Epoch'][0]))])
    pt2 = np.array([data_dict_attitude['Lat'][0], data_dict_attitude['Long'][0], data_dict_attitude['Alt'][0]]).T
    distance_2d = np.array([distance.distance(pt1[i][:2], pt2[i][:2]).m for i in range(len(data_dict_attitude['Epoch'][0]))])
    distance_3d = np.sqrt(np.square(distance_2d) + np.square(pt1[:, 2] - pt2[:, 2])) / stl.m_to_km
    light_time_delay = distance_3d * stl.m_to_km / stl.lightSpeed  # in seconds

    # --- [1] Adjust/inteprolate the DCM timebase with the light delay+calculated slope ---
    Bkeys = ['Bx', 'By', 'Bz']
    Bkeys_new = ['B_E', 'B_N', 'B_U']
    MAGdata = np.array([data_dict_MAG[key][0] for key in Bkeys]).T
    T0 = dt.datetime(2022,11,20,17,20)
    T0_MAG = stl.EpochTo_T0_Rocket(data_dict_MAG['Epoch'][0],T0=T0)
    T0_attitude = stl.EpochTo_T0_Rocket(data_dict_attitude['Epoch'][0], T0=T0) - light_time_delay
    T0_attitude = (1+slope_cal[rocket_str])*T0_attitude + intercept_cal[rocket_str]

    for i in range(1,4):
        for j in range(1,4):
            DCMkey = f'a{i}{j}'
            data_dict_temp[DCMkey][0] = np.interp(T0_MAG, T0_attitude,data_dict_attitude[DCMkey][0])
    for akey in ['Lat','Long','Alt']:
        data_dict_temp[akey][0] = np.interp(T0_MAG,T0_attitude,data_dict_attitude[akey][0])

    # --- [2] Apply the DCM to despin the MAG dta ---
    DCM = np.array([
        [
            [data_dict_temp['a11'][0][i],data_dict_temp['a12'][0][i],data_dict_temp['a13'][0][i]],
            [data_dict_temp['a21'][0][i],data_dict_temp['a22'][0][i],data_dict_temp['a23'][0][i]],
            [data_dict_temp['a31'][0][i],data_dict_temp['a32'][0][i],data_dict_temp['a33'][0][i]],
        ]
        for i in range(len(data_dict_MAG['Epoch'][0]))
    ])

    # apply the DCM
    B_ENU = np.array([np.matmul(DCM[i],MAGdata[i]) for i in range(len(data_dict_MAG['Epoch'][0]))])

    # --- [3] Calculate the CHAOS magnetic field for comparison ---
    stl.prgMsg('Calculating CHAOS')
    B_MODEL = stl.CHAOS(data_dict_temp['Lat'][0], data_dict_temp['Long'][0], data_dict_temp['Alt'][0] / stl.m_to_km, data_dict_MAG['Epoch'][0])
    stl.Done(start_time)

    if outputData:
        stl.prgMsg('Creating output file')

        # === [4] Construct the output data dict ===
        for key in ['Epoch','Lat','Long','Alt','L-Shell','mLat','mLong','ILat','ILong']:
            if key in data_dict_MAG.keys():
                data_dict_output = {
                    **data_dict_output,
                    **{key:deepcopy(data_dict_MAG[key])}
                }

        # === Calclate the Bmag ===
        B_mag = np.array([np.linalg.norm(vec) for vec in B_ENU])

        data_dict_output = {**data_dict_output,
                            **{
                                'B_E': [B_ENU[:, 0], deepcopy(data_dict_MAG['Bx'][1])],
                                'B_N': [B_ENU[:, 1], deepcopy(data_dict_MAG['Bx'][1])],
                                'B_U': [B_ENU[:, 2], deepcopy(data_dict_MAG['Bx'][1])],
                                'B_model_E': [B_MODEL[:, 0], deepcopy(data_dict_MAG['Bx'][1])],
                                'B_model_N': [B_MODEL[:, 1], deepcopy(data_dict_MAG['Bx'][1])],
                                'B_model_U': [B_MODEL[:, 2], deepcopy(data_dict_MAG['Bx'][1])],
                                'Bmag': [B_mag, deepcopy(data_dict_MAG['Bmag'][1])]
                            }}

        for key in ['B_E', 'B_N', 'B_U', 'B_model_E', 'B_model_N', 'B_model_U']:
            data_dict_output[f'{key}'][1]['LABLAXIS'] = key


        # remove nans in the despun data
        bad_idxs = []
        for i in range(3):
            bad_idxs.append(np.where(np.isnan(B_ENU[:,i]))[0])
        bad_idxs = np.array(list(set(np.array(bad_idxs).flatten()))) # get only the unique data in an array

        if bad_idxs.tolist() != []:
            for key in data_dict_output.keys():
                data_dict_output[key][0] = np.delete(data_dict_output[key][0],bad_idxs)

        # Write out the File
        fileoutName_despin = f'ACESII_{ACESII.fliers_dict[rocket_str]}_l2_RingCore_ENU_fullCal'
        outputPath = f'{DataPaths.ACES_data_folder}/L2/{wInstr}/{rocket_str}//{fileoutName_despin}.cdf'
        stl.outputDataDict(outputPath, data_dict_output, instrNam='RingCore')
        stl.Done(start_time)




#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(MAG_L1_to_L2_despin,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)

