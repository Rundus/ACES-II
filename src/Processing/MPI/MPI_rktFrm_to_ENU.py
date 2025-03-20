# --- MPI_rktFrm_to_ENU.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

from src.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False # Just print the names of files

# --- Select the Rocket ---
# 5 -> ACES II Low Flier - ONLY ONE AVAILABLE
wRocket = 5

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = []

# --- OutputData ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.my_imports import *
import io
input_path_modifier = 'L3\\MPI\\'


#######################
# --- MAIN FUNCTION ---
#######################
def MPI_rktFrm_to_ENU(wflyer, justPrintFileNames):

    # --- FILE I/O ---
    rocket_folder_path = DataPaths.ACES_data_folder
    rocketID = ACESII.payload_IDs[wflyer]

    # Set the paths for the file names
    data_repository = f'{rocket_folder_path}{input_path_modifier}\\{ACESII.fliers[wflyer]}\\'
    input_files = glob(data_repository + '*.cdf')
    input_names = [ifile.replace(data_repository, '') for ifile in input_files]
    input_names_searchable = [ifile.replace(input_path_modifier.lower() + '_', '').replace('_v00', '') for ifile in input_names]

    if justPrintFileNames:
        for i, file in enumerate(input_files):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(os.path.getsize(file) / (10 ** 6), 1)))
        return

    print('\n')
    print(stl.color.UNDERLINE + f'Converting to L3 cdf data' + stl.color.END)

    #######################
    # --- LOAD THE DATA ---
    #######################
    data_dict_attitude = stl.loadDictFromFile(r'C:\Data\ACESII\attitude\low\ACESII_36364_Attitude_Solution.cdf')
    data_dict_MPI = stl.loadDictFromFile(input_files[0])

    # --- prepare the output ---
    no_of_MPIs = 4
    data_dict_output = {f'MPI{j+1}_Epoch':deepcopy(data_dict_MPI[f'MPI{j+1}_Epoch']) for j in range(no_of_MPIs)}


    ##################################################################
    # --- Calculate the plasma density from Ion Saturation Current ---
    ##################################################################
    stl.prgMsg('Interpolating Attitude Data onto MPI')

    # Step 1: Interpolate Attitude DCM data into MPI datasets
    DCMs = []
    for i in range(no_of_MPIs):
        target_epoch = data_dict_MPI[f'MPI{i+1}_Epoch'][0]
        data_dict_attitude_interp = stl.InterpolateDataDict(InputDataDict=deepcopy(data_dict_attitude),
                                InputEpochArray=deepcopy(data_dict_attitude['Epoch'][0]),
                                wKeys=['a11', 'a12', 'a13', 'a21', 'a22', 'a23', 'a31', 'a32', 'a33', 'Epoch'],
                                targetEpochArray=target_epoch,)

        DCMs.append(np.array(
            [
                [[data_dict_attitude_interp['a11'][0][i], data_dict_attitude_interp['a12'][0][i], data_dict_attitude_interp['a13'][0][i]],
                 [data_dict_attitude_interp['a21'][0][i], data_dict_attitude_interp['a22'][0][i], data_dict_attitude_interp['a23'][0][i]],
                 [data_dict_attitude_interp['a31'][0][i], data_dict_attitude_interp['a32'][0][i], data_dict_attitude_interp['a33'][0][i]]]
                for i in range(len(data_dict_attitude_interp['Epoch'][0]))
            ]
        )
        )

    stl.Done(start_time)

    ########################
    # --- APPLY THE DCMs ---
    ########################
    stl.prgMsg('Rotating into ENU coordinates')

    for idx in range(no_of_MPIs):

        # form the velocity vector
        v_vec_rktFrm = np.array([[data_dict_MPI[f'MPI{i+1}_Vx'][0][j], data_dict_MPI[f'MPI{i+1}_Vy'][0][j], 0] for j in range(len(data_dict_MPI[f'MPI{idx+1}_Epoch'][0]))])
        wDCM = DCMs[idx]

        # rotate the velocities
        v_vec_ENU = np.array([np.matmul(wDCM[i], v_vec_rktFrm[i]) for i in range(len(v_vec_rktFrm))])

        data_dict_output = {**data_dict_output,
                            **{f'MPI{idx + 1}_VE': [v_vec_ENU[:,0],
                                                    {'LABLAXIS': f'MPI{idx + 1}_VE', 'DEPEND_0': f'MPI{idx + 1}_Epoch',
                                                     'DEPEND_1': None,
                                                     'DEPEND_2': None,
                                                     'FILLVAL': ACESII.epoch_fillVal,
                                                     'FORMAT': 'E12.2',
                                                     'UNITS': 'm/s',
                                                     'VALIDMIN': v_vec_ENU[:,0].min(),
                                                     'VALIDMAX': v_vec_ENU[:,0].max(),
                                                     'VAR_TYPE': 'data', 'SCALETYP': 'linear'}],
                               f'MPI{idx + 1}_VN': [v_vec_ENU[:,1],
                                                    {'LABLAXIS': f'MPI{idx + 1}_VN', 'DEPEND_0': f'MPI{idx + 1}_Epoch',
                                                     'DEPEND_1': None,
                                                     'DEPEND_2': None,
                                                     'FILLVAL': ACESII.epoch_fillVal,
                                                     'FORMAT': 'E12.2',
                                                     'UNITS': 'm/s',
                                                     'VALIDMIN': v_vec_ENU[:,1].min(),
                                                     'VALIDMAX': v_vec_ENU[:,1].max(),
                                                     'VAR_TYPE': 'data', 'SCALETYP': 'linear'}],
                               f'MPI{idx + 1}_VU': [v_vec_ENU[:, 2],
                                                    {'LABLAXIS': f'MPI{idx + 1}_VU', 'DEPEND_0': f'MPI{idx + 1}_Epoch',
                                                     'DEPEND_1': None,
                                                     'DEPEND_2': None,
                                                     'FILLVAL': ACESII.epoch_fillVal,
                                                     'FORMAT': 'E12.2',
                                                     'UNITS': 'm/s',
                                                     'VALIDMIN': v_vec_ENU[:, 2].min(),
                                                     'VALIDMAX': v_vec_ENU[:, 2].max(),
                                                     'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]
                               }
                            }

    stl.Done(start_time)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        stl.prgMsg('Creating output file')
        fileoutName_fixed = f'ACESII_{rocketID}_l3_MPI_ENU.cdf'
        outputPath = data_repository + f'\\{fileoutName_fixed}'
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
target_files = f'{DataPaths.ACES_data_folder}{input_path_modifier}{ACESII.fliers[wRocket-4]}\*.cdf'

if len(glob(target_files)) == 0:
    print(stl.color.RED + 'There are no .txt files in the specified directory' + stl.color.END)
else:
    MPI_rktFrm_to_ENU(wRocket-4, justPrintFileNames)

