# --- MPI_txt_to_L3.py ---
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
def MPI_txt_to_L3(wflyer, justPrintFileNames):

    # --- FILE I/O ---
    rocket_folder_path = DataPaths.ACES_data_folder
    rocketID = ACESII.payload_IDs[wflyer]

    # Set the paths for the file names
    data_repository = f'{rocket_folder_path}{input_path_modifier}\\{ACESII.fliers[wflyer]}\\'
    input_files = glob(data_repository + '*.txt')
    input_names = [ifile.replace(data_repository, '') for ifile in input_files]
    input_names_searchable = [ifile.replace(input_path_modifier.lower() + '_', '').replace('_v00', '') for ifile in input_names]

    if justPrintFileNames:
        for i, file in enumerate(input_files):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(os.path.getsize(file) / (10 ** 6), 1)))
        return

    print('\n')
    print(stl.color.UNDERLINE + f'Converting to L3 cdf data' + stl.color.END)

    # --- prepare the output ---
    data_dict_output = {}

    # --- get the data from the .txt files ---
    file_data = []
    for idx, file in enumerate(input_files):
        with io.open(file, mode="r", encoding="utf-8") as f:
            next(f)
            temp_data = []
            for line in f:
                temp_data.append([float(val) for val in line.split()])

            temp_data = np.array(temp_data).T

            # correct the output from T0 = 1900-01-01
            T0 = dt.datetime(1900,1,1)
            temp_epoch_us = temp_data[0]*1E6 # convert to micro seconds
            Epoch = np.array([T0+ dt.timedelta(microseconds=val) for val in temp_epoch_us])
            data_dict_output = {**data_dict_output, **{f'MPI{idx+1}_Epoch': [Epoch,{'LABLAXIS': f'MPI{idx+1}_Epoch', 'DEPEND_0': None,
                                                                   'DEPEND_1': None,
                                                                   'DEPEND_2': None,
                                                                   'FILLVAL': ACESII.epoch_fillVal,
                                                                   'FORMAT': 'E12.2',
                                                                   'UNITS': 'ns',
                                                                   'VALIDMIN': temp_data[0].min(),
                                                                   'VALIDMAX': temp_data[0].max(),
                                                                   'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}],
                                f'MPI{idx+1}_Vx': [temp_data[1],{'LABLAXIS': f'MPI{idx+1}_Vx', 'DEPEND_0': f'MPI{idx+1}_Epoch',
                                                                   'DEPEND_1': None,
                                                                   'DEPEND_2': None,
                                                                   'FILLVAL': ACESII.epoch_fillVal,
                                                                   'FORMAT': 'E12.2',
                                                                   'UNITS': 'm/s',
                                                                   'VALIDMIN': 0,
                                                                   'VALIDMAX': 0,
                                                                   'VAR_TYPE': 'data', 'SCALETYP': 'linear'}],
                                f'MPI{idx+1}_Vy': [temp_data[2],{'LABLAXIS': f'MPI{idx+1}_Vy', 'DEPEND_0': f'MPI{idx+1}_Epoch',
                                                                   'DEPEND_1': None,
                                                                   'DEPEND_2': None,
                                                                   'FILLVAL': ACESII.epoch_fillVal,
                                                                   'FORMAT': 'E12.2',
                                                                   'UNITS': 'm/s',
                                                                   'VALIDMIN': 0,
                                                                   'VALIDMAX': 0,
                                                                   'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]
                                }
                                }


    ##################################################################
    # --- Calculate the plasma density from Ion Saturation Current ---
    ##################################################################

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        stl.prgMsg('Creating output file')
        fileoutName_fixed = f'ACESII_{rocketID}_l3_MPI_rktFrm.cdf'
        outputPath = data_repository + f'\\{fileoutName_fixed}'
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
target_files = f'{DataPaths.ACES_data_folder}{input_path_modifier}{ACESII.fliers[wRocket-4]}\*.txt'

if len(glob(target_files)) == 0:
    print(stl.color.RED + 'There are no .txt files in the specified directory' + stl.color.END)
else:
    MPI_txt_to_L3(wRocket-4, justPrintFileNames)

