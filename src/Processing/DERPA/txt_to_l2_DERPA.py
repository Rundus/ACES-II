# --- txt_to_l2_DERPA.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
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

# --- Select DERPA 1 or DERPA 2 [1/2]---
wDERPA = 1 # 1 or 2

# --- OutputData ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.my_imports import *
import io
input_path_modifier = 'science\\DERPA_raw\\'


#######################
# --- MAIN FUNCTION ---
#######################
def txt_to_l2_DERPA(wflyer, justPrintFileNames, wDERPA):

    # --- FILE I/O ---
    rocket_folder_path = DataPaths.ACES_data_folder
    rocketID = ACESII.payload_IDs[wflyer]

    # Set the paths for the file names
    data_repository = f'{rocket_folder_path}{input_path_modifier}\\{ACESII.fliers[wflyer]}\\'
    DERPA_flyers = ['HI', 'LO']
    input_files = glob(data_repository + f'*{DERPA_flyers[wRocket-4]}{wDERPA}*')
    input_names = [ifile.replace(data_repository, '') for ifile in input_files]
    input_names_searchable = [ifile.replace(input_path_modifier.lower() + '_', '').replace('_v00', '') for ifile in input_names]

    if justPrintFileNames:
        for i, file in enumerate(input_files):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(os.path.getsize(file) / (10 ** 6), 1)))
        return

    print('\n')
    print(stl.color.UNDERLINE + f'Converting to L2 cdf data' + stl.color.END)

    # --- prepare the output ---
    data_dict_output = {}

    # --- get the data from the .txt files ---
    for idx, file in enumerate(input_files):

        with io.open(file, mode="r", encoding="utf-8") as f:
            temp_data = []
            for line in f:
                temp_data.append([val for val in line.split()])

            # Handle the '?'
            temp_data = np.array(temp_data).T
            temp_data[np.where(temp_data=='?')] = 0
            temp_data = np.array(temp_data,dtype='float64')

            # correct the output from T0 = 2020-11-20
            T0 = dt.datetime(2022,11,20,17,20)
            temp_epoch_seconds = temp_data[0] # convert to micro seconds
            Epoch = np.array([T0+ dt.timedelta(seconds=val) for val in temp_epoch_seconds])

            # Store all the RAW data
            if 'peak' in file:
                data_dict_output = {**data_dict_output, **{f'Epoch_peak': [Epoch, {'LABLAXIS': f'Epoch',
                                                                              'DEPEND_0': None,
                                                                              'UNITS': None,
                                                                              'VAR_TYPE': 'data'}],
                                                           f'peak_voltage': [temp_data[1],
                                                                            {'LABLAXIS': f'Fit Voltage Peak',
                                                                             'DEPEND_0': 'Epoch',
                                                                             'UNITS': 'V',
                                                                             'VAR_TYPE': 'data'}]}
                                    }

            elif 'skinu' in file:
                data_dict_output = {**data_dict_output, **{f'Epoch_skinu': [Epoch, {'LABLAXIS': f'Epoch',
                                                                              'DEPEND_0': None,
                                                                              'UNITS': None,
                                                                              'VAR_TYPE': 'data'}],
                                                           f'skin_current': [temp_data[1],
                                                                            {'LABLAXIS': f'Skin Current',
                                                                             'DEPEND_0': 'Epoch',
                                                                             'UNITS': 'uA',
                                                                             'VAR_TYPE': 'data'}]}
                                    }

            elif 'temp' in file:
                data_dict_output = {**data_dict_output, **{f'Epoch_temp': [Epoch,{'LABLAXIS': f'Epoch',
                                                                             'DEPEND_0': None,
                                                                            'UNITS': None,
                                                                            'VAR_TYPE': 'data'}],
                                                           f'temperature': [temp_data[1],
                                                                                   {'LABLAXIS': f'Electron Temperature',
                                                                                    'DEPEND_0': 'Epoch',
                                                                                    'UNITS': 'eV',
                                                                                    'VAR_TYPE': 'data'}],
                                                           }
                                    }




    # --- Get everything on a SINGLE timebase ---
    # Description: The Skin current starts at 50s for every DERPA and the high/low flyers.
    # Choose this Epoch as the referene Epoch and interpolate the rest

    # Reference Epoch
    data_dict_output = {**data_dict_output, **{'Epoch': [np.array(deepcopy(data_dict_output['Epoch_skinu'][0])), deepcopy(data_dict_output['Epoch_skinu'][1])]}}
    from scipy.interpolate import CubicSpline
    seconds_current = stl.EpochTo_T0_Rocket(InputEpoch=data_dict_output['Epoch_skinu'][0],T0=dt.datetime(2022,11,20,17,20))
    bad_idx = np.abs(data_dict_output['Epoch'][0] - (T0 + dt.timedelta(seconds=95))).argmin()


    # get new temperature
    seconds_temp = stl.EpochTo_T0_Rocket(InputEpoch=data_dict_output['Epoch_temp'][0], T0=dt.datetime(2022, 11, 20, 17, 20))
    cs = CubicSpline(seconds_temp, data_dict_output['temperature'][0])
    data_dict_output['temperature'][0] = np.array(cs(seconds_current))

    data_dict_output['temperature'][0][:bad_idx] = 0

    # get new peak voltage
    seconds_voltage = stl.EpochTo_T0_Rocket(InputEpoch=data_dict_output['Epoch_peak'][0], T0=dt.datetime(2022, 11, 20, 17, 20))
    cs = CubicSpline(seconds_voltage, data_dict_output['peak_voltage'][0])
    data_dict_output['peak_voltage'][0] = np.array(cs(seconds_current))
    data_dict_output['peak_voltage'][0][:bad_idx] = 0

    # remove all the epoch keys
    data_dict_output.pop('Epoch_temp')
    data_dict_output.pop('Epoch_peak')
    data_dict_output.pop('Epoch_skinu')


    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        stl.prgMsg('Creating output file')
        fileoutName = f'ACESII_{rocketID}_l2_ERPA{wDERPA}.cdf'
        outputPath = f'C:\Data\ACESII\L2\{ACESII.fliers[wRocket-4]}' + f'\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
target_files = f'{DataPaths.ACES_data_folder}{input_path_modifier}{ACESII.fliers[wRocket-4]}\*.txt'

if len(glob(target_files)) == 0:
    print(stl.color.RED + 'There are no .txt files in the specified directory' + stl.color.END)
else:
    txt_to_l2_DERPA(wRocket-4, justPrintFileNames,wDERPA)

