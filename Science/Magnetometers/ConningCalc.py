# --- template.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: determine the coning angle for each payload
# throughout flight



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import itertools
# --- --- --- --- ---

import time
from ACESII_code.class_var_func import Done, setupPYCDF

start_time = time.time()
# --- --- --- --- ---



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket =5

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]
wFile_mag = 1

modifier = ''
inputPath_modifier = 'attitude' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_mag = 'mag' # e.g. 'L1' or 'L1'. It's the name of the broader input folde
outputPath_modifier = 'attitude' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import os
from tqdm import tqdm
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, prgMsg, L2_TRICE_Quick, outputCDFdata
from glob import glob
from os.path import getsize

setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)


def main(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)


    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf')
    inputFiles_mag = glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wflyer]}{modifier}\*Ringcore*')
    outputFiles = glob(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    input_names_mag = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles_mag]
    output_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    input_names_searchable_mag = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier_mag.lower() + '_', '').replace('_v00', '') for ifile in input_names_mag]
    output_names_searchable = [ofile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(outputPath_modifier.lower() +'_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '')
    dataFile_name_mag = inputFiles_mag[wFile_mag].replace(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wflyer]}{modifier}\\', '')
    fileoutName = dataFile_name.replace(inputPath_modifier.lower(), outputPath_modifier.lower())

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in dataFile_name:
            descriptiorNam = ['EEPAA', 'LEESA', 'IEPAA', 'Langmuir_Probe']
            wInstr = [index, instr, descriptiorNam[index]]
        else:
            wInstr = [0, 'attitude', '']

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB   Made L2: {:3s} '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the ACS file ---
        prgMsg(f'Loading data from {dataFile_name}')
        data_dict = {}
        with pycdf.CDF(inputFiles[wFile]) as inputDataFile:
            for key, val in inputDataFile.items():
                data_dict = {**data_dict, **{key : [inputDataFile[key][...], {key:val for key,val in inputDataFile[key].attrs.items()  }  ]  }  }

        data_dict['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch'][0][i]) for i in (range(len(data_dict['Epoch'][0])))])

        Done(start_time)

        # --- get the data from the mag file ---
        prgMsg(f'Loading data from {inputFiles_mag[wFile_mag]}')
        data_dict_mag = {}
        with pycdf.CDF(inputFiles_mag[wFile_mag]) as inputDataFile:
            for key, val in inputDataFile.items():
                data_dict_mag = {**data_dict_mag, **{key: [inputDataFile[key][...], {key: val for key, val in inputDataFile[key].attrs.items()}]}}

        data_dict_mag['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_mag['Epoch'][0][i]) for i in (range(len(data_dict_mag['Epoch'][0])))])

        Done(start_time)

        ######################################
        # --- Determine Magnetic Alignment ---
        ######################################

        # Determine the angle between the payload spin axis and the ambient magnetic field using the REAL magnetometer data

        # Since the mag_payload data is already in payload coordniantes, we just need the angle between B_mag and the spin axis (X-axis in payload frame) or +Yaxis (in RingCore frame)
        spin_axis = np.array([1,0,0])

        # create normalized B_components array
        Bcomps = np.array(
            [
                [data_dict_mag['Bx'][0][i],data_dict_mag['By'][0][i],data_dict_mag['Bz'][0][i]]
                for i in range(len(data_dict_mag['Bx'][0]))
            ]
        )

        b_norm = np.array([v/np.linalg.norm(v) for v in Bcomps])

        spinAxis_B_alignment = np.array(
            [180 - np.degrees(np.arccos(np.dot(spin_axis,vec))) for vec in b_norm]
        )

        # add variable to data dict
        data_dict = {**data_dict, **{'Epoch_B':
                                         [data_dict_mag['Epoch'][0], {'LABLAXIS': 'Epoch_B',
                                                            'DEPEND_0': 'Epoch_B', 'DEPEND_1': None,
                                                            'DEPEND_2': None,
                                                            'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                            'UNITS': 'ns',
                                                            'VALIDMIN': data_dict_mag['Epoch'][0].min(),
                                                            'VALIDMAX': data_dict_mag['Epoch'][0].max(),
                                                            'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]}}

        data_dict = {**data_dict, **{'spinAxis_B_alignment':
                                         [spinAxis_B_alignment, {'LABLAXIS': 'spinAxis_B_alignment',
                                                         'DEPEND_0': 'Epoch_B', 'DEPEND_1': None,
                                                         'DEPEND_2': None,
                                                         'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                         'UNITS': 'deg',
                                                         'VALIDMIN': spinAxis_B_alignment.min(), 'VALIDMAX': spinAxis_B_alignment.max(),
                                                         'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('Creating output file')

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            outputCDFdata(outputPath, data_dict, outputModelData,globalAttrsMod,wInstr[2])

        Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 0:  # ACES II Integration High
    rocketFolderPath = Integration_data_folder
    wflyer = 0
elif wRocket == 1: # ACES II Integration Low
    rocketFolderPath = Integration_data_folder
    wflyer = 1
elif wRocket == 2:  # TRICE II High
    rocketFolderPath = TRICE_data_folder
    wflyer = 0
elif wRocket == 3: # TRICE II Low
    rocketFolderPath = TRICE_data_folder
    wflyer = 1
elif wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        main(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            main(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            main(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)