# --- groundTrackVelocty.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Code to estimate the groundtrack velocity - intentially doesn't
# account for any velocity in the altitude direction

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

from ACESII_code.myImports import *

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
wRocket = 5

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

modifier = ''
inputPath_modifier = 'attitude' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'attitude' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none

def groundTrackVelocity(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf')
    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    input_names_searchable = [ifile.replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    fileoutName = input_names[wFile]
    print(fileoutName)

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in fileoutName:
            descriptiorNam = ['EEPAA', 'LEESA', 'IEPAA', 'Langmuir_Probe']
            wInstr = [index, instr, descriptiorNam[index]]
        else:
            wInstr = [0, 'attitude', '']

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print(color.CYAN + 'Calculating Ground Track Velocity' + color.END)

        # --- get the data from the tmCDF file ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict = loadDictFromFile(inputFiles[wFile], {})
        data_dict['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch'][0][i]) for i in (range(len(data_dict['Epoch'][0])))])
        Done(start_time)

        #######################
        # --- DO STUFF HERE ---
        #######################
        from ACESII_code.class_var_func import lat_to_meter,long_to_meter
        lat_km = np.array([lat_to_meter*lat for lat in data_dict['Latgd'][0]])
        long_km = np.array([long_to_meter(data_dict['Long'][0][i], data_dict['Latgd'][0][i]) for i in range(len(data_dict['Long'][0]))])
        Epoch_dt = np.array([ pycdf.lib.tt2000_to_datetime(data_dict['Epoch'][0][i]) for i in range(1,len(data_dict['Epoch'][0]))])
        Epoch = np.array([tme/1E9 for tme in data_dict['Epoch'][0]])

        # --- estimate ground track ---
        groundtrackVel = []

        for i in range(len(Epoch)-1):
            groundtrackVel.append(
                np.sqrt(
                    (lat_km[i+1] - lat_km[i])**2 + (long_km[i+1] - long_km[i])**2
                )/ (Epoch[i+1] - Epoch[i])
            )

        # add one value to the end of the data so it matches data_dict['Epoch]
        groundtrackVel.append(0)
        groundtrackVel = np.array(groundtrackVel)

        data_dict = {**data_dict, **{'Ground_Track_Vel':
                                         [groundtrackVel, {'LABLAXIS': 'Ground Track Velocity',
                                                                 'DEPEND_0': 'Epoch', 'DEPEND_1': None,
                                                                 'DEPEND_2': None,
                                                                 'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                                 'UNITS': 'km/s',
                                                                 'VALIDMIN': groundtrackVel.min(),
                                                                 'VALIDMAX': groundtrackVel.max(),
                                                                 'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('Creating output file')

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod, wInstr[2])




# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

if wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        groundTrackVelocity(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            groundTrackVelocity(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            groundTrackVelocity(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)