# --- File_template_ACESII.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:



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
justPrintFileNames = False
wRocket = 5
modifier = ''
inputPath_modifier = 'L2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\ExB_drift' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
BFile = 7
EFile = 2

# ---------------------------
outputData = True
# ---------------------------

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none


def ExB_drift(wRocket, rocketFolderPath, justPrintFileNames):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wRocket-4]
    globalAttrsMod = rocketAttrs.globalAttributes[wRocket-4]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    ModelData = L2_TRICE_Quick(wRocket-4)
    inputFiles, input_names, input_names_searchable = getInputFiles(rocketFolderPath=rocketFolderPath, wRocket=wRocket, inputPath_modifier=inputPath_modifier)
    fileoutName = f'ACESII_{rocketID}_ExB_drift.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return

    print('\n')


    # --- --- --- --- --- -
    # --- LOAD THE DATA ---
    # --- --- --- --- --- -

    # --- get the data from the file ---
    prgMsg(f'Loading data from {inputPath_modifier} Files')
    data_dict_B, globAttrs = loadDictFromFile(inputFilePath=inputFiles[BFile], getGlobalAttrs=True)
    data_dict_E = loadDictFromFile(inputFilePath= inputFiles[EFile])
    Done(start_time)

    # --- get data into vector form ---
    coordcompNames_B, coordSetName, coordSet = getCoordinateKeys(data_dict_B)
    coordcompNames_E, coordSetName, coordSet = getCoordinateKeys(data_dict_E)

    Evec = np.array([data_dict_E[coordcompNames_E[0]][0], data_dict_E[coordcompNames_E[1]][0], data_dict_E[coordcompNames_E[2]][0]]).T
    Bvec = (1E-9)*np.array([data_dict_B[coordcompNames_B[0]][0], data_dict_B[coordcompNames_B[1]][0], data_dict_B[coordcompNames_B[2]][0]]).T
    Bmag = (1E-9)*data_dict_B['Bmag'][0]

    # --- --- --- --- --- -
    # --- CALCULATE ExB ---
    # --- --- --- --- --- -
    ExB = (1/1000)*np.cross(Evec, Bvec)
    ExBvec = np.array([ExB[i]/(Bmag[i]**2)for i in range(len(Bmag))])

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        prgMsg('Creating output file')
        exampleAttrs = {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal,
                        'FORMAT': 'E12.2',
                        'UNITS': None, 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data', 'SCALETYP': 'linear',
                        'LABLAXIS': None}
        data_dict_output = {'ExB_e':[ExBvec[:, 0], exampleAttrs],
                            'ExB_p': [ExBvec[:, 1], exampleAttrs],
                            'ExB_r': [ExBvec[:, 2], exampleAttrs],
                            'Epoch':deepcopy(data_dict_B['Epoch']),
                            'Alt':deepcopy(data_dict_B['Alt']),
                            'ILat':deepcopy(data_dict_B['ILat']),
                            'ILong':deepcopy(data_dict_B['ILat'])}
        data_dict_output['ExB_e'][1]['UNITS'] = 'km/s'
        data_dict_output['ExB_p'][1]['UNITS'] = 'km/s'
        data_dict_output['ExB_r'][1]['UNITS'] = 'km/s'

        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wRocket-4]}\\{fileoutName}'

        outputCDFdata(outputPath, data_dict_output, globalAttrsMod=globAttrs, instrNam='ExB')

        Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    ExB_drift(wRocket, rocketFolderPath, justPrintFileNames)
