# --- template.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

from src.my_imports import *

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
wFiles = [2]


modifier = ''
inputPath_modifier = 'mag'
outputPath_modifier = 'l3'



# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from mat73 import loadmat
from copy import deepcopy

def mat_to_cdf(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    inputrocketFolderPath = rocketFolderPath + 'ACESII_matlab'
    outputrocketFolderPath = rocketFolderPath

    # --- ACES II Flight/Integration Data ---
    rocketID = rocketAttrs.rocketID[wflyer]

    inputFiles = glob(f'{inputrocketFolderPath}\{inputPath_modifier}\{ACESII.fliers[wRocket-4]}{modifier}\*.mat')
    outputFiles = glob(f'{outputrocketFolderPath}{outputPath_modifier}\{ACESII.fliers[wRocket-4]}\*.cdf')

    input_names = [ifile.replace(f'{inputrocketFolderPath}\{inputPath_modifier}\{ACESII.fliers[wRocket-4]}{modifier}\\', '') for ifile in inputFiles]
    output_names = [ofile.replace(f'{outputrocketFolderPath}{outputPath_modifier}\{ACESII.fliers[wRocket-4]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    output_names_searchable = [ofile.replace('ACES_', '').replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace(outputPath_modifier.lower() +'_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = inputFiles[wFile].replace(f'{inputrocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}{modifier}\\', '')

    fileoutName = dataFile_name.replace(f'{inputrocketFolderPath}\\{inputPath_modifier}\{ACESII.fliers[wRocket-4]}{modifier}\\', "").replace('.mat','.cdf').replace('l1','l2')


    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB   Made cdf: {:3s} '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
        return

    print('\n')
    print(stl.color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + stl.color.END)
    print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

    ##########################
    # --- DEFINE VARIABLES ---
    ##########################
    stl.prgMsg('Converting variables to .cdf format')
    exampleVar = {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': ACESII.epoch_fillVal,
                    'FORMAT': 'I5', 'UNITS': None, 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data',
                    'SCALETYP': 'linear', 'LABLAXIS': None}
    exampleEpoch = {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': ACESII.epoch_fillVal,
              'FORMAT': 'I5', 'UNITS': 'ns', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data',
              'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000', 'TIME_SCALE': 'Terrestrial Time',
              'REFERENCE_POSITION': 'Rotating Earth Geoid', 'SCALETYP': 'linear', 'LABLAXIS': 'Epoch'}

    # DETERMINE THE FILETYPE
    loadThisFile = inputFiles[wFile]

    try:
        # Load in data to data_dict
        mat = loadmat(loadThisFile)

    except Exception as e:
        print('\n')
        print(e,'\n')
        stl.prgMsg('Loading with scipy.io instead')
        mat = scipy.io.loadmat(loadThisFile)

    # --- convert data into dictionary ---
    data_dict = {}
    for key,val in mat.items():
        if key not in ['__header__','__version__','__globals__','__function_workspace__','None']:
            if key.lower() in ['epoch','Time']:
                data_dict = {**data_dict,**{'epoch':[np.array(val),deepcopy(exampleEpoch)]}}
            else:
                data_dict = {**data_dict, **{key: [np.array(val), deepcopy(exampleVar)]}}


    # --- Handle the Time variable ---
    # If "Time" variable is actually "Time since launch" do the conversion from SECONDS
    # if data_dict['Epoch'][0][0] < 100000000000000000:
    #     data_dict['Epoch'][0] = np.array([data_dict['Epoch'][0][i]*(10**(9)) + rocketAttrs.Launch_Times[wRocket - 4] for i in range(len(data_dict['Epoch'][0]))])


    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    stl.prgMsg('Creating output file')

    outputPath = f'{outputrocketFolderPath}{outputPath_modifier}\\{fileoutName}'

    # --- --- --- --- --- ---
    # --- WRITE OUT DATA ---
    # --- --- --- --- --- ---
    stl.outputCDFdata(outputPath=outputPath,data_dict=data_dict)

    stl.Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
mat_to_cdf(wRocket)