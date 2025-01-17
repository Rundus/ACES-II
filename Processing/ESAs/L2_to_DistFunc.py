# --- L0_to_L1.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert electrostatic analyzer data from diffNFlux to Distribution Function



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import itertools
# --- --- --- --- ---
import time
from spaceToolsLib import Done,setupPYCDF,outputCDFdata,loadDictFromFile
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
wRocket = 4

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [[0, 2, 4], [0, 4]]

IsElectron = True
wIon = 0
inputPath_modifier = 'L2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'L3\DistFunc' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
from warnings import filterwarnings # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.
filterwarnings("ignore")
from tqdm import tqdm
from missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from class_var_func import L2_ACES_Quick,L2_TRICE_Quick
from spaceToolsLib import color,prgMsg,cm_to_m,q0,IonMasses,m_e
from glob import glob
from os.path import getsize
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)


def Distribution_Function(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    outputFolderPath = rocketFolderPath

    # --- ACES II Flight/Integration Data ---
    rocketAttrs,b,c = ACES_mission_dicts()
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'Distribution_Function'
    L2ModelData = L2_ACES_Quick(wflyer)

    # Set the paths for the file names
    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')
    outputFiles = glob(f'{outputFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*.cdf')

    input_names = [file.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\\', '') for file in inputFiles]
    output_names = [file.replace(f'{outputFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '') for file in outputFiles]

    input_names_searchable = [file.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace('l2_', '').replace('_v00','').replace('__', '_') for file in input_names]
    output_names_searchable = [file.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace('distFunc_', '').replace('_v00', '') for file in output_names]

    dataFile_name = input_names[wFile].replace(f'{rocketFolderPath}\{fliers[wflyer]}\\', '')
    fileoutName = dataFile_name.replace('l2', 'distFunc')

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in dataFile_name:
            wInstr = [index, instr]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:70s}{:5.1f} MB   Made ESACurrents: {:3s} '.format(i, input_names_searchable[i],round(getsize(file) / (10 ** 6), 1), anws[0]))
    else:
        print('\n')
        print(color.UNDERLINE + f'Calculating Distribution Function for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the tmCDF file ---
        prgMsg('Loading data from L2Files')
        data_dict = loadDictFromFile(inputFilePath=inputFiles[wFile])

        Done(start_time)

        # --- --- --- --- --- --- --- --- ---
        # --- Calculate Instrument Data ---
        # --- --- --- --- --- --- --- --- ---

        prgMsg('Calculating the Distribution Function')

        # --- CALCULATE DISTRIBUTION FUNCTION ---
        diffNFlux = data_dict['Differential_Number_Flux'][0]
        oneCountLevel = data_dict['oneCountLevel'][0]
        pitchAngle = data_dict['Pitch_Angle'][0]
        Energies = data_dict['Energy'][0]

        # define empty numpy array
        sizes = [len(diffNFlux),len(diffNFlux[0]), len(diffNFlux[0][0])]
        ranges = [range(sizes[0]), range(sizes[1]), range(sizes[2])]
        distFunc = np.zeros(shape=(sizes[0], sizes[1], sizes[2]))
        distFunc_oneCount = np.zeros(shape=(sizes[0], sizes[1], sizes[2]))

        if IsElectron:
            m = m_e
        else:
            m = IonMasses[wIon]

        print(f'\nNum. of iterations: {sizes[0]*sizes[1]*sizes[2]}\n')

        # --- Calculate DistFunc in SI units ---
        for tme, ptch, engy in tqdm(itertools.product(*ranges)):
            if diffNFlux[tme][ptch][engy] <= rocketAttrs.epoch_fillVal:
                distFunc[tme][ptch][engy] = rocketAttrs.epoch_fillVal
            else:
                distVal = (cm_to_m*cm_to_m/(q0*q0 ))*(((m**2)*diffNFlux[tme][ptch][engy]) / (2 * Energies[engy]))
                if distVal < 0:
                    distFunc[tme][ptch][engy] = 0
                else:
                    distFunc[tme][ptch][engy] = distVal

            distFunc_oneCount[tme][ptch][engy] = (cm_to_m*cm_to_m/(q0*q0 ))*(((m**2)*oneCountLevel[tme][ptch][engy]) / (2 * Energies[engy]))

        distFunc = np.array(distFunc)

        del data_dict['Differential_Number_Flux'],data_dict['Differential_Energy_Flux']

        Done(start_time)


        if outputData:
            # --- --- --- --- --- --- ---
            # --- WRITE OUT THE DATA ---
            # --- --- --- --- --- --- ---
            prgMsg('Creating output file')

            outputPath = f'{outputFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            data_dict = {**data_dict, **{'Distribution_Function':
                                             [distFunc, {'LABLAXIS': 'Distribution_Function',
                                                       'DEPEND_0': 'Epoch',
                                                       'DEPEND_1': 'Pitch_Angle',
                                                       'DEPEND_2': 'Energy',
                                                       'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                       'UNITS': 'm!A-6!Ns!A3!N',
                                                       'VALIDMIN': distFunc.min(), 'VALIDMAX': distFunc.max(),
                                                       'VAR_TYPE': 'data', 'SCALETYP': 'log'}]}}

            data_dict = {**data_dict, **{'oneCountLevel':
                                             [oneCountLevel, {'LABLAXIS': 'Distribution_Function',
                                                         'DEPEND_0': 'Epoch',
                                                         'DEPEND_1': 'Pitch_Angle',
                                                         'DEPEND_2': 'Energy',
                                                         'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                         'UNITS': 'm!A-6!Ns!A3!N',
                                                         'VALIDMIN': distFunc.min(), 'VALIDMAX': distFunc.max(),
                                                         'VAR_TYPE': 'support_data', 'SCALETYP': 'log'}]}}

            outputCDFdata(outputPath=outputPath,data_dict=data_dict,ModelData=L2ModelData,globalAttrsMod=globalAttrsMod,instrNam=wInstr[1])

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

if len(glob(f'{rocketFolderPath}L2\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        Distribution_Function(wRocket, 0, rocketFolderPath, justPrintFileNames, wflyer)
    elif not wFiles[wRocket-4]:
        for fileNo in (range(len(glob(f'{rocketFolderPath}L2\{fliers[wflyer]}\*.cdf')))):
            Distribution_Function(wRocket, fileNo, rocketFolderPath, justPrintFileNames, wflyer)
    else:
        for filesNo in wFiles[wRocket-4]:
            Distribution_Function(wRocket, filesNo, rocketFolderPath, justPrintFileNames, wflyer)