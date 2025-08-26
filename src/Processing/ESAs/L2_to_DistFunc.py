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
from src.my_imports import *

def L2_to_DistFunc(wRocket, wFile, rocketFolderPath, justPrintFileNames):

    outputFolderPath = rocketFolderPath

    # --- ACES II Flight/Integration Data ---
    # Set the paths for the file names
    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}\*.cdf')
    outputFiles = glob(f'{outputFolderPath}{outputPath_modifier}\{ACESII.fliers[wRocket-4]}\*.cdf')

    input_names = [file.replace(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}\\', '') for file in inputFiles]
    output_names = [file.replace(f'{outputFolderPath}{outputPath_modifier}\{ACESII.fliers[wRocket-4]}\\', '') for file in outputFiles]

    input_names_searchable = [file.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace('l2_', '').replace('_v00','').replace('__', '_') for file in input_names]
    output_names_searchable = [file.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace('distFunc_', '').replace('_v00', '') for file in output_names]

    dataFile_name = input_names[wFile].replace(f'{rocketFolderPath}\{ACESII.fliers[wRocket-4]}\\', '')
    fileoutName = dataFile_name.replace('l2', 'distFunc')

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa']):
        if instr in dataFile_name:
            wInstr = [index, instr]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:70s}{:5.1f} MB   Made ESACurrents: {:3s} '.format(i, input_names_searchable[i],round(getsize(file) / (10 ** 6), 1), anws[0]))
        return

    print('\n')
    print(stl.color.UNDERLINE + f'Calculating Distribution Function for {dataFile_name}' + stl.color.END)
    print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

    # --- get the data ---
    stl.prgMsg('Loading data from L2Files')
    data_dict_esa = stl.loadDictFromFile(inputFilePath=inputFiles[wFile])
    data_dict_meff_i = stl.loadDictFromFile(rf'C:\Data\ACESII\science\ion_mass_effective\\{ACESII.fliers[wRocket-4]}\\*.cdf*')
    Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}

    # --- --- --- --- --- --- --- --- ---
    # --- Calculate Instrument Data ---
    # --- --- --- --- --- --- --- --- ---
    stl.prgMsg('Calculating the Distribution Function')

    # --- CALCULATE DISTRIBUTION FUNCTION ---
    diffNFlux = data_dict_esa['Differential_Number_Flux'][0]
    oneCountLevel = data_dict_esa['oneCountLevel'][0]
    Energies = data_dict_esa['Energy'][0]

    # define empty numpy array
    sizes = [len(diffNFlux),len(diffNFlux[0]), len(diffNFlux[0][0])]
    ranges = [range(sizes[0]), range(sizes[1]), range(sizes[2])]
    distFunc = np.zeros(shape=(sizes[0], sizes[1], sizes[2]))
    distFunc_oneCount = np.zeros(shape=(sizes[0], sizes[1], sizes[2]))


    print(f'\nNum. of iterations: {sizes[0]*sizes[1]*sizes[2]}\n')

    # --- Calculate DistFunc in SI units ---
    for tme, ptch, engy in tqdm(itertools.product(*ranges)):
        if diffNFlux[tme][ptch][engy] <= ACESII.epoch_fillVal:
            distFunc[tme][ptch][engy] = ACESII.epoch_fillVal
        else:
            distVal = (stl.cm_to_m*stl.cm_to_m/(stl.q0*stl.q0 ))*(((stl.m**2)*diffNFlux[tme][ptch][engy]) / (2 * Energies[engy]))
            if distVal < 0:
                distFunc[tme][ptch][engy] = 0
            else:
                distFunc[tme][ptch][engy] = distVal

        distFunc_oneCount[tme][ptch][engy] = (stl.cm_to_m*stl.cm_to_m/(stl.q0*stl.q0 ))*(((stl.m**2)*oneCountLevel[tme][ptch][engy]) / (2 * Energies[engy]))

    Done(start_time)




    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        stl.prgMsg('Creating output file')
        outputPath = f'{outputFolderPath}{outputPath_modifier}\{ACESII.fliers[wRocket-4]}\\{fileoutName}'

        data_dict_output = {**data_dict_output, **{'Distribution_Function':
                                         [distFunc, {'LABLAXIS': 'Distribution_Function',
                                                   'DEPEND_0': 'Epoch',
                                                   'DEPEND_1': 'Pitch_Angle',
                                                   'DEPEND_2': 'Energy',
                                                   'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                   'UNITS': 'm!A-6!Ns!A3!N',
                                                   'VALIDMIN': distFunc.min(), 'VALIDMAX': distFunc.max(),
                                                   'VAR_TYPE': 'data', 'SCALETYP': 'log'}]}}

        data_dict = {**data_dict_output, **{'oneCountLevel':
                                         [oneCountLevel, {'LABLAXIS': 'Distribution_Function',
                                                     'DEPEND_0': 'Epoch',
                                                     'DEPEND_1': 'Pitch_Angle',
                                                     'DEPEND_2': 'Energy',
                                                     'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                     'UNITS': 'm!A-6!Ns!A3!N',
                                                     'VALIDMIN': distFunc.min(), 'VALIDMAX': distFunc.max(),
                                                     'VAR_TYPE': 'support_data', 'SCALETYP': 'log'}]}}

        outputCDFdata(outputPath=outputPath,data_dict=data_dict)

        Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

rocketFolderPath = DataPaths.ACES_data_folder


if len(glob(f'{rocketFolderPath}L2\{ACESII.fliers[wRocket-4]}\*.cdf')) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    if justPrintFileNames:
        L2_to_DistFunc(wRocket, 0, rocketFolderPath, justPrintFileNames)
    elif not wFiles[wRocket-4]:
        for fileNo in (range(len(glob(f'{rocketFolderPath}L2\{ACESII.fliers[wRocket-4]}\*.cdf')))):
            L2_to_DistFunc(wRocket, fileNo, rocketFolderPath, justPrintFileNames)
    else:
        for filesNo in wFiles[wRocket-4]:
            L2_to_DistFunc(wRocket, filesNo, rocketFolderPath, justPrintFileNames)