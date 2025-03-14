# --- template.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
from src.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False # Just print the names of files

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
wFiles = [0]

# --- OutputData ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src. import XXX as fToggles

#######################
# --- MAIN FUNCTION ---
#######################
def main(wflyer, wFile, justPrintFileNames):

    # --- FILE I/O ---
    rocket_folder_path = DataPaths.ACES_data_folder
    rocketID = ACESII.payload_IDs[wflyer]
    globalAttrsMod = deepcopy(ACESII.global_attributes[wflyer])
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'LangmuirData'

    # Set the paths for the file names
    data_repository = f'{rocket_folder_path}{fToggles.inputPath_modifier}\\{ACESII.fliers[wflyer]}{fToggles.modifier}\\'
    input_files = glob(data_repository+'*langmuir*')
    input_names = [ifile.replace(data_repository, '') for ifile in input_files]
    input_names_searchable = [ifile.replace(fToggles.inputPath_modifier.lower() + '_', '').replace('_v00', '') for ifile in input_names]
    dataFile_name = input_files[wFile].replace(data_repository,'')

    if justPrintFileNames:
        for i, file in enumerate(input_files):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(os.path.getsize(file) / (10 ** 6), 1)))
        return

    print('\n')
    print(stl.color.UNDERLINE + f'Converting to {fToggles.outputPath_modifier} data for {dataFile_name}' + stl.color.END)
    print('[' + str(wFile) + ']   ' + str(round(os.path.getsize(input_files[wFile]) / (10 ** 6), 1)) + 'MiB')

    # --- get the data from the L2 file ---
    stl.prgMsg(f'Loading data from {fToggles.inputPath_modifier} Files')
    data_dict = stl.loadDictFromFile(input_files[wFile])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}


    ##################################################################
    # --- Calculate the plasma density from Ion Saturation Current ---
    ##################################################################
    stl.prgMsg('Calculating Fixed ni')


    varAttrs = {'LABLAXIS': 'plasma density', 'DEPEND_0': 'Epoch',
                                                                   'DEPEND_1': None,
                                                                   'DEPEND_2': None,
                                                                   'FILLVAL': ACESII.epoch_fillVal,
                                                                   'FORMAT': 'E12.2',
                                                                   'UNITS': '!Ncm!A-3!N',
                                                                   'VALIDMIN': 0,
                                                                   'VALIDMAX': 0,
                                                                   'VAR_TYPE': 'data', 'SCALETYP': 'linear'}




    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        stl.prgMsg('Creating output file')
        fileoutName_fixed = f'ACESII_{rocketID}_langmuir_fixed.cdf'
        outputPath = f'{rocket_folder_path}{fToggles.outputPath_modifier}\{ACESII.fliers[wflyer]}\\{fileoutName_fixed}'
        stl.outputCDFdata(outputPath, data_dict_fixed, instrNam= 'Langmuir Probe',globalAttrsMod=globalAttrsMod)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
target_files = f'{DataPaths.ACES_data_folder}{fToggles.inputPath_modifier}\{ACESII.fliers[wRocket-4]}{fToggles.modifier}\*.cdf'

if len(glob(target_files)) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    for file_idx in wFiles:
        main(wRocket-4, file_idx, justPrintFileNames)

