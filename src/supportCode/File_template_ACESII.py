# --- File_template_ACESII.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False
wRocket = 4
modifier = ''
inputPath_modifier = '' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = '' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# ---------------------------
outputData = False
# ---------------------------

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none


def main(wRocket, rocketFolderPath, justPrintFileNames):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wRocket-4]
    globalAttrsMod = rocketAttrs.globalAttributes[wRocket-4]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    ModelData = L2_TRICE_Quick(wRocket-4)

    inputFiles, input_names, input_names_searchable = getInputFiles(rocketFolderPath=rocketFolderPath,wRocket=wRocket,inputPath_modifier=inputPath_modifier)

    fileoutName = f'ACESII_{rocketID}_{wInstr[1]}_InvertedV_removed.cdf'

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
    # load data here
    Done(start_time)

    exampleAttrs = {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal,
                    'FORMAT': 'E12.2',
                    'UNITS': None, 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data', 'SCALETYP': 'linear',
                    'LABLAXIS': None}





    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        prgMsg('Creating output file')

        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wRocket-4]}\\{fileoutName}'

        outputCDFdata(outputPath, data_dict, ModelData, globalAttrsMod, wInstr[1])

        Done(start_time)





# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
elif wRocket == 5: # ACES II Low
    rocketFolderPath = ACES_data_folder



if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    main(wRocket, rocketFolderPath, justPrintFileNames)
