# --- cdf_to_mat.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: script to convert .cdfs to .mat by looking into a directory or naming a specific file



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
inputPath_modifier = 'attitude' # e.g. 'L1' or 'L1'. It's the name of the broader input folder inside data\ACESII
outputPath_modifier = 'attitude' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder inside data\ACESII\ACESII_matlab


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.io import savemat


def cdf_to_mat(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):
    outputrocketFolderPath = rocketFolderPath + 'ACESII_matlab'

    # --- ACESII ---
    if wRocket in [0, 1, 4, 5]:
        # --- ACES II Flight/Integration Data ---
        rocketAttrs, b, c = ACES_mission_dicts()
        rocketID = rocketAttrs.rocketID[wflyer]
        globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
        globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
        outputModelData = L0_ACES_Quick(wflyer)

    # --- TRICE II ---
    elif wRocket in [2, 3]:
        globalAttrsMod = {}
        rocketAttrs, b, c = TRICE_mission_dicts()
        rocketID = rocketAttrs.rocketID[wflyer]
        outputModelData = L0_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf')
    outputFiles = glob(f'{outputrocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    output_names = [ofile.replace(f'{outputrocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    output_names_searchable = [ofile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(outputPath_modifier.lower() +'_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '')
    fileoutName = dataFile_name.replace(inputPath_modifier.lower(), outputPath_modifier.lower()).replace('.cdf','.mat')

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB   Made L2: {:3s} '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the tmCDF file ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict = loadDictFromFile(inputFiles[wFile])

        # convert epoch to tt2000
        data_dict['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch'][0][i]) for i in range(len(data_dict['Epoch'][0]))])
        # data_dict['Epoch_monitors'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_monitors'][0][i]) for i in range(len(data_dict['Epoch_monitors'][0]))])
        Done(start_time)

        #####################################
        # --- prepare dictionary for .mat ---
        #####################################
        # remove the attributes attached to the dictonary
        mdic = {}
        for key, val in data_dict.items():
            mdic = {**mdic, **{key:val[0]}}

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        prgMsg('Creating output file')

        outputPath = f'{outputrocketFolderPath}\{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'
        # --- delete output file if it already exists ---
        if os.path.exists(outputPath):
            os.remove(outputPath)

        # output the data
        savemat(outputPath, mdict=mdic)

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
        cdf_to_mat(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            cdf_to_mat(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            cdf_to_mat(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)