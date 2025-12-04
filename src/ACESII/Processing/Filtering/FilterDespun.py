# --- FilterDespun.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Butterworth Filter the E and B data



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
wRocket = 5
modifier = ''
inputPath_modifier = 'L2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder inside data\ACESII
outputPath_modifier = 'L2' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder inside data\ACESII\ACESII_matlab
wFiles = [10]
# ----------------------------
lowcutoff, highcutoff = 0.7, 0.7
fs = 128
order = 4
# filtertype = 'HighPass'
filtertype = 'HighPass'
# -----------------------------
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.class_var_func import butter_filter


def FilterDespun(wRocket, wFile,rocketFolderPath, justPrintFileNames):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wRocket-4]

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}{modifier}\*.cdf')
    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}{modifier}\\', '') for ifile in inputFiles]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names[i], round(getsize(file) / (10 ** 6), 1)))
        return


    # --- --- --- --- --- -
    # --- LOAD THE DATA ---
    # --- --- --- --- --- -

    # --- get the data from the file ---
    prgMsg(f'Loading data from {inputPath_modifier} Files')
    data_dict,GlobalAttrs = loadDictFromFile(inputFiles[wFile],getGlobalAttrs=True)
    Done(start_time)

    # --- --- --- --- --- ---
    # --- FILTER THE DATA ---
    # --- --- --- --- --- ---

    # LP Probe Data
    if 'langmuir' in  input_names[wFile].lower():
        data = deepcopy(data_dict['ni'][0])
        data_dict['ni'][0] = butter_filter(data, lowcutoff, highcutoff, fs, order, filtertype)
        fileoutName = input_names[wFile].replace('.cdf', f'_{filtertype}_low{lowcutoff}_high{highcutoff}.cdf')

    else: # E-Field and B-Field data
        compNames, coordSys, coordSet = getCoordinateKeys(data_dict)
        for comp in compNames:
            data = deepcopy(data_dict[comp][0])
            data_dict[comp][0] = butter_filter(data, lowcutoff, highcutoff, fs, order, filtertype)

        fileoutName = input_names[wFile].replace(coordSys, f'{coordSys}_{filtertype}_low{lowcutoff}_high{highcutoff}')

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        prgMsg('Creating output file')

        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wRocket-4]}\\{fileoutName}'
        outputCDFdata(outputPath, data_dict, globalAttrsMod=GlobalAttrs)

        Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        FilterDespun(wRocket, 0, rocketFolderPath, justPrintFileNames)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}\*.cdf')))):
            FilterDespun(wRocket, fileNo, rocketFolderPath, justPrintFileNames)
    else:
        for filesNo in wFiles:
            FilterDespun(wRocket, filesNo, rocketFolderPath, justPrintFileNames)


