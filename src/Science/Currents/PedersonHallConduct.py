# --- template.py ---
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
wFiles = []

modifier = ''
inputPath_modifier_attitude = 'attitude' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_density = 'science\Langmuir' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_mag = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\PedersenHallConduct' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

outputData = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none

def PedersonHallConduct(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)


    inputFiles_attitude = glob(f'{rocketFolderPath}{inputPath_modifier_attitude}\{fliers[wflyer]}{modifier}\*.cdf')
    inputFiles_density = glob(f'{rocketFolderPath}{inputPath_modifier_density}\{fliers[wflyer]}{modifier}\*.cdf')
    inputFiles_mag = glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wflyer]}{modifier}\*RingCore_rktFrm*')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier_attitude}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles_attitude]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier_attitude.lower() +'_', '').replace('_v00', '') for ifile in input_names]

    fileoutName = f'ACESII_{rocketID}_{outputPath_modifier.lower()}_PedersenHallConduct.cdf'


    if justPrintFileNames:
        for i, file in enumerate(inputFiles_attitude):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Calculating Hall and Pederson Conductivities' + color.END)

        # --- get the data from the attitude file ---
        prgMsg(f'Loading data from Attitude Files')
        data_dict_attitude = loadDictFromFile(inputFiles_attitude[0],{})
        data_dict_attitude['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_attitude['Epoch'][0][i]) for i in (range(len(data_dict_attitude['Epoch'][0])))])
        Done(start_time)

        # --- get the data from the density file ---
        prgMsg(f'Loading data from Density Files')
        data_dict_density = loadDictFromFile(inputFiles_density[0], {})
        data_dict_density['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_density['Epoch'][0][i]) for i in (range(len(data_dict_density['Epoch'][0])))])
        Done(start_time)

        # --- get the data from the integration_tad_files file ---
        prgMsg(f'Loading data from Magnetometer Files')
        data_dict_mag = loadDictFromFile(inputFiles_mag[0], {})
        data_dict_mag['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_mag['Epoch'][0][i]) for i in (range(len(data_dict_mag['Epoch'][0])))])
        Done(start_time)

        #######################
        # --- DO STUFF HERE ---
        #######################








        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('Creating output file')

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            outputCDFdata(outputPath, data_dict, ModelData, globalAttrsMod,wInstr[1])

            Done(start_time)







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
        PedersonHallConduct(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            PedersonHallConduct(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            PedersonHallConduct(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)