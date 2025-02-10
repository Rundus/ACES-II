# --- magFiltering.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: I use Autoplot to create nice spectrograms of the magnetometer data
# I use this file to clean up/filter the mag data

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
wFiles = [0] # try to use _Spun.cdf data

modifier = ''
inputPath_modifier = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modelIGRF_modifier = 'science/IGRF_interpolated'
outputPath_modifier = 'science/magSpectrograms' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# if using ENU coordniates
useENUcoordinates = False

# additional toggel
subtractModelIGRF = False

##################
# FILTER TOGGLES #
##################
filtOrder = 8 # order of filter
low_cutoff_Freq = 2 # cut off frequency where gain has reached -3dB
high_cutoff_Freq = 0
dataSampleFreq = 128 # sample per second of the data
filterType = 'highpass'
plotFilteredData = True

# output the data
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.class_var_func import butter_filter


def magSpectrograms(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'Science'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*RingCore_*')
    inputFiles_modelIGRF = glob(f'{rocketFolderPath}{inputPath_modelIGRF_modifier}\{fliers[wflyer]}{modifier}\*IGRF*')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]

    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '')

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['Tesseract', 'RingCore']):
        if instr.lower() in dataFile_name.lower():
            wInstr = [index, instr, instr]

    fileoutName = f'ACESII_{rocketID}_{wInstr[1]}_Bandpass.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the tmCDF file ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict = loadDictFromFile(inputFiles[wFile], {})
        data_dict['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch'][0][i]) for i in (range(len(data_dict['Epoch'][0])))])
        Done(start_time)

        #########################
        # --- FITLER THE DATA ---
        #########################

        if useENUcoordinates:
            magAxes = ['B_east', 'B_north', 'B_up']
        else:
            magAxes = ['Bx', 'By', 'Bz']


        if subtractModelIGRF:
            prgMsg('Subtracting off IGRF Field')
            data_dict_IGRF = loadDictFromFile(inputFiles_modelIGRF[0],{})
            IGRFAxes= ['IGRF_East','IGRF_North','IGRF_Up']

            for j,axes in enumerate(magAxes):
                data_dict[axes][0] = np.array([data_dict[axes][0][i] - data_dict_IGRF[IGRFAxes[j]][0][i] for i in range(len(data_dict_IGRF['Epoch'][0]))])

            Done(start_time)

        # remove fillvals
        badIndicies = []

        for i in range(len(data_dict['Epoch'][0])):
            if data_dict['Epoch'][0][i] <= 0:
                badIndicies.append(i)
            elif np.isnan(data_dict[magAxes[0]][0][i]):
                badIndicies.append(i)
            elif np.isnan(data_dict[magAxes[1]][0][i]):
                badIndicies.append(i)
            elif np.isnan(data_dict[magAxes[2]][0][i]):
                badIndicies.append(i)

        data_dict['Epoch'][0] = np.delete(data_dict['Epoch'][0], badIndicies)
        data_dict[magAxes[0]][0] = np.delete(data_dict[magAxes[0]][0], badIndicies)
        data_dict[magAxes[1]][0] = np.delete(data_dict[magAxes[1]][0], badIndicies)
        data_dict[magAxes[2]][0] = np.delete(data_dict[magAxes[2]][0], badIndicies)
        rawData = np.array([data_dict[axes][0] for axes in magAxes])


        FilteredData = np.array([
            butter_filter(data=data_dict[axes][0],
                          lowcutoff=low_cutoff_Freq,
                          highcutoff=high_cutoff_Freq,
                          fs=dataSampleFreq,
                          order=filtOrder,
                          filtertype=filterType)
            for axes in magAxes])

        time = np.array([pycdf.lib.tt2000_to_datetime(tme) for tme in data_dict['Epoch'][0]])

        if plotFilteredData:
            targetTimeLower = pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 24, 55, 00))
            indexLower = np.abs(data_dict['Epoch'][0] - targetTimeLower).argmin()
            targetTimeUpper = pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 25, 9, 00))
            indexUpper = np.abs(data_dict['Epoch'][0] - targetTimeUpper).argmin()

            Epoch = data_dict['Epoch'][0][indexLower:indexUpper]
            Epoch_reduced = np.array([pycdf.lib.tt2000_to_datetime(tme) for tme in Epoch])
            B1 = data_dict[magAxes[0]][0][indexLower:indexUpper]
            B2 = data_dict[magAxes[1]][0][indexLower:indexUpper]
            B3 = data_dict[magAxes[2]][0][indexLower:indexUpper]
            B = [B1, B2, B3]
            B_filtered = [thing[indexLower:indexUpper] for thing in FilteredData]

            wComp = 1
            fig, ax = plt.subplots(2)
            fig.suptitle(f'HighCutOff: {high_cutoff_Freq} Hz\n'
                         f'LowCutOFf: {low_cutoff_Freq} Hz\n'
                         f'Order: {filtOrder}')
            ax[0].plot(Epoch_reduced, B[wComp])
            ax[0].set_ylabel(f'{magAxes[wComp]}_raw')
            ax[1].plot(Epoch_reduced, B_filtered[wComp])
            ax[1].set_ylabel(f'{magAxes[wComp]}_filtered')
            ax[1].set_xlabel('Time')
            plt.show()

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            # Add filtered data to data dict
            bandpass_limits = np.array([low_cutoff_Freq, high_cutoff_Freq])
            data_dict = {**data_dict, **{f'Bandpass_Limits':
                                             [bandpass_limits, {'LABLAXIS': f'Bandpass_Limits',
                                                                    'DEPEND_0':'Bandpass_Limits',
                                                                    'DEPEND_1': None,
                                                                    'DEPEND_2': None,
                                                                    'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                                    'UNITS': 'nT',
                                                                    'VALIDMIN': bandpass_limits.min(),
                                                                    'VALIDMAX': bandpass_limits.max(),
                                                                    'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}

            for i, data in enumerate(FilteredData):
                data_dict = {**data_dict, **{f'{magAxes[i]}_filtered':
                                                 [data, {'LABLAXIS': f'{magAxes[i]}_filtered',
                                                                        'DEPEND_0': 'Epoch',
                                                                        'DEPEND_1': None,
                                                                        'DEPEND_2': None,
                                                                        'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                                        'UNITS': 'nT',
                                                                        'VALIDMIN': data.min(),
                                                                        'VALIDMAX': data.max(),
                                                                        'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
            # --- --- --- --- --- --- ---
            # --- WRITE OUT THE DATA ---
            # --- --- --- --- --- --- ---
            prgMsg('Creating output file')

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            globalAttrsMod['Descriptor'] = wInstr[1]
            outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod, wInstr[1])

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
        magSpectrograms(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            magSpectrograms(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            magSpectrograms(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)