# --- L1_to_rebin_padPair_collection_cal1.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Get the "true" pitch angle that was calculated by the RingCore integration_tad_files and
# then use it to re-bin the ESA data and collect ChiSquare-Calibrations data
# reduce dataset to only include data after ACS despun rocket since this is where the magPitchAngles are nominal.
# Value is chosen based off of the attitude data variable "SpinAlignmentToB"

# The ChiSquare Calibration values are chosen only AFTER the
# ACS has despun the rocket for nominal data collection, otherwise the
# magnetic field will be very out of line with the ESA detector



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
wFiles = [[1, 3, 5], [1, 4]]

inputPath_modifier = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_magPitch = 'calibration\ESA_magPitch_calibration' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'calibration\ESA_magPitch_calibration' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
modifier = ''

# --- PITCH BIN RESORTING TOGGLES ---
degWidth = 5 # controls the size of the acceptance range for re-bining

# --- CHI SQUARE CALIBRATION TOGGLES ---
outputPath_modifier_chiCal = 'calibration/ESA_ChiSquare_calibration'

# --- Output Data TOGGLES ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from itertools import product


def L1_to_L1magCalESA(wRocket, wFile, rocketFolderPath, justPrintFileNames):

    # --- ACES II Flight/Integration Data ---
    rocketID = ACESII.payload_IDs[wRocket-4]
    globalAttrsMod = ACESII.global_attributes[wRocket-4]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L1'
    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}{modifier}\*.cdf')
    inputFiles_magPitch = glob(f'{rocketFolderPath}{inputPath_modifier_magPitch}\{ACESII.fliers[wRocket-4]}{modifier}\*magPitch*')
    outputFiles = glob(f'{rocketFolderPath}{outputPath_modifier}\{ACESII.fliers[wRocket-4]}\*.cdf')
    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}{modifier}\\', '') for ifile in inputFiles]
    input_names_magPitch = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}{modifier}\\', '') for ifile in inputFiles]
    output_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{ACESII.fliers[wRocket-4]}\\', '') for ofile in outputFiles]
    input_names_searchable = [ifile.replace('ACESII_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]


    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}{modifier}\\', '')

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in dataFile_name:
            wInstr = [index, instr]

    fileoutName = f'ACESII_{rocketID}_l1_{wInstr[1]}_magCal.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(stl.color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + stl.color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the l1 ESA file ---
        stl.prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict_esa = stl.loadDictFromFile(inputFiles[wFile])
        data_dict_esa['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_esa['Epoch'][0][i]) for i in range(len(data_dict_esa['Epoch'][0]))])
        stl.Done(start_time)

        # --- get the data from the MagPitch file ---
        stl.prgMsg(f'Loading data from MagPitch Files')

        this_magPitchFile = ''
        for file in inputFiles_magPitch:
            if wInstr[1] in file and 'magPitch' in file:
                this_magPitchFile = file

        data_dict_magPitch = stl.loadDictFromFile(this_magPitchFile)
        data_dict_magPitch['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_magPitch['Epoch'][0][i]) for i in (range(len(data_dict_magPitch['Epoch'][0])))])
        stl.Done(start_time)

        #############################
        # --- Re-Bin the ESA data ---
        #############################
        stl.prgMsg('Rebinnig the ESA data')

        # reduce dataset to only include data after ACS despun rocket since this is where the magPitchAngles are nominal. Value is chosen based off of the attitude data variable "SpinAlignmentToB"
        targetIndexTimes = [
            np.abs(data_dict_esa['Epoch'][0] - pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 21, 48, 300))).argmin(),
            np.abs(data_dict_esa['Epoch'][0] - pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 23, 46, 750))).argmin()
        ]

        esaRanges = [range(len(data_dict_esa[wInstr[1]][0])), range(len(data_dict_esa[wInstr[1]][0][0])), range(len(data_dict_esa[wInstr[1]][0][0][0]))]
        esaData = data_dict_esa[wInstr[1]][0]
        magPitch = data_dict_magPitch['Mag_Calculated_Pitch_Angle'][0]
        padAngles = data_dict_esa['Pitch_Angle'][0]

        # --- prepare data for re-binning ---
        esaDataSorted = [[[[] for engy in esaRanges[2]] for ptch in esaRanges[1]] for tme in esaRanges[0]]

        # --- --- --- --- --- --- --- --- --- --
        # --- COLLECT CHI^2 CALIBRATION DATA ---
        # --- --- --- --- --- --- --- --- --- --

        # apended data will look like: [uncalPoint, prinPoint, uncalPad_index, prinPad_index]
        calData = []

        # for two datapoints at the same time/energy but different pitch pads whose "True" pitch angle are <5deg differenet,
        # store this data in calData

        # --- rebin the data and collect calibration ---
        for tme, ptch, engy in tqdm(product(*esaRanges)):

            pitchVal = magPitch[tme][ptch][engy]
            esaVal = esaData[tme][ptch][engy]

            if wInstr[1] != 'iepaa':
                if (0 < ptch < 20): # 0deg - 180 deg case
                    pad = padAngles[ptch]
                elif ptch == 20: # 190 deg case
                    pad = 170
                elif ptch == 0: # -10 deg case
                    pad = 10
            else:
                pad = padAngles[ptch]

            if np.abs(pad - pitchVal) <= degWidth: # if you're within the right pad, place it where it was
                esaDataSorted[tme][ptch][engy].append(esaVal)

            elif np.abs(pad - pitchVal) > degWidth: # if pitchVal is further than "degWidth" away from its assigned pad # find which bin to put the value in   # put it in the right bin

                # add step here b/c we don't want the -10deg bin putting data into 20 deg or higher, same for 190deg bin
                # Thus, we ONLY allow -10deg bin to give to 0deg bin and 190deg bin to 180deg
                index = np.abs(padAngles - pitchVal).argmin()
                if ptch not in [0, 20]: # not the -10deg, 190 deg case
                    esaDataSorted[tme][index][engy].append(esaVal)
                elif ptch == 0 and index == 1:
                    esaDataSorted[tme][index][engy].append(esaVal)
                elif ptch == 20 and index == 19:
                    esaDataSorted[tme][index][engy].append(esaVal)

                # ChiSquare CALIBRATION collection: check if pitchval_uncal needs to be moved but pitchval_prin doesn't
                # logic: Since we are looping through all vals, we may have a situation where
                # the uncal value right now is going to be paired with a principal values
                # that also needs to be moved. We don't want to pair these. Instead, ensure
                # the paired prin value is not going to move i.e. it will be within its pad range
                if np.abs(magPitch[tme][index][engy] - padAngles[index]) <= degWidth:

                    # check if values occurs during the time after ACS despun rocket for maiden voyage.
                    if tme >= targetIndexTimes[wRocket - 4]:
                        # [uncalPoint, prinPoint, uncalPad_index, prinPad_index]
                        if ptch not in [0, 20]:  # not the -10deg, 190 deg case
                            calData.append([esaData[tme][ptch][engy], esaData[tme][index][engy], ptch, index])
                        elif ptch == 0 and index == 1:
                            calData.append([esaData[tme][ptch][engy], esaData[tme][index][engy], ptch, index])
                        elif ptch == 20 and index == 19:
                            calData.append([esaData[tme][ptch][engy], esaData[tme][index][engy], ptch, index])


        stl.Done(start_time)


        # --- --- --- --- --- --- --- --- --- --- -
        # --- REDUCE DATALISTS TO SINGLE VALUES ---
        # --- --- --- --- --- --- --- --- --- --- -
        stl.prgMsg('Reducing DataSet')
        for tme, ptch, engy in product(*esaRanges): # take the average of all the lists in the data
            if len(esaDataSorted[tme][ptch][engy]) != 0:
                # remove zeros from the average
                dataToAverage = np.array(esaDataSorted[tme][ptch][engy])
                dataToAverage= dataToAverage[np.where(dataToAverage !=0)]

                if len(dataToAverage) == 0:
                    esaDataSorted[tme][ptch][engy] = 0
                else:
                    esaDataSorted[tme][ptch][engy] = int(round(sum(dataToAverage) / len(dataToAverage)))

                    if 0 < int(round(sum(dataToAverage)/len(dataToAverage))) <= 2:
                        print(tme, ptch, engy, int(round(sum(dataToAverage)/len(dataToAverage))), esaDataSorted[tme][ptch][engy])

            else:
                esaDataSorted[tme][ptch][engy] = ACESII.epoch_fillVal

        esaDataSorted = np.array(esaDataSorted)

        stl.Done(start_time)

        # --- --- --- --- --- --- --- ---
        # --- PREPARE DATA FOR OUTPUT ---
        # --- --- --- --- --- --- --- ---

        if outputData:
            data_dict = {}
            data_dict = {**data_dict, **{f'{wInstr[1]}':
                                             [esaDataSorted, {'LABLAXIS': f'{wInstr[1]}',
                                                          'DEPEND_0': 'Epoch', 'DEPEND_1': 'Pitch_Angle',
                                                          'DEPEND_2': 'Energy',
                                                          'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2',
                                                          'UNITS': 'counts',
                                                          'VALIDMIN': esaDataSorted.min(), 'VALIDMAX': esaDataSorted.max(),
                                                          'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
            # get Energy and Pitch angle into data_dict
            for key, val in data_dict_esa.items():
                if key not in [f'{wInstr[1]}']:
                    data_dict = {**data_dict, **{key:data_dict_esa[key]}}

            # --- --- --- --- --- --- --- --- ---
            # --- WRITE OUT THE L1 MagCal DATA ---
            # --- --- --- --- --- --- --- --- ---
            stl.prgMsg('Creating output file')
            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
            globalAttrsMod['Descriptor'] = ACESII.ESA_names[wInstr[0]]
            stl.outputCDFdata(outputPath, data_dict,globalAttrsMod= globalAttrsMod,instrNam= wInstr[1])

            stl.Done(start_time)

            # --- --- --- --- --- --- --- --- --- -
            # --- WRITE OUT THE CHISQR CAL DATA ---
            # --- --- --- --- --- --- --- --- --- -
            # store the raw ChiSquare calibration data
            del data_dict
            calData = np.array(calData) if len(calData) > 0 else np.array([[0, 0, 0, 0]])

            data_dict = {f'ChiSquare_calData': [calData, {'LABLAXIS': 'ChiSquareCalData',
                                                        'DEPEND_0': None, 'DEPEND_1': None,
                                                        'DEPEND_2': None,
                                                        'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2',
                                                        'UNITS': 'counts',
                                                        'VALIDMIN': calData.min(),
                                                        'VALIDMAX': calData.max(),
                                                        'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}

            fileoutName = f'ACESII_{rocketID}_{wInstr[1]}_ChiSquareData.cdf'
            outputPath = f'{rocketFolderPath}{outputPath_modifier_chiCal}\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
            globalAttrsMod['Descriptor'] = ACESII.ESA_names[wInstr[0]]

            stl.outputCDFdata(outputPath, data_dict,globalAttrsMod= globalAttrsMod, instrNam= wInstr[1])


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = DataPaths.ACES_data_folder

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}\*.cdf')) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    if justPrintFileNames:
        L1_to_L1magCalESA(wRocket, 0, rocketFolderPath, justPrintFileNames)
    elif not wFiles[wRocket-4]:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}\*.cdf')))):
            L1_to_L1magCalESA(wRocket, fileNo, rocketFolderPath, justPrintFileNames)
    else:
        for filesNo in wFiles[wRocket-4]:
            L1_to_L1magCalESA(wRocket, filesNo, rocketFolderPath, justPrintFileNames)