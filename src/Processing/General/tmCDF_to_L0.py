# --- tmCDF_to_L0.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Read in CDF telemetry files and processes them to L0. ALSO, peel out magnetometer data and
# store it for further processing


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files in
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
# [#1,#2,...etc] --> only specific files
wFiles = [5]

inputPath_modifier = 'tmCDF' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'L0' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# --- magnetometers ---
getMAGdata = True # get the mag data
justgetMAGdata = True # Only get the MAG data
useAntoniosDataInstead = True # use Antonios data instead of my dataset
wMag = 0 # 0 --> RingCore, 1--> CHIMERA/Tesseract (depends on flyer)
modifer_mag = 'L0'
modifier_antonio_data = 'L0'

# Output the data
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from warnings import filterwarnings # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.
filterwarnings("ignore")

def tmCDF_to_L0(wRocket, wFile, rocketFolderPath, justPrintFileNames,wflyer):

    # --- --- --- --- --- --- ---
    # --- FILE SPECIFIC DATA ---
    # --- --- --- --- --- --- ---

    # --- ACES II Flight/Integration Data ---
    if wRocket in [0,1,4,5]:
        rocketAttrs, data_dicts, deConvolveKeys = ACES_mission_dicts()
        rocketID = rocketAttrs.rocketID[wflyer]
        globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
        globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + f'{outputPath_modifier}'
        L0ModelData = L0_TRICE_Quick(wflyer)

    # --- TRICE II ---
    elif wRocket in [2,3]:
        globalAttrsMod = {}
        rocketAttrs = TRICE_mission_dicts()
        rocketID = rocketAttrs.rocketID[wflyer]
        L0ModelData = L0_TRICE_Quick(wflyer)

    # Set the paths for the file names
    tmCDFFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')
    L0Files = glob(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*.cdf')

    tmCDF_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\\', '') for ifile in tmCDFFiles]
    L0_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '').replace(".cdf", '') for ofile in L0Files]

    tmCDF_names_searchable = [cdfname.replace('.cdf', '').replace('_', '').replace('36359', '').replace('36364', '') for cdfname in tmCDF_names]
    L0_names_searchable = [fname.replace('ACES_','').replace('l0_','').replace('eepaa_','').replace('leesa_','').replace('iepaa_','').replace('_v00','').replace('_','').replace('eepaa','').replace('leesa','').replace('iepaa','').replace('36359','').replace('36364','') for fname in L0_names]

    dataFile_name = tmCDFFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\\', '')

    # Output Naming Format: EEPAA, LEESA, IEPAA
    if wRocket in [0, 1]: # Integration Files
        dF_name = dataFile_name.replace('36_359_', '').replace('36_364_', '').replace('.cdf', '')
        fileoutName = [f'{rocketAttrs.missionNam}_{rocketID}_l0_{rocketAttrs.InstrNames_LC[i]}_{dF_name}.cdf' for i in (range(len(rocketAttrs.InstrNames_LC)))]
    elif wRocket in [4, 5]: #ACESII flight files
        dF_name = dataFile_name.replace('36359_', '').replace('36364_', '').replace('ACESII_','').replace('flight_',"").replace('_v00','').replace('.cdf', '')
        fileoutName = [f'{rocketAttrs.missionNam}_{rocketID}_l0_{rocketAttrs.InstrNames_LC[i]}_{dF_name}.cdf' for i in (range(len(rocketAttrs.InstrNames_LC)))]
        if justgetMAGdata:
            if wRocket==4:
                wInstr = [wMag, 'RingCore'] if wMag == 0 else [wMag, 'Chimera']
            elif wRocket == 5:
                wInstr = [wMag, 'RingCore'] if wMag == 0 else [wMag, 'Tesseract']
            if 'mag_cal' in dF_name:
                fileoutName = f'{rocketAttrs.missionNam}_{rocketID}_l0_{wInstr[1]}_magCal_20221013_magFrm.cdf'
            else:
                fileoutName = f'{rocketAttrs.missionNam}_{rocketID}_l0_{wInstr[1]}_magFrm.cdf'
    else: #TRICE
        dF_name = dataFile_name.replace(rocketID,'').replace(rocketAttrs.missionNam, '').replace('_', '').replace('k0','').replace('.cdf', '')
        fileoutName = [f'{rocketAttrs.missionNam}_{rocketID}_l0_{rocketAttrs.InstrNames_LC[i]}_{dF_name}.cdf' for i in (range(len(rocketAttrs.InstrNames_LC)))]

    # --- --- --- --- --- --- --- ------
    # --- PRODUCE TELEM DATA FILES ---
    # --- --- --- --- --- --- --- ------
    if justPrintFileNames:
        for i, file in enumerate(tmCDFFiles):
            anws = ["yes" if tmCDF_names_searchable[i].replace('.cdf', "") in L0_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB   Made CDF: {:3s} '.format(i, tmCDF_names[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(tmCDFFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the tmCDF file ---
        prgMsg('Getting data from tmCDFfile')
        with pycdf.CDF(tmCDFFiles[wFile]) as tmCDFDataFile:
            epoch_from_file = tmCDFDataFile['Epoch'][...]
            atmCDF_epoch = [pycdf.lib.datetime_to_tt2000(epoch_from_file[i]) for i in range(len(epoch_from_file))]
            atmCDF_sfid = tmCDFDataFile['sfid'][...]
            atmCDF_mf = tmCDFDataFile['minorframe'][...]
        Done(start_time)

        tmCDF_epoch = atmCDF_epoch
        tmCDF_sfid = atmCDF_sfid
        tmCDF_mf = atmCDF_mf

        # --- --- --- --- --- ----
        # --- COLLECTING DATA ---
        # --- --- --- --- --- ----

        if not justgetMAGdata:
            # --- pre-allocate memory for ESA data---
            No_of_esa_measurements = int(len(tmCDF_sfid)/4)
            for i in (range(rocketAttrs.NumOfInstr - 1)):
                data_dicts[i]['Epoch'][0] = np.zeros(shape=(No_of_esa_measurements), dtype='int64')
                data_dicts[i]['sfid'][0] = np.zeros(shape=(No_of_esa_measurements), dtype='uint8')
                data_dicts[i]['Sector_Counts'][0] = np.zeros(shape=(No_of_esa_measurements, rocketAttrs.num_of_sector_counts[i]), dtype='int32')

        ######################
        # --- GET MAG DATA ---
        ######################

        if getMAGdata:
            prgMsg('Collecting MAG data')

            if not useAntoniosDataInstead:
                if wMag == 0:

                    # Get the data
                    RingCore_data_temp = np.array([tmCDF_mf[iv][120 - 1] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1)%2 == 1])
                    RingCore_sfid_temp = np.array([tmCDF_sfid[iv] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1)%2 == 1])
                    RingCore_epoch_temp = np.array([tmCDF_epoch[iv] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1)%2 == 1])
                    Done(start_time)

                    # package all the MAG data into a data_dict
                    ringCore_data = {
                        'Epoch':[],
                        'HouseKeeping_ID': [],
                        'HouseKeeping_Data': [],
                        'TimeOffset_High': [],
                        'TimeOffset_Low': [],
                        'Bx_high': [],
                        'Bx_low': [],
                        'By_high': [],
                        'By_low': [],
                        'Bz_high': [],
                        'Bz_low': [],
                    }

                    prgMsg('Collecting mag Data')

                    for i in tqdm(range(len(RingCore_sfid_temp))):

                        # I'm CHOOSING TO PICK THE EPOCH AT WORD #0 and WORD #11 in a Major Frame
                        if RingCore_sfid_temp[i] == 0: # Epoch
                            ringCore_data['Epoch'].append(RingCore_epoch_temp[i])
                        elif RingCore_sfid_temp[i] == 2: # BxHigh
                            ringCore_data['Bx_high'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 4: #BxLow
                            ringCore_data['Bx_low'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 6: # ByHigh
                            ringCore_data['By_high'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 8: #ByLow
                            ringCore_data['By_low'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 10: # BzHigh
                            ringCore_data['Bz_high'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 12: # BzLow
                            ringCore_data['Bz_low'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 14: # HOUSEKEEPING ID
                            ringCore_data['HouseKeeping_ID'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 16: # Housekeeping data
                            ringCore_data['HouseKeeping_Data'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 18:
                            ringCore_data['TimeOffset_High'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 20:
                            ringCore_data['Epoch'].append(RingCore_epoch_temp[i])
                            ringCore_data['TimeOffset_Low'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 22:  # BxHigh
                            ringCore_data['Bx_high'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 24:  # BxLow
                            ringCore_data['Bx_low'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 26:  # ByHigh
                            ringCore_data['By_high'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 28:  # ByLow
                            ringCore_data['By_low'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 30:  # BzHigh
                            ringCore_data['Bz_high'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 32:  # BzLow
                            ringCore_data['Bz_low'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 34:  # HOUSEKEEPING ID
                            ringCore_data['HouseKeeping_ID'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 36:  # Housekeeping data
                            ringCore_data['HouseKeeping_Data'].append(RingCore_data_temp[i])
                        elif RingCore_sfid_temp[i] == 38:
                            ringCore_data['TimeOffset_High'].append(RingCore_data_temp[i])
                            ringCore_data['TimeOffset_Low'].append(rocketAttrs.epoch_fillVal)

                    # --- --- --- --- --- --- --- --
                    # --- calculate 24-Bit Words ---
                    # --- --- --- --- --- --- --- --

                    # create a container to hold all the processed data
                    data_dict_mag = {
                        'Bx': [[], {'LABLAXIS':'Bx', 'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'I5', 'UNITS': 'ADC', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data', 'SCALETYP': 'linear'}],
                        'By': [[], {'LABLAXIS':'By', 'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'I5', 'UNITS': 'ADC', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data', 'SCALETYP': 'linear'}],
                        'Bz': [[], {'LABLAXIS':'Bz', 'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'I5', 'UNITS': 'ADC', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data', 'SCALETYP': 'linear'}],
                        'Epoch': [ringCore_data['Epoch'], {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'I5', 'UNITS': 'ns', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000', 'TIME_SCALE': 'Terrestrial Time', 'REFERENCE_POSITION': 'Rotating Earth Geoid', 'SCALETYP': 'linear'}],
                        'HouseKeeping_ID': [[], {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'I5', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}],
                        'HouseKeeping_Data': [[], {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'I5', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}],
                        'TimeOffset': [[], {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'I5', 'UNITS': 'ns', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]
                    }

                    labels = ['Bx', 'By', 'Bz']

                    for i in range(len(data_dict_mag['Epoch'][0])):

                        # get the housekeeping data
                        hkID = '{0:016b}'.format(ringCore_data['HouseKeeping_ID'][i])
                        hkData = '{0:016b}'.format(ringCore_data['HouseKeeping_Data'][i])
                        data_dict_mag['HouseKeeping_ID'][0].append(int(hkID[9:16], 2))
                        data_dict_mag['HouseKeeping_Data'][0].append(int(hkData, 2))

                        # TimeOffset - 18-bit unsigned word 2+16. The 18-bit word converts by 61.035ns/count
                        # to gets you time in nanoseconds SINCE LAST MAJORFRAME
                        word_high_offset = ringCore_data['TimeOffset_High'][i]
                        word_low_offset = ringCore_data['TimeOffset_Low'][i]

                        if word_low_offset == rocketAttrs.epoch_fillVal:
                            data_dict_mag['TimeOffset'][0].append(rocketAttrs.epoch_fillVal)
                        else:
                            data_dict_mag['TimeOffset'][0].append(int('{0:016b}'.format(word_high_offset)[15:16] + '{0:016b}'.format(word_low_offset), 2) * 61.035 ) # converts to ns since majorframe

                        # GET THE MANGETIC FIELD - It's a 24-bit SIGNED word, so we must apply 2's complement to the whole word if it goes negative
                        for j in range(3):
                            # If the 24th bit is 1 --> positive value, else negative. Convert to decimal number
                            # Take the high value, only use first 7 bits, then concatenate it with the all of the low value bits
                            word_high = ringCore_data[labels[j] + '_high'][i]
                            word_low = ringCore_data[labels[j] + '_low'][i]
                            sign = '{0:016b}'.format(word_high)[-8]

                            if sign == '1': # 1 is negative (opposite to what antonio said)
                                word_high_twosComp = '{0:016b}'.format(word_high)
                                word_high_twosComp = word_high_twosComp.replace('1', '2').replace('0', '1').replace('2', '0') # the 2's complement equivalent
                                word_low_twosComp = '{0:016b}'.format(word_low)
                                word_low_twosComp = word_low_twosComp.replace('1', '2').replace('0', '1').replace('2', '0')
                                fullword = word_high_twosComp[-7:] + word_low_twosComp
                                fullword = -1 * int(fullword, 2)

                            elif sign == '0': # 0 is postivie (opposite to what antonio said)
                                fullword = '{0:016b}'.format(word_high)[-7:] + '{0:016b}'.format(word_low)
                                fullword = int(fullword, 2)

                            data_dict_mag[labels[j]][0].append(fullword)

                    ###################################
                    # --- REMOVE OVERSAMPLING ISSUE ---
                    ###################################

                    # use repeated hk index to determine bad data i.e. the hk index must go from 0-127 without repeats.
                    # there ARE repeats due to the oversampling of the telemetry
                    bad_indicies = []

                    # find all repeated indicies
                    for i in range(len(data_dict_mag['HouseKeeping_ID'][0])-1):
                        if data_dict_mag['HouseKeeping_ID'][0][i+1] == data_dict_mag['HouseKeeping_ID'][0][i]:
                            bad_indicies.append(i+1)

                    # convert everything to numpy arrays and remove repeated indicies
                    for key, val in data_dict_mag.items():

                        if len(data_dict_mag[key][0]) > 400:
                            data_dict_mag[key][0] = np.array(data_dict_mag[key][0])

                            data_dict_mag[key][0] = np.delete(data_dict_mag[key][0],bad_indicies)
                elif wMag == 1:

                    if wRocket == 5: # Low Flyer has Ring Core and Tesseract
                        # Get the data
                        Tesseract_data = np.array([tmCDF_mf[iv][16 - 1] for iv in (range(len(tmCDF_sfid)))])
                        Tesseract_sfid = np.array(tmCDF_sfid)
                        Tesseract_epoch = np.array(tmCDF_epoch)

                        # Process the data
                        data_dict_mag = {}
                        data_dict_mag = {**data_dict_mag, **{'Tesseract_data':[Tesseract_data, {'DEPEND_0': 'Tesseract_epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'data','SCALETYP':'linear'}]}}
                        data_dict_mag = {**data_dict_mag, **{'Tesseract_sfid': [Tesseract_sfid, {'DEPEND_0': 'Tesseract_epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'data','SCALETYP':'linear'}]}}
                        data_dict_mag = {**data_dict_mag, **{'Tesseract_epoch': [Tesseract_epoch, {'DEPEND_0': 'Tesseract_epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': rocketAttrs.epoch_fillVal,'FORMAT': 'I5','UNITS': 'ns','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','MONOTON':'INCREASE','TIME_BASE':'J2000','TIME_SCALE':'Terrestrial Time','REFERENCE_POSITION':'Rotating Earth Geoid','SCALETYP':'linear'}]}}
            else:
                antonioMagFile = glob(f'{rocketFolderPath}{modifier_antonio_data}\{fliers[wflyer]}\*Antonio_*')[0]
                data_dict_mag = loadDictFromFile(antonioMagFile,{})
            Done(start_time)

            if outputData:
                prgMsg('Writing out MAG data')

                outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

                outputCDFdata(outputPath, data_dict_mag, L0ModelData, globalAttrsMod, wInstr[1])

                Done(start_time)

        if not justgetMAGdata:
            # --- variables with special conditions ---
            prgMsg('Collecting variables with special conditions')
            for j in tqdm(range(rocketAttrs.NumOfInstr)):
                if j == 0:
                    data_dicts[j]['Boom_Monitor'][0] = [tmCDF_mf[iv][148 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 1]
                    data_dicts[j]['Epoch_monitors'][0] = [tmCDF_epoch[iv] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 17]
                    data_dicts[j]['28V_Monitor'][0] = [tmCDF_mf[iv][141 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 17]
                    data_dicts[j]['EXP_Current'][0] = [tmCDF_mf[iv][121 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1)%4 == 1]
                elif j == 1:
                    data_dicts[j]['Boom_Monitor'][0] = [tmCDF_mf[iv][148 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 6]
                    data_dicts[j]['Epoch_monitors'][0] = [tmCDF_epoch[iv] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 18]
                    data_dicts[j]['28V_Monitor'][0] = [tmCDF_mf[iv][141 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 18]
                    data_dicts[j]['EXP_Current'][0] = [tmCDF_mf[iv][121 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) % 4 == 0]
                elif j == 2:
                    data_dicts[j]['Boom_Monitor'][0] = [tmCDF_mf[iv][148 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 2]
                    data_dicts[j]['Epoch_monitors'][0] = [tmCDF_epoch[iv] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 33]
                    data_dicts[j]['28V_Monitor'][0] = [tmCDF_mf[iv][141 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 33]
                    data_dicts[j]['EXP_Current'][0] = [tmCDF_mf[iv][121 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) % 4 == 2]
                elif j == 3:
                    data_dicts[j]['Boom_Monitor_1'][0] = [tmCDF_mf[iv][148 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 3]
                    data_dicts[j]['Epoch_monitor_1'][0] = [tmCDF_epoch[iv] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 3]

                    data_dicts[j]['Boom_Monitor_2'][0] = [tmCDF_mf[iv][148 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 4]
                    data_dicts[j]['Epoch_monitor_2'][0] = [tmCDF_epoch[iv] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 4]

                    data_dicts[j]['28V_Monitor'][0] = [tmCDF_mf[iv][141 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 19]
                    data_dicts[j]['Epoch_monitors'][0] = [tmCDF_epoch[iv] for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) == 19]

                    data_dicts[j]['EXP_Current'][0] = [tmCDF_mf[iv][121 - 1] >> 6 for iv in (range(len(tmCDF_sfid))) if (tmCDF_sfid[iv] + 1) % 4 == 3]

            Done(start_time)

            # --- Major Frame counter ---
            prgMsg('Counting major frames')
            for ii in (range(len(tmCDF_sfid))): # NOTE: Major frame counter 1 is the first 16 binary digits of the majorframe number and major frame counter 2 is the next 16 binary digits, which totals a 32 binary word.
                if (tmCDF_sfid[ii] + 1) not in [4, 12, 20, 28, 36]:
                    for j in (range(rocketAttrs.NumOfInstr)):  # MajF counter 1
                        data_dicts[j]['major_frame_counter'][0].append([tmCDF_mf[ii][72 - 1]])
                else:
                    for j in (range(rocketAttrs.NumOfInstr)):  # MajF counter 2
                        data_dicts[j]['major_frame_counter'][0].append([tmCDF_mf[ii][136 - 1], tmCDF_mf[ii][72 - 1]])
            Done(start_time)

            # --- LP Data ---
            prgMsg('Collecting LP Data')
            data_dicts[3]['Channel_Counts'][0] = [[tmCDF_mf[iv][word - 1] for word in rocketAttrs.LP_words] for iv in (range(len(tmCDF_sfid))) ]
            data_dicts[3]['Epoch'][0] = [tmCDF_epoch[iv] for iv in (range(len(tmCDF_sfid)))]
            data_dicts[3]['sfid'][0] = [tmCDF_sfid[iv] for iv in (range(len(tmCDF_sfid)))]
            Done(start_time)

            # --- --- --- --- --- --- --- --- --- --- --- ---
            # --- COLLECT MINORFRAME DATA FOR THE ESAs  ---
            # --- --- --- --- --- --- --- --- --- --- --- ---

            for iii in range(rocketAttrs.NumOfInstr - 1):
                words = rocketAttrs.ESA_words[iii]
                deConvolveKey = deConvolveKeys[iii]
                deConvolveItems = deConvolveKey.items()
                SectorDataCounter = 0

                for iv in tqdm(range(No_of_esa_measurements)):
                    esaMeasurement = tmCDF_mf[4 * iv: 4 * (iv + 1)]
                    data_dicts[iii]['sfid'][0][iv] = tmCDF_sfid[4 * iv]
                    data_dicts[iii]['Epoch'][0][iv] = tmCDF_epoch[4 * iv]

                    for index in range(len(esaMeasurement)):
                        mfdata = [esaMeasurement[index][k - 1] for k in words]
                        keyset = [val for key, val in deConvolveItems if key in [i for i in range(1 + index * len(words), 1 + (1 + index) * len(words))]]

                        for v in (range(len(mfdata))):  # place all the data
                            if keyset[v] == 'Sector_Counts':
                                data_dicts[iii]['Sector_Counts'][0][iv][SectorDataCounter] = mfdata[v]
                                SectorDataCounter += 1

                                if SectorDataCounter == rocketAttrs.num_of_sector_counts[iii]:
                                    SectorDataCounter = 0

                            elif keyset[v] == 'status_word_2':
                                data_dicts[iii]['HV_div16'][0].append((mfdata[v] >> 14) & 1)
                                data_dicts[iii]['HV_enable'][0].append((mfdata[v] >> 7) & 1)
                                data_dicts[iii]['sweep_step'][0].append((mfdata[v] >> 8) & 63)
                                data_dicts[iii]['TP5_enable'][0].append((mfdata[v] >> 5) & 1)

                            elif keyset[v] == '625kHz_Clock_Input':
                                data_dicts[iii][keyset[v]][0].append(mfdata[v])
                                if mfdata[v] / rocketAttrs.ESA_CLK_INPUT >= 2:
                                    data_dicts[iii]['Count_Interval'][0].append(65535)
                                else:
                                    data_dicts[iii]['Count_Interval'][0].append(mfdata[v] / rocketAttrs.ESA_CLK_INPUT)
                            elif keyset[v] != None:
                                data_dicts[iii][keyset[v]][0].append(mfdata[v])

            # --- --- --- --- ---
            # --- PROCESSING ---
            # --- --- --- --- ---
            prgMsg('Processing Data')
            for jjj in range(rocketAttrs.NumOfInstr):

                # --- Major Frame Counter ---
                MF_counter, get_first_major_frame = [], []
                MFDat = data_dicts[jjj]['major_frame_counter'][0]
                topNum, botNum = 0, 0

                for data in MFDat:
                    if len(data) == 1:
                        botNum = data[0]
                    elif len(data) == 2:
                        botNum = data[1]
                        topNum = data[0]

                        if get_first_major_frame == []: # used to start the minor frame counter
                            if topNum != 0:
                                get_first_major_frame.append( (topNum << 16) - 1 ) # -1 to remove the MF that I start on
                            else:
                                get_first_major_frame.append(botNum - 1)

                    num = (topNum << 16) + botNum
                    MF_counter.append(num)

                data_dicts[jjj]['major_frame_counter'][0] = MF_counter

                # --- minor Frame Counter ---
                data_dicts[jjj]['minor_frame_counter'][0] = np.array([i + (get_first_major_frame[0]*40) for i in (range(1,len(tmCDF_sfid)+1))])

                # --- Convert all data to numpy arrays ---
                for key, dataVal in data_dicts[jjj].items():
                    data_dicts[jjj][key][0] = np.array(dataVal[0])


            Done(start_time)
            if outputData:
                # --- --- --- --- --- --- ---
                # --- WRITE OUT THE DATA ---
                # --- --- --- --- --- --- ---
                writeOut_path = [f'{rocketFolderPath}L0\{fliers[wflyer]}\\{fileoutName[i]}' for i in range(rocketAttrs.NumOfInstr)]

                # --- loop through instruments ---
                for v in range(rocketAttrs.NumOfInstr):
                    prgMsg(f'Writing out {rocketAttrs.InstrNames[v]} data')
                    outputPath = writeOut_path[v]

                    # --- delete output file if it already exists ---
                    if os.path.exists(outputPath):
                        os.remove(outputPath)

                    # --- open the output file ---
                    with pycdf.CDF(outputPath,'') as L0File:
                        L0File.readonly(False)

                        # --- write out global attributes ---
                        inputGlobDic = L0ModelData.cdfFile.globalattsget()
                        for key, val in inputGlobDic.items():
                            if key == 'Descriptor':
                                globalAttrsMod[key] = rocketAttrs.InstrNames_Full[v]

                            if key in globalAttrsMod:
                                L0File.attrs[key] = globalAttrsMod[key]
                            else:
                                L0File.attrs[key] = val

                        # --- WRITE OUT DATA ---
                        for varKey, varVal in data_dicts[v].items():
                            if varKey == 'Epoch' or varKey == 'Epoch_monitors':
                                L0File.new(varKey, data=varVal[0], type = 33)
                            elif varKey == 'Sector_Counts':
                                L0File.new(varKey, data=varVal[0], type=pycdf.const.CDF_UINT2)
                            else:
                                L0File.new(varKey, data=varVal[0])

                            # --- Write out the attributes and variable info ---
                            for attrKey, attrVal in data_dicts[v][varKey][1].items():
                                if attrKey == 'VALIDMIN':
                                    L0File[varKey].attrs[attrKey] = varVal[0].min()
                                elif attrKey == 'VALIDMAX':
                                    L0File[varKey].attrs[attrKey] = varVal[0].max()
                                elif attrVal != None:
                                    L0File[varKey].attrs[attrKey] = attrVal

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

if len(glob(f'{rocketFolderPath}tmCDF\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        tmCDF_to_L0(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}tmCDF\{fliers[wflyer]}\*.cdf')))):
            tmCDF_to_L0(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            tmCDF_to_L0(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)