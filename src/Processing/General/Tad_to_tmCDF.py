# --- Tad_to_tmCDF.py ---
# --- Author: C. Feltman (inspired by S. Bounds Tad_to_telem code) ---
# DESCRIPTION: Read in .tad telemetry files and convert them to .cdf files with Epoch, minorframe and sfid data
# In order to run correctly, run Scott's data



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

# Just print the names of files in
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
wFiles = [2]

# Break data into chunkNum pieces to write out
chunkNum = 10

inputPath_modifier = 'tad' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'tmCDF' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder



#########################
# --- SPECIAL TOGGLES ---
#########################


# Reduce the back-end of the data based on the wPercent variable.
reduceFile = False

# select how much of the back-end of the data not to include, given in decimal percents
wPercent = [0.04, 0.08] # nominal 4% and 8% for BSS file


# To fix the gSync Epoch dropouts, I get the previously non-fillval epoch values and add the median delta T to it
deltaT_median = 249984 # nanoseconds, in TT2000

# --- --- --- --- --- --- ---
# --- HIGH FLYER TOGGLES ---
# --- --- --- --- --- --- ---
# requires reduceFile = True

# If processing scott's copy of the high flyer, remove the sfid 40 point.
plugHolesInScottsData = False

# If processing BSS high flyer, use Scott's tmCDF file to fill in a value of sfid == 60
plugHolesInBSSData_High = False

# If processing BSS low flyer, fill the one sfid = -1 value with a fillval
plugHolesInBSSData_Low = False





# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from warnings import filterwarnings # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.
filterwarnings("ignore")
from struct import unpack
from cdflib import cdfwrite
from ACESII_code.class_var_func import newCDFvar,color,tmCDF_TRICE_Quick




# # # ----------------------------------------------------------------------
# wPercent = [0.04, 0.08]
# wflyer = 0
#
# # High
# # writeOut_path = r"D:\Data\ACESII\tmCDF\high\ACESII_36359_flight_11202022_processedtm_20221120T170834_v00.cdf"
# writeOut_path = r"D:\Data\ACESII\tmCDF\high\ACESII_36359_flight_playback_BSS_Card0_rawtm_20221120T170902_v00.cdf"
#
#
# specialCaseRawtm = False
# specialCaseBSS = True
#
# # --- --- --- --- --- --- --- --- ---
# # --- REDUCE AND CLEAN UP DATA ---
# # --- --- --- --- --- --- --- --- ---
# with pycdf.CDF(writeOut_path) as unprocessedCDF:
#     tmCDF_epoch = unprocessedCDF['Epoch'][...]
#     tmCDF_sfid = unprocessedCDF['sfid'][...]
#     tmCDF_mf = unprocessedCDF['minorframe'][...]
#
# # Don't use the last % of data in order to speed up the code
# if wPercent[wflyer] > 0:
#     dontUse_amount = int(len(tmCDF_epoch) * wPercent[wflyer])
#     atmCDF_epoch = tmCDF_epoch[0:-1 * dontUse_amount]
#     atmCDF_sfid = tmCDF_sfid[0:-1 * dontUse_amount]
#     atmCDF_mf = tmCDF_mf[0:-1 * dontUse_amount]
# else:
#     atmCDF_epoch = tmCDF_epoch
#     atmCDF_sfid = tmCDF_sfid
#     atmCDF_mf = tmCDF_mf
#
# # Strip away the ends of the data so it starts/ends with a complete majorframe
# mf_startpoint, mf_endpoint = np.where(atmCDF_sfid == 0), np.where(atmCDF_sfid == 39)
# mf_start, mf_end = mf_startpoint[0][0], mf_endpoint[0][-1] + 1  # pythonic + 1
#
# tmCDF_epoch = np.array([pycdf.lib.datetime_to_tt2000(atmCDF_epoch[i]) for i in (range(mf_start, mf_end))])
# tmCDF_sfid = atmCDF_sfid[mf_start:mf_end]
# tmCDF_mf = atmCDF_mf[mf_start:mf_end]
# Done(start_time)
#
#
# # a0 = 722237393164357000 # epoch where sfid 60 occurs in BSS
# #
# # index = np.abs(tmCDF_epoch - a0).argmin()
# # sfidepoch60 = []
# # for band in range(-20,20):
# #     sfidepoch60.append([index,index+band, tmCDF_sfid[index + band], tmCDF_epoch[index + band]])
# #
# # sfidepoch60 = np.array(sfidepoch60)
# #
# #
# # a0 = 722237125199575984 # epoch where sfid 17 -1 30 31 ... occurs in BSS
# #
# # index = np.abs(tmCDF_epoch - a0).argmin()
# # sfidepochm1 = []
# # for band in range(-20, 20):
# #     sfidepochm1.append([index, index+band, tmCDF_sfid[index + band], tmCDF_epoch[index + band]])
# #
# # sfidepochm1 = np.array(sfidepochm1)
#
#
# prgMsg('Plugging Holes in Data')
#
#
# # --- --- --- --- --- --- --- --- --- ---
# # --- CORRECT GSYNCS AND SFID'S >= 40 ---
# # --- --- --- --- --- --- --- --- --- ---
# if specialCaseRawtm:
#     index40 = np.where(atmCDF_sfid == 40)  # should be at 3764720
#     tmCDF_sfid = np.delete(atmCDF_sfid, index40, axis=0)
#     tmCDF_mf = np.delete(atmCDF_mf, index40, axis=0)
#     tmCDF_epoch = np.delete(atmCDF_epoch, index40, axis=0)
#
# if specialCaseBSS:
#     spindex = 3653578 # index at which the gsync error occurs and 13 values go missing. Starting at 17 to 30
#     tmCDF_sfid = np.insert(tmCDF_sfid, spindex, [-1 for i in range(11)] )
#     tmCDF_mf = np.insert(tmCDF_mf, [spindex + i for i in range(11)], [[-1 for i in range(150)] for j in range(11)], axis= 0 )
#     tmCDF_epoch = np.insert(tmCDF_epoch, spindex, [-9223372036854775808 for i in range(11)])
#
#
#
# fix_indicies = []
# sfid40 = []
# neg1_indicies = []
# counter = tmCDF_sfid[0]
# spcounter = 0
# for index in tqdm(range(len(tmCDF_sfid))):
#
#     if tmCDF_sfid[index] != counter:
#         fix_indicies.append([index, tmCDF_sfid[index], counter, tmCDF_epoch[index]])
#
#         if tmCDF_sfid[index] == -1:
#             tmCDF_epoch[index] = tmCDF_epoch[index - 1] + deltaT_median
#             # tmCDF_sfid[index] = counter
#             for band in range(-10, 10):
#             # band = 0
#                 neg1_indicies.append([index+band, tmCDF_sfid[index + band], counter+band, tmCDF_epoch[index+band]])
#
#         if tmCDF_sfid[index] >= 40:
#             for band in range(-20,20):
#                 sfid40.append([index+band, tmCDF_sfid[index + band], counter + band, tmCDF_epoch[index + band]])
#
#     counter += 1
#     if counter >= 40:
#         counter = 0
#
# fix_indicies = np.array(fix_indicies)
# sfid40 = np.array(sfid40)
# neg1_indicies = np.array(neg1_indicies)
#
# Done(start_time)
# # ----------------------------------------------------------------------
#









def unpackData(sRec, nNomDataLen,nDataLen):
    'Unpacks the record into integer values for all data in a minor frame'
    nData = np.array(unpack(nDataLen, sRec[16:]))
    nSync = np.array(unpack("I", sRec[12:16]))
    newData = [0] * nNomDataLen
    newData[0::2] = nData >> 16 # bit shift down by 16 binary values
    newData[1::2] = nData & 65535  # bitwise AND
    newData[nNomDataLen:] = [(nSync >> 16), nSync & 65535]
    return newData

def Gettimestamp(year, nFrmHdr):
    'Converts frame header into a datetime timestamp (year, frmhdr, type)'
    'where type:"Text" returns a text string, else it returns a datetime tuple'
    doy = (nFrmHdr[0] >> 24 & 15) * 100 + (nFrmHdr[0] >> 20 & 15) * 10 + (nFrmHdr[0] >> 16 & 15)
    hours = (nFrmHdr[0] >> 12 & 15) * 10 + (nFrmHdr[0] >> 8 & 15)
    mins = (nFrmHdr[0] >> 4 & 15) * 10 + (nFrmHdr[0] & 15)
    secs = (nFrmHdr[1] >> 28 & 15) * 10 + (nFrmHdr[1] >> 24 & 15)
    msec = (nFrmHdr[1] >> 20 & 15) * 100 + (nFrmHdr[1] >> 16 & 15) * 10 + (nFrmHdr[1] >> 12 & 15)
    usec = (nFrmHdr[1] >> 8 & 15) * 100 + (nFrmHdr[1] >> 4 & 15) * 10 + (nFrmHdr[1] & 15)

    if hours > 23:
        doy += 1
        hours = hours - 24
    if doy == 0:
        doy = 1
    d = dt.datetime.strptime(year.decode("utf-8") + ' ' + str(doy), '%Y %j')
    t = dt.time(hours, mins, secs, msec * 1000 + usec, tzinfo=None)
    tt = dt.datetime.combine(d, t)

    return pycdf.lib.datetime_to_tt2000(tt)

def tadToTelem(wRocket,wFile,chunkNum,rocketFolderPath,justPrintFileNames,wflyer):

    # --- --- --- --- --- --- ---
    # --- FILE SPECIFIC DATA ---
    # --- --- --- --- --- --- ---
    # ACESII Flight/Integration Data
    if wRocket in [0,1,4,5]:
        rocketAttrs,b,c = ACES_mission_dicts()
        tmCDF_folder_path = rf'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}'
        rocketID = rocketAttrs.rocketID[wflyer]
        globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
        tmModelData = tmCDF_TRICE_Quick(wflyer)

    # TRICE II Flight Data
    elif wRocket in [2,3]:
        tmCDF_folder_path = rf'{rocketFolderPath}\{outputPath_modifier}\{fliers[wflyer]}'
        globalAttrsMod = {}
        rocketAttrs,b,c = TRICE_mission_dicts()
        rocketID = rocketAttrs.rocketID[wflyer]
        tmModelData = tmCDF_TRICE_Quick(wflyer)


    tadFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.tad')
    tmCDFFiles = glob(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*.cdf')
    tad_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\\', '').replace(".tad", '') for ifile in tadFiles]
    tmCDF_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '').replace(".cdf",'') for ofile in tmCDFFiles]
    wdataFile_path = tadFiles[wFile]
    wdataFile_name = tadFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\\', '')
    dataFile_name = wdataFile_name

    # --- --- --- --- --- --- --- ------
    # --- PRODUCE TELEM DATA FILES ---
    # --- --- --- --- --- --- --- ------
    if justPrintFileNames:
        for i, file in enumerate(tadFiles):
            anws = ["yes" if tad_names[i].replace('.tad', "") in tmCDF_names else "no" ]
            print('[{:.0f}] {:55s}{:5.1f} MB   Made CDF: {:3s} '.format(i, tad_names[i], round(getsize(file) / (10 ** 6), 1),anws[0]))
    else:
        print(color.UNDERLINE  + f'Converting Byte data for {dataFile_name[0:]}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(tadFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # Get the number of times we'll need to loop. This processes takes only a few seconds
        with open(wdataFile_path,'rb') as tadFile:
            LenCounter = 0
            tadFileHdr = tadFile.read(10 + 12 + 22 + 260 + 12 + 4 + 4 + 4) # 328 byte header
            year = tadFileHdr[22:44].strip(b'\x00').split(b'/')[2][0:4]
            bytes = tadFile.read(rocketAttrs.g_nNomRecLen)
            ts = Gettimestamp(year, unpack("3I", bytes[:12]))  # gets the TT2000 value for the timestamp in the header

            while True:
                bytes = tadFile.read(rocketAttrs.g_nNomRecLen)
                LenCounter += 1
                if not bytes:  # if you run out of data break the loop
                    break


        # --- Open the tad file to read it and output all its data ---
        with open(wdataFile_path, 'rb') as tadFile:
            records = 0
            tadFileHdr = tadFile.read(10 + 12 + 22 + 260 + 12 + 4 + 4 + 4)

            if wRocket in [0, 1]:
                fileout = f'{dataFile_name.replace(".tad", "")}.cdf'
            else:
                fileout = f'{tad_names[wFile]}_rawtm.cdf'

            # --- --- --- --- --- ---
            # --- OUTPUT FILES ---
            # --- --- --- --- --- ---
            writeOut_path = f'{tmCDF_folder_path}\\{fileout}'
            outputFile = cdfwrite.CDF(writeOut_path, cdf_spec=tmModelData.info, delete=True)

            # --- write out global attributes ---
            inputGlobDic = tmModelData.cdfFile.globalattsget()
            globalAttrs = {}
            for key, val in inputGlobDic.items():
                if key in globalAttrsMod:
                    globalAttrs[key] = {0: globalAttrsMod[key]}
                else:
                    globalAttrs[key] = {0: val}

            outputFile.write_globalattrs(globalAttrs)

            # --- Initialize the data ---
            zEpoch = np.array([0])
            zsfid = np.array([0])
            zMF = np.array([[0 for i in (range(rocketAttrs.g_nWords))]])

            modParams = {'FILLVAL': rocketAttrs.epoch_fillVal}
            newCDFvar('Epoch', zEpoch, modParams, tmModelData).writeToFile(outputFile)

            modParams = {'VALIDMIN': 0, 'VALIDMAX': 39, 'FILLVAL': -1, 'Data_Type': 2}
            newCDFvar('sfid', zsfid, modParams, tmModelData).writeToFile(outputFile)

            modParams = {'VALIDMIN': 0, 'VALIDMAX': 65535, 'Dim_Sizes': [rocketAttrs.g_nWords], 'FILLVAL': -1}
            newCDFvar('minorframe', zMF, modParams, tmModelData).writeToFile(outputFile)

            outputFile.close()

            # --- --- --- --- ---
            # --- Main loop ---
            # --- --- --- --- ---
            emin, emax = (10**(20)),0 #these are just used to keep track of the global min/max of the file
            lostgysn = []
            chunkSizes = []
            sfid40 = []

            # loop through the chunks and output data
            for i in range(chunkNum):

                # --- create placeholder variables ---
                chunk_size = round(LenCounter/chunkNum)
                chunkSizes.append(chunk_size)
                Epoch = np.zeros(shape=(chunk_size))
                MF = np.zeros(shape=(chunk_size, rocketAttrs.g_nWords))
                sfid = np.zeros(shape=(chunk_size))

                # --- --- --- --- --- --- ---
                # --- COLLECT THE DATA ---
                # --- --- --- --- --- --- ---
                # Loop through the tad file by reading g_nNomRecLen sized pieces a "chunk_size" number of times

                for loopNum in tqdm(range(chunk_size)):

                    bytes = tadFile.read(rocketAttrs.g_nNomRecLen)

                    # If there's no more bytes to read, exit the loop and write out the data
                    if not bytes:
                        break

                    # If you lose gSync fill in the data location with garbage data
                    elif bytes[12:16] != rocketAttrs.gSync:
                        lostgysn.append([i, loopNum])
                        Epoch[loopNum] = rocketAttrs.epoch_fillVal
                        MF[loopNum] = [-1 for j in (range(rocketAttrs.g_nWords))]
                        sfid[loopNum] = -1

                    # Data didn't lose gSync
                    else:

                        # See if time data is corrupted
                        try:
                            Epoch[loopNum] = Gettimestamp(year, unpack("3I", bytes[:12]))
                        except:
                            Epoch[loopNum] = rocketAttrs.epoch_fillVal

                        MF[loopNum] = unpackData(bytes, rocketAttrs.nNomDataLen,rocketAttrs.nDataLen)
                        (nMF,) = unpack("<H", bytes[10:12])
                        sfid[loopNum] = nMF

                        if nMF >= 40:
                            sfid40.append([chunkNum, loopNum, nMF])

                asfid = np.array(sfid)
                aEpoch = np.array(Epoch)
                aMF = np.array(MF)
                records += len(aEpoch)


                # --- --- --- --- --- --- ---
                # --- WRITE OUT THE DATA ---
                # --- --- --- --- --- --- ---

                with pycdf.CDF(writeOut_path) as outputFile:
                    outputFile.readonly(False)

                    if i == 0: # Overwrite the initialized data
                        outputFile['Epoch'][0:] = aEpoch
                        outputFile['sfid'][0:] = asfid
                        outputFile['minorframe'][0:] = aMF
                    else:
                        outputFile['Epoch'][len(outputFile['Epoch']):] = aEpoch
                        outputFile['sfid'][len(outputFile['sfid']):] = asfid
                        outputFile['minorframe'][len(outputFile['minorframe']):] = aMF

                    if len(aEpoch) != 0:
                        epochSet = sorted(set(aEpoch))
                        if len(epochSet) > 1:
                            if epochSet[1] < emin: #use index 1 to avoid the very negative fillval
                                outputFile['Epoch'].attrs['VALIDMIN'] = epochSet[1]
                                emin = epochSet[1]
                            if epochSet[-1] > emax:
                                outputFile['Epoch'].attrs['VALIDMAX'] = epochSet[-1]
                                emax = epochSet[-1]

                    if i == (chunkNum - 1):
                        Done(start_time)
                        print(f'Number of records written: {records}\n')


            print('--- No of lost Gsyncs ---')
            print(len(lostgysn))

            print('--- No of data that has sfid >=40 ---')
            print(len(sfid40))


            # --- --- --- --- --- ---
            # --- RESIZING DATA ---
            # --- --- --- --- --- ---

            if reduceFile: # special clean-up processes for flight data
                prgMsg('Resizing Data')

                # --- Initialize the output file ---
                fileout = f'{tad_names[wFile]}_fixedtm.cdf'
                outputPath = f'{tmCDF_folder_path}\\{fileout}'

                outputFile = cdfwrite.CDF(outputPath, cdf_spec=tmModelData.info, delete=True)

                # --- write out global attributes ---
                inputGlobDic = tmModelData.cdfFile.globalattsget()
                globalAttrs = {}
                for key, val in inputGlobDic.items():
                    if key in globalAttrsMod:
                        globalAttrs[key] = {0: globalAttrsMod[key]}
                    else:
                        globalAttrs[key] = {0: val}

                outputFile.write_globalattrs(globalAttrs)

                # --- Initialize the data ---
                zEpoch = np.array([0])
                zsfid = np.array([0])
                zMF = np.array([[0 for i in (range(rocketAttrs.g_nWords))]])

                modParams = {'FILLVAL': rocketAttrs.epoch_fillVal}
                newCDFvar('Epoch', zEpoch, modParams, tmModelData).writeToFile(outputFile)

                modParams = {'VALIDMIN': 0, 'VALIDMAX': 39, 'FILLVAL': -1, 'Data_Type': 2}
                newCDFvar('sfid', zsfid, modParams, tmModelData).writeToFile(outputFile)

                modParams = {'VALIDMIN': 0, 'VALIDMAX': 65535, 'Dim_Sizes': [rocketAttrs.g_nWords], 'FILLVAL': -1}
                newCDFvar('minorframe', zMF, modParams, tmModelData).writeToFile(outputFile)
                outputFile.close()

                with pycdf.CDF(outputPath) as outputFile:
                    outputFile.readonly(False)

                    # --- --- --- --- --- --- --- --- ---
                    # --- REDUCE AND CLEAN UP DATA ---
                    # --- --- --- --- --- --- --- --- ---
                    with pycdf.CDF(writeOut_path) as unprocessedCDF:
                        tmCDF_epoch = unprocessedCDF['Epoch'][...]
                        tmCDF_sfid = unprocessedCDF['sfid'][...]
                        tmCDF_mf = unprocessedCDF['minorframe'][...]

                    # Don't use the last % of data in order to speed up the code
                    if wPercent[wflyer] > 0:
                        dontUse_amount = int(len(tmCDF_epoch) * wPercent[wflyer])
                        atmCDF_epoch = tmCDF_epoch[0:-1 * dontUse_amount]
                        atmCDF_sfid = tmCDF_sfid[0:-1 * dontUse_amount]
                        atmCDF_mf = tmCDF_mf[0:-1 * dontUse_amount]
                    else:
                        atmCDF_epoch = tmCDF_epoch
                        atmCDF_sfid = tmCDF_sfid
                        atmCDF_mf = tmCDF_mf

                    # Strip away the ends of the data so it starts/ends with a complete majorframe
                    mf_startpoint, mf_endpoint = np.where(atmCDF_sfid == 0), np.where(atmCDF_sfid == 39)
                    mf_start, mf_end = mf_startpoint[0][0], mf_endpoint[0][-1] + 1  # pythonic + 1

                    tmCDF_epoch = np.array([pycdf.lib.datetime_to_tt2000(atmCDF_epoch[i]) for i in (range(mf_start, mf_end))])
                    tmCDF_sfid = atmCDF_sfid[mf_start:mf_end]
                    tmCDF_mf = atmCDF_mf[mf_start:mf_end]
                    Done(start_time)

                    # --- fix Scott's High Flyer ---
                    if plugHolesInScottsData:

                        prgMsg('Plugging Holes in Data')

                        # --- --- --- --- --- --- --- --- --- ---
                        # --- CORRECT GSYNCS AND SFID'S >= 40 ---
                        # --- --- --- --- --- --- --- --- --- ---

                        # special case to clean up the high flyer data. Remove the single data point where sfid == 40
                        index40 = np.where(tmCDF_sfid == 40)  # should be at 3764720
                        tmCDF_sfid = np.delete(tmCDF_sfid, index40, axis=0)
                        tmCDF_mf = np.delete(tmCDF_mf, index40, axis=0)
                        tmCDF_epoch = np.delete(tmCDF_epoch, index40, axis=0)

                        fix_indicies = []
                        counter = tmCDF_sfid[0]
                        for index in tqdm(range(len(tmCDF_sfid))):

                            if tmCDF_sfid[index] != counter:
                                fix_indicies.append([tmCDF_sfid[index], counter])

                                if tmCDF_sfid[index] == -1:  # fabricate the gSync'd Epochs. No worry about starting with a -1 b/c of the code above
                                    tmCDF_epoch[index] = tmCDF_epoch[index - 1] + deltaT_median
                                    tmCDF_sfid[index] = counter
                                    tmCDF_mf[index] = [-1 for i in range(150)]

                            counter += 1
                            if counter >= 40:
                                counter = 0

                        Done(start_time)

                    # --- Fix BSS High Flyer ---
                    if plugHolesInBSSData_High: # special case that inserts 12 fillvals at index 3653578 and inserts a value from Scotts data at sfid == 60, index 4725425
                        prgMsg('Plugging Holes in Data')
                        # --- PROBLEM #1 sfid 60---

                        path_to_processed_scotts_file = rf"{ACES_data_folder}\{outputPath_modifier}\{wflyer[0]}\ACESII_36359_flight_11202022_fixedtm_20221120T170834.cdf"

                        with pycdf.CDF(path_to_processed_scotts_file) as scottsFile:
                            epoch_from_file = scottsFile['Epoch'][...]
                            tmCDF_epoch_scott = np.array([pycdf.lib.datetime_to_tt2000(epoch_from_file[i]) for i in range(len(epoch_from_file))])
                            tmCDF_sfid_scott = np.array(scottsFile['sfid'][...])
                            tmCDF_mf_scott = np.array(scottsFile['minorframe'][...])


                        a0 = 722237393164357000 # epoch where sfid 60 occurs in BSS
                        index60 = np.abs(tmCDF_epoch_scott - a0).argmin()

                        # replace sfid 60 point
                        spindex = 4725425  # where sfid 60 occurs in BSS
                        tmCDF_sfid[spindex] = tmCDF_sfid_scott[index60]
                        tmCDF_mf[spindex] = tmCDF_mf_scott[index60]
                        tmCDF_epoch[spindex] = tmCDF_epoch_scott[index60]

                        # --- PROBLEM #2 insert data for the missing 12 values ---
                        # The 12 values look like 15 16 17 [NONE] 30 31 32
                        spindex = 3653578

                        # remove the bad value
                        tmCDF_sfid = np.delete(tmCDF_sfid,[spindex],axis=0)
                        tmCDF_mf = np.delete(tmCDF_mf, [spindex], axis=0)
                        tmCDF_epoch = np.delete(tmCDF_epoch, [spindex], axis=0)

                        # insert the new value
                        tmCDF_sfid = np.insert(tmCDF_sfid, spindex, [18 + i for i in range(11+1)] )
                        tmCDF_mf = np.insert(tmCDF_mf, [spindex + i for i in range(11+1)], [[-1 for i in range(150)] for j in range(11+1)],axis= 0 )
                        tmCDF_epoch = np.insert(tmCDF_epoch, spindex, [rocketAttrs.epoch_fillVal for i in range(11+1)] )

                        Done(start_time)

                    # --- Fix BSS Low Flyer ---
                    if plugHolesInBSSData_Low:
                        prgMsg('Plugging Holes in Data')

                        fix_indicies = []
                        counter = tmCDF_sfid[0]
                        for index in tqdm(range(len(tmCDF_sfid))):

                            if tmCDF_sfid[index] != counter:
                                fix_indicies.append([tmCDF_sfid[index], counter])

                                if tmCDF_sfid[index] == -1:  # fabricate the gSync'd Epochs. No worry about starting with a -1 b/c of the code above
                                    tmCDF_epoch[index] = tmCDF_epoch[index - 1] + deltaT_median
                                    tmCDF_sfid[index] = counter
                                    tmCDF_mf[index] = [-1 for i in range(150)]

                            counter += 1
                            if counter >= 40:
                                counter = 0

                        Done(start_time)

                    # --- --- --- --- --- --- --- ---
                    # --- OUTPUT PROCESSED FILE ---
                    # --- --- --- --- --- --- --- ---

                    # Write out the output data in chunks
                    prgMsg('Writing out the processed data')
                    tmCDF_epoch_split = np.array(np.array_split(tmCDF_epoch, chunkNum))
                    tmCDF_sfid_split = np.array(np.array_split(tmCDF_sfid, chunkNum))
                    tmCDF_mf_split = np.array(np.array_split(tmCDF_mf, chunkNum))

                    for i in tqdm(range(chunkNum)):
                        if i == 0: # Overwrite the initialized data
                            outputFile['Epoch'][0:] = tmCDF_epoch_split[i]
                            outputFile['sfid'][0:] = tmCDF_sfid_split[i]
                            outputFile['minorframe'][0:] = tmCDF_mf_split[i]
                        else:
                            outputFile['Epoch'][len(outputFile['Epoch']):] = tmCDF_epoch_split[i]
                            outputFile['sfid'][len(outputFile['sfid']):] = tmCDF_sfid_split[i]
                            outputFile['minorframe'][len(outputFile['minorframe']):] = tmCDF_mf_split[i]

                    if len(aEpoch) != 0:
                        epochSet = sorted(set(tmCDF_epoch))
                        if len(epochSet) > 1:
                            outputFile['Epoch'].attrs['VALIDMIN'] = epochSet[1]
                            outputFile['Epoch'].attrs['VALIDMAX'] = epochSet[-1]

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

if len(glob(rf'{rocketFolderPath}tad\{fliers[wflyer]}\*.tad')) == 0:
    print(color.RED + 'There are no .tad files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        tadToTelem(wRocket, 0, chunkNum, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}tad\*.tad')))):
            tadToTelem(wRocket, fileNo, chunkNum, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            tadToTelem(wRocket, filesNo, chunkNum, rocketFolderPath, justPrintFileNames,wflyer)