# --- L1_magCal_to_L1_ChiSquare_cal2.py ---
# --- Author: C. Feltman ---
# CODE DESCRIPTION: performs the ChiSqaure padpair calibration on the integration_tad_files-aligned
# ESA data. This constitutes the FINAL calibration needed for the ESA and so the
# output file will be put into the "L1" folder as "ACESII_XXXXX_l1_instrNam_fullCal.cdf

# ChiCal DESCRIPTION: use a single parameter fit between uncalibrated pads and
# principal pads, excluding any data that has a 0 count value in principal or uncal data.

# NOTE: The Low Flyer did not have any calibration data worthy of use. The statistics were never larger than 50 pairs
# and the data itself is nearly vertical


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

# --- --- --- --- ---
from src.my_imports import *
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files
justPrintFileNames = False
justPrintChiFileNames = False

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [[0, 1, 2], [0, 1]]

inputPath_modifier = 'calibration\ESA_magPitch_calibration' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_chiSquare = 'calibration\ESA_ChiSquare_calibration' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'L1' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
modifier = ''

# --- CHI SQUARE CALIBRATION TOGGLES ---
viewPadPairData_numOfPairs = False
plotPadPairData = False
outlierThresh = 30 # (nominally 30) determine threshold to remove datapoints that are considered outliers
ChiTheshold = [0.42, 2] # bounds that ChiSquare must be to be allowed to adjust the data. Inclusive

# order to ChiSquare fits [uncalpal,prinpal]. NOTE: ONLY FOR HIGH FLYER since low flyer was too spin-aligned
ChiOrder36359_EEPAA = [[10, 20], [0, 10], [-10, 20], [30, 20], [40, 30], [50, 40], [60, 50], [70, 60], [80, 90], [100, 90], [180,170],[160,170],[190,160],[150,160],[150,160],[140,150],[130,140],[120,130],[110,120]]
ChiOrder36359_LEESA = [[0, 10], [20, 10], [-10, 20], [150, 160], [170, 160], [180, 170], [190, 160]]

# --- Output Data ---
outputData = True



# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import spaceToolsLib as stl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def L1_magCalESA_to_L1_ChiSquareCald_cal2(wRocket, wFile, rocketFolderPath, justPrintFileNames):

    # --- ACES II Flight/Integration Data ---
    rocketID = ACESII.payload_IDs[wRocket-4]
    globalAttrsMod = ACESII.global_attributes[wRocket-4]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L1'

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}{modifier}\*magCal*')
    inputFiles_chiSquare = glob(f'{rocketFolderPath}{inputPath_modifier_chiSquare}\{ACESII.fliers[wRocket-4]}{modifier}\*.cdf')

    outputFiles = glob(f'{rocketFolderPath}{outputPath_modifier}\{ACESII.fliers[wRocket-4]}\*_fullcalibration.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}{modifier}\\', '') for ifile in inputFiles]
    input_names_chiSquare = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier_chiSquare}\{ACESII.fliers[wRocket-4]}{modifier}\\', '') for ifile in inputFiles_chiSquare]

    output_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{ACESII.fliers[wRocket-4]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    input_names_searchable_chiSquare = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names_chiSquare]

    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}{modifier}\\','')

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in dataFile_name:
            wInstr = [index, instr]

    fileoutName = f'ACESII_{rocketID}_l1_{wInstr[1]}_fullCal.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return
    elif justPrintChiFileNames:
        for i, file in enumerate(inputFiles_chiSquare):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable_chiSquare[i], round(getsize(file) / (10 ** 6), 1)))
        return



    print('\n')
    print(stl.color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + stl.color.END)
    print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

    # --- get the data from the l1 ESA file ---
    stl.prgMsg(f'Loading data from {inputPath_modifier} Files')
    data_dict_esa = stl.loadDictFromFile(inputFiles[wFile])
    data_dict_esa['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_esa['Epoch'][0][i]) for i in range(len(data_dict_esa['Epoch'][0]))])
    stl.Done(start_time)

    # --- get the data from the MagPitch file ---
    stl.prgMsg(f'Loading data from ChiSquare Files')

    this_chiSquareFile = ''
    for file in inputFiles_chiSquare:
        if wInstr[1] in file:
            this_chiSquareFile = file

    data_dict_chiSquare = stl.loadDictFromFile(this_chiSquareFile)

    stl.Done(start_time)

    # --- --- --- --- --- --- --- -
    # --- PERFORM CHISQUARE CAL ---
    # --- --- --- --- --- --- --- -

    stl.prgMsg('Preparing ChiSquare Data')

    # --- viewing the calibration data ---
    # description: organize the calibration data into bins that show all
    # the PRINCIPAL pads and how many datapoints fall into those pads. Then do the reverse

    # for each pad, create len(Pitch_Angle) - 1 empty bins to store all pairs from other bins that go into it.
    # This creates a len(Pitch_Angle)xlen(Pitch_Angle) matrix , where each element is a list full of pairs, each of format [uncalPoint,prinPoint]
    # these values are stored in "padPairs"

    padAngle = data_dict_esa['Pitch_Angle'][0]

    calData = list(data_dict_chiSquare['ChiSquare_calData'][0])

    padPairs = [[[] for i in range(len(padAngle))] for j in range(len(padAngle))]

    # fill in the padPairs matrix full of [uncal,prin] points. Here I have enforced that the plotted xaxis be "principal" pad and y-axis be "uncalibratedpad"

    # calData is stored in the format: [uncalData, prinData]

    for i in range(len(calData)):

        if not (int(calData[i][0]) == 0 or int(calData[i][1]) == 0): # Don't Include any points with a 0 count value
            padPairs[ int(calData[i][2]) ][  int(calData[i][3]) ].append([ int(calData[i][0]), int(calData[i][1])])

    # plot the padPair data
    if viewPadPairData_numOfPairs:
        vmins = [100, 10]
        vmaxes = [50000, 50000]

        viewData = [[] for i in range(len(padAngle))]

        for i in range(len(padPairs)): # first loop over len(Pitch_Angle)

            for j in range(len(padPairs[i])): # second loop over len(Pitch_Angle)
                viewData[i].append(len(padPairs[i][j]))

        # plot the number of pairs for each pad for each source

        if wInstr[1] != 'iepaa':
            axis = [-15 + 10*i for i in range(len(padAngle) + 1)]
            ticksteps = [-10, 200, 10]
        else:
            axis = [-15 + 30 * i for i in range(len(padAngle) + 1)]
            ticksteps = [0, 190, 30]

        X, Y = np.meshgrid(axis, axis)
        Z = viewData

        print('\n')
        for i in range(len(viewData)):
            prepstring = f'[{padAngle[i]}]'
            for thing in viewData[i]:
                prepstring = prepstring + f'{thing:8}'

            print(prepstring)
        print([f'[{pad}]' for pad in padAngle])

        fig, ax = plt.subplots()
        cmap = ax.pcolormesh(X, Y, Z, vmin=vmins[wRocket-4], vmax=vmaxes[wRocket-4], cmap='turbo',norm='log')
        ax.set_xticks(range(*ticksteps))
        ax.set_yticks(range(*ticksteps))
        ax.set_ylabel('Uncalibrated Pad')
        ax.set_xlabel('Principal Pad')
        plt.title(f'ACES-II {rocketID} {wInstr[1]}')
        colorbar = plt.colorbar(mappable=cmap,label='# of pairs')
        plt.grid()
        plt.show()
    stl.Done(start_time)

    # --- --- --- --- --- --- ---
    # --- FIT USING CHISQUARE ---
    # --- --- --- --- --- --- ---
    if wRocket == 4 and wInstr[1] != 'iepaa':
        stl.prgMsg('Correcting Data using ChiSquare')

        ChiOrder = ChiOrder36359_EEPAA if wInstr[1] == 'eepaa' else ChiOrder36359_LEESA
        esaData_fullCal = np.array(data_dict_esa['counts'][0])

        # --- loop through the pad-pairs ---
        counter = 0
        for padPair in ChiOrder:
            counter += 1
            uncalPad = np.abs(padAngle - padPair[0]).argmin()
            prinPad = np.abs(padAngle - padPair[1]).argmin()
            dataChiCal = padPairs[uncalPad][prinPad]
            uncalData = np.array([dataChiCal[i][0] for i in range(len(dataChiCal))])
            prinData = np.array([dataChiCal[i][1] for i in range(len(dataChiCal))])

            # remove outliers from the data
            uncalData_temp = []
            prinData_temp = []

            for i in range(len(uncalData)):

                if np.abs(uncalData[i] - prinData[i]) < outlierThresh: # remove big outliers

                    if not (uncalData[i] >= 65535 or prinData[i] >= 65535):

                        if (uncalData[i] > 0 and prinData[i] > 0): # remove saturation values and any negative values
                            uncalData_temp.append(uncalData[i])
                            prinData_temp.append(prinData[i])

            uncalData, prinData = np.array(uncalData_temp), np.array(prinData_temp)

            # if the # of padPairs aren't zero after removing outliers
            if not (len(uncalData) < 2 or len(prinData) < 2):
                def calFunc(x, A):
                    return A * x
                params, cov = curve_fit(calFunc, uncalData, prinData)
                uncalData_fitted = np.array([int(round(calFunc(x, params[0]))) for x in uncalData])
                uncalData_errors, prinData_errors = uncalData, prinData

                nu = len(prinData) - len(params)
                chiSquare = stl.calcChiSquare(prinData, uncalData_fitted, prinData_errors, uncalData_errors, nu)

                # --- --- --- --- --- --- --- -
                # --- APPLY FIT TO ESA DATA ---
                # --- --- --- --- --- --- --- -
                if chiSquare >= ChiTheshold[0] and chiSquare <= ChiTheshold[1]:

                    for tme in range(len(esaData_fullCal)):
                        for engy in range(len(esaData_fullCal[0][0])):

                            if esaData_fullCal[tme][uncalPad][engy] >= 0:

                                if esaData_fullCal[tme][uncalPad][engy] == 0:
                                    calibratedValue = 0
                                elif esaData_fullCal[tme][uncalPad][engy] <= 3:
                                    if params[0] <= 2/3:
                                        calibratedValue = 0
                                    else:
                                        calibratedValue = 3
                                else:
                                    calibratedValue = round(esaData_fullCal[tme][uncalPad][engy] * params[0])

                                esaData_fullCal[tme][uncalPad][engy] = calibratedValue

                            else:
                                esaData_fullCal[tme][uncalPad][engy] = ACESII.epoch_fillVal


                if plotPadPairData:

                    fig, ax = plt.subplots()
                    ax.scatter(uncalData_fitted,prinData)
                    ax.set_ylabel('Principal Data [counts]')
                    ax.set_xlabel('Uncalibrated Data Fitted [counts]')
                    ax.set_title(f'{counter}\n'
                                 f'Uncalibrated Pad: {padAngle[uncalPad]}$^\circ$\n'
                                 f'Principal Pad: {padAngle[prinPad]}$^\circ$ \n $\chi$ {chiSquare}\n'
                                 f'No. of Points: {len(prinData)}')

                    # --- Fit the calibration curve ---
                    N = 100
                    xData_fitted = np.linspace(uncalData.min(), prinData.max(), N)
                    yData_fitted = np.array([calFunc(x, params[0]) for x in xData_fitted])

                    ax.plot(xData_fitted, yData_fitted,color='red')
                    plt.legend([f'A: {round(params[0],3)}\n'])
                    plt.tight_layout()

                    plt.show()

        stl.Done(start_time)
    else:
        esaData_fullCal = np.array(data_dict_esa['counts'][0])


    # --- --- --- --- --- --- --- ---
    # --- PREPARE DATA FOR OUTPUT ---
    # --- --- --- --- --- --- --- ---
    if outputData:
        stl.prgMsg('Outputting Data')
        data_dict_esa['counts'][0] = esaData_fullCal
        data_dict_esa['counts'][1]['VALIDMIN'] = esaData_fullCal.min()
        data_dict_esa['counts'][1]['VALIDMAX'] = esaData_fullCal.max()

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        globalAttrsMod['Descriptor'] = ACESII.ESA_names[wInstr[0]]
        stl.outputCDFdata(outputPath, data_dict_esa, globalAttrsMod= globalAttrsMod,instrNam= wInstr[1])

        stl.Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = DataPaths.ACES_data_folder

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}\*.cdf')) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    if justPrintFileNames:
        L1_magCalESA_to_L1_ChiSquareCald_cal2(wRocket, 0, rocketFolderPath, justPrintFileNames)
    elif not wFiles[wRocket-4]:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}\*.cdf')))):
            L1_magCalESA_to_L1_ChiSquareCald_cal2(wRocket, fileNo, rocketFolderPath, justPrintFileNames)
    else:
        for filesNo in wFiles[wRocket-4]:
            L1_magCalESA_to_L1_ChiSquareCald_cal2(wRocket, filesNo, rocketFolderPath, justPrintFileNames)