# --- MAG_L0_to_L1.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Magnetometer data to Science Units AND adjust Epoch-range with time-offset

# Truncates all data to only that past 17:20:00. Look at tmCDF_to_L0 to see when Mag data Epoch
# is sampled (should be sfid == 1 & 21)

# Added Correction: We forgot to take into account that the coil system at Iowa and at Wallops
# are left-handed to minimize gradients.  We simply flipped the integration_tad_files axis Y-Axis the
# wrong way rather than flip the reference data axis. So we add a minus sign to the mag Y-axis in
# our code.

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
wFiles = [1] # always set to [1]

modifier = ''
inputPath_modifier = 'L0' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'L1' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# calibration Toggles
applyScales = True
applyFifthOrderCal = True
applyRotateToRocket = True
applyFlipBy = False
applyVectorCals = True

# Fix outliers through interpolation
fixOutliers = True
percentThreshold = 10000 # change the value if the i+1 value has a percent difference of this percent
targetTimes = [[dt.datetime(2022,11,20,17,22,00,000000), dt.datetime(2022,11,20,17,27,00,000000)],
               [dt.datetime(2022,11,20,17,23,45,000000), dt.datetime(2022,11,20,17,27,55,000000)]]# the science region

outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from numpy.polynomial import Polynomial
from scipy.interpolate import CubicSpline

def MAG_L0_to_L1(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L1'
    outputModelData = L1_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*RingCore*')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]

    input_names_searchable = [ifile.replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('.cdf','') for ifile in input_names]

    fileoutName = f'ACESII_{rocketID}_{outputPath_modifier.lower()}_{input_names_searchable[wFile]}.cdf'

    fileoutName = fileoutName.replace('_magFrm', '_rktFrm') if applyRotateToRocket else fileoutName


    # determine which instrument the file corresponds to:
    mags = ['RingCore', 'Tesseract']
    for index, instr in enumerate(mags):
        if instr.lower() in inputFiles[wFile].lower():
            wInstr = [index, instr, mags[index]]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:35s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {input_names_searchable[wFile]}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the tmCDF file ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict_mag = loadDictFromFile(inputFiles[wFile], {})
        data_dict_mag['Epoch'][0] = [pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_mag['Epoch'][0]]
        Done(start_time)

        #######################
        # --- TRUNCATE DATA ---
        #######################
        prgMsg('Truncating Data')

        # --- do a quality check on the data to ensure epoch is > 0 or at least fillval ---
        # NOTE: There's some weirdness where the data doesn't recognize a negative epoch as a fillval. Below code
        # removes the weirdness
        badIndicies = []
        for i, tme in enumerate(data_dict_mag['Epoch'][0]):
            if tme < 0:
                badIndicies.append(i)

        for key, val in data_dict_mag.items():
            data_dict_mag[key][0] = np.delete(data_dict_mag[key][0], badIndicies)

        # --- remove data previous to 17:20:00 ---
        targetTime = pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 20, 00, 000000))
        truncIndex = np.abs(np.array(data_dict_mag['Epoch'][0]) - targetTime).argmin()
        for key, val in data_dict_mag.items():
            data_dict_mag[key][0] = data_dict_mag[key][0][truncIndex:]
        Done(start_time)

        B = np.array([[data_dict_mag['Bx'][0][i], data_dict_mag['By'][0][i], data_dict_mag['Bz'][0][i]] for i in range(len(data_dict_mag['Epoch'][0]))])

        # #################################
        # --- CONVERT TO SCIENCE UNITS ---
        # #################################
        prgMsg('Converting and Calibrating Data')

        if applyScales:
            scales = np.array(rocketAttrs.ringCoreScaleFactors[wflyer])  # Linear correction coefficents
            B = np.array([[data_dict_mag['Bx'][0][i]*scales[0], data_dict_mag['By'][0][i]*scales[1], data_dict_mag['Bz'][0][i]*scales[2]] for i in range(len(data_dict_mag['Epoch'][0]))])

            # --- assign data ---
            data_dict_mag['Bx'][0], data_dict_mag['By'][0], data_dict_mag['Bz'][0] = np.array(B[:, 0]), np.array(B[:, 1]), np.array(B[:, 2])


        if applyFifthOrderCal:
            fithOrderCorrection = np.array(rocketAttrs.ringCore5thOrderCorrections[wflyer]) # non-linear correction 5th, order coeffients. The terms descend in exponent order i.e. fithOrderCorrection[0] --> x^5 term, etc

            polys = [np.polynomial.Polynomial(fithOrderCorrection[0][::-1]),
                     np.polynomial.Polynomial(fithOrderCorrection[1][::-1]),
                     np.polynomial.Polynomial(fithOrderCorrection[2][::-1])]
            B = np.array([[polys[0](data_dict_mag['Bx'][0][i]), polys[1](data_dict_mag['By'][0][i]), polys[2](data_dict_mag['Bz'][0][i])] for i in range(len(data_dict_mag['Epoch'][0]))])

            # --- assign data ---
            data_dict_mag['Bx'][0], data_dict_mag['By'][0], data_dict_mag['Bz'][0] = np.array(B[:, 0]), np.array(B[:, 1]), np.array(B[:, 2])

        if applyRotateToRocket:
            # The integration_tad_files frame can be converted to rocket-body frame by a 90deg rotation about Z-axis
            from ACESII_code.class_var_func import Rz
            B = []

            for i in range(len(data_dict_mag['Epoch'][0])):
                B_temp = np.array([data_dict_mag['Bx'][0][i],data_dict_mag['By'][0][i],data_dict_mag['Bz'][0][i]])
                B_rot = np.matmul(Rz(-90), B_temp)
                B.append(B_rot)

            B = np.array(B)

            # --- assign data ---
            data_dict_mag['Bx'][0], data_dict_mag['By'][0], data_dict_mag['Bz'][0] = np.array(B[:, 0]), np.array(B[:, 1]), np.array(B[:, 2])


        if applyFlipBy:
            # The integration_tad_files Y-axis was flipped during calibration. Here we add back in the minus sign
            data_dict_mag['By'][0] = -1*np.array(data_dict_mag['By'][0])

        if applyVectorCals:
            # description: Apply the calibration matrix which is a vector cal with the model: orthogonality, offset, and Euler rotation.
            # It is a coupled matrix containing all cals found via robust weighted regression.
            # Columns 1:3 are multiplied by the spun data then column 4 is added (basically a line equation in 3-D)

            calMatrix = rocketAttrs.ringCoreCalMatrix[wRocket-4]
            orthoEuler = np.array([calMatrix[i][0:3] for i in range(len(calMatrix))])
            offset = np.array([calMatrix[i][3] for i in range(len(calMatrix))])
            B = np.array([np.matmul(orthoEuler, np.array([data_dict_mag['Bx'][0][i], data_dict_mag['By'][0][i], data_dict_mag['Bz'][0][i]])) + offset for i in range(len(data_dict_mag['Epoch'][0]))])

            # --- assign data ---
            data_dict_mag['Bx'][0], data_dict_mag['By'][0], data_dict_mag['Bz'][0] = np.array(B[:, 0]), np.array(B[:, 1]), np.array(B[:, 2])

        Done(start_time)

        if fixOutliers:
            # There's an outlier in the High Flyer Bx data. Here we'll identify it and fix it
            # NOTE: This assumes the outlier is alone! Not surounded by other outliers
            ttimes = targetTimes[wRocket-4]
            scienceRegionIndicies = [np.abs(data_dict_mag['Epoch'][0] -  pycdf.lib.datetime_to_tt2000(ttimes[0])).argmin(),
                                     np.abs(data_dict_mag['Epoch'][0] - pycdf.lib.datetime_to_tt2000(ttimes[1])).argmin()]

            labels = ['Bx', 'By', 'Bz']
            for label in labels:

                for i in range(len(data_dict_mag[label][0])-1):

                    if scienceRegionIndicies[0] <= i <= scienceRegionIndicies[1]:

                        percentChange = 100*(data_dict_mag[label][0][i+1] - data_dict_mag[label][0][i])/(np.abs(data_dict_mag[label][0][i]))


                        if np.abs(percentChange) > percentThreshold:

                            # spline interpolate around the problem point, which is the i+1th point
                            yData_fitpoints = [data_dict_mag[label][0][j] for j in [i-k for k in range(1,40)]]
                            xData_fitpoints = [data_dict_mag['Epoch'][0][j] for j in [i-k for k in range(1,40)]]
                            splCub = CubicSpline(xData_fitpoints[::-1], yData_fitpoints[::-1])
                            data_dict_mag[label][0][i] = splCub(data_dict_mag['Epoch'][0][i])

        # get the magnitude
        data_dict_mag['Bmag'][0] = np.array([np.linalg.norm([data_dict_mag['Bx'][0][i], data_dict_mag['By'][0][i], data_dict_mag['Bz'][0][i]]) for i in range(len(data_dict_mag['Epoch'][0]))])

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            if applyScales == False and applyVectorCals == False and applyFifthOrderCal == False:
                data_dict_mag['Bx'][1]['UNITS'] = 'ADC'
                data_dict_mag['By'][1]['UNITS'] = 'ADC'
                data_dict_mag['Bz'][1]['UNITS'] = 'ADC'
            else:
                data_dict_mag['Bx'][1]['UNITS'] = 'nT'
                data_dict_mag['By'][1]['UNITS'] = 'nT'
                data_dict_mag['Bz'][1]['UNITS'] = 'nT'

            prgMsg('Creating output file')

            for key, val in data_dict_mag.items():
                data_dict_mag[key][0] = np.array(data_dict_mag[key][0])

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'
            outputCDFdata(outputPath, data_dict_mag, outputModelData, globalAttrsMod, wInstr[1])
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
        MAG_L0_to_L1(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            MAG_L0_to_L1(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            MAG_L0_to_L1(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)