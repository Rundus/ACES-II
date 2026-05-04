# --- RingCore_L1_DetermineAdjustments_before_despin.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: interpolate attitude data (alt,lat,long,DCM) onto ringCore Epoch and use
# it to Despin the ringCore data. The interpolation is done by drawing a line between
# successive points in the attitude data, seeing how many integration_tad_files epoch values fall
# between the two attitude epoch values, then evaluating the line at the magnetomter epoch values
# and storing the new attitude data.


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
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
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4
wFiles = [0]

inputPath_modifier = 'l1'
inputPath_modifier_attitude = 'attitude' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier_IGRF = 'science\IGRF_interpolated' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
outputPath_modifier_L2 = 'L2' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder


# data is reduced to this region
targetTimes = [[dt.datetime(2022,11,20,17,24,55,000000), dt.datetime(2022,11,20,17,25,11,000000)],
               [dt.datetime(2022,11,20,17,24,53,000000), dt.datetime(2022,11,20,17,25,10,000000)]]# reduce data l2 to only regions of science interest


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# --- constant offset over whole flight ---
SECTION_determineOffsets = True
# interactive plot
interactivePlot = True

# grid search

# time scalar offset + b
Nn = 1
deltaB = 0
offsetGridSearchRange_b = [104944444.44444445 - deltaB,104944444.44444445 + deltaB] # + offsets applied to attitude solution timeseries. Only use positive offsets since ACS was always faster than MAG. Format: [linspace LOW, linspace HIGH]

# time slope m*t
M = 1
deltaT = 0
offsetGridSearchRange_m = [1-deltaT, 1.0 + deltaT] # + offsets applied to attitude solution timeseries. Only use positive offsets since ACS was always faster than MAG. Format: [linspace LOW, linspace HIGH]

# vector amplitude A*
K = 1
deltaA = 0
offsetGridSearchRange_A = [1, 1]

# vector +offset
L = 20
wVecComp = 0 # which component of the RKT data to look at
deltaVecComp = 382
offsetGridSearchRange_VecComp = [-deltaVecComp, -375]

# --- results of determineScalarOffset ---
offsetResults_intercept = [104944444.44444445, 117000000]
offsetResults_slope = [1, 1]
offsetResults_Amplitude = [1, 1]
offsetResults_VectorScalar = [np.array([-234.47368421052633,0,0]),np.array([0,0,0])]

# HF results: [104944444.44444445, 390227194712994.5] format: [time in ns, BestChiSquare]
# LF results: [105789473.68421052, 1715649142403403.8]

# --- determine the grid search parameters to remove coning using FFT power ---
SECTION_grid_Search_coning = False
plotInteractiveDeltaB = True
plotFFT = False
wAxis = 'Bx'


outputData = True
# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from pyIGRF import igrf_value
from numpy import array,matmul,dot,abs
from numpy.linalg import inv,norm
from myspaceToolsLib.filter import butter_filter
from scipy.interpolate import CubicSpline
rocketAttrs, b, c = ACES_mission_dicts()
from scipy.fft import rfft, fftfreq
from matplotlib.widgets import Slider
def buildDataDict(labels, dataToOutput):
    data_dict = {}
    for i, data in enumerate(dataToOutput):
        if labels[i] != 'Epoch':  # non-epoch case
            data_dict = {**data_dict, **{labels[i]: [data, {'LABLAXIS': labels[i], 'DEPEND_0': 'Epoch',
                                                            'DEPEND_1': None,
                                                            'DEPEND_2': None,
                                                            'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                            'UNITS': 'nT',
                                                            'VALIDMIN': data.min(),
                                                            'VALIDMAX': data.max(),
                                                            'VAR_TYPE': 'data',
                                                            'SCALETYP': 'linear'}]}}
        else:  # Epoch case

            # data = array([datum + offset[wRocket-4] for datum in data])
            data_dict = {**data_dict, **{labels[i]: [data, {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None,
                                                            'FILLVAL': rocketAttrs.epoch_fillVal,
                                                            'FORMAT': 'I5', 'UNITS': 'ns',
                                                            'VALIDMIN': data.min(), 'VALIDMAX': data.max(),
                                                            'VAR_TYPE': 'support_data',
                                                            'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000',
                                                            'TIME_SCALE': 'Terrestrial Time',
                                                            'REFERENCE_POSITION': 'Rotating Earth Geoid',
                                                            'SCALETYP': 'linear'}]}}
    return data_dict

def RingCore_Despin(wRocket,wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    ModelData = L0_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*RingCore*')
    inputFiles_attitude = glob(f'{rocketFolderPath}{inputPath_modifier_attitude}\{fliers[wflyer]}\*.cdf')
    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\\', '') for ifile in inputFiles]
    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]


    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Interpolating L1 data to DeSpun L2 data' + color.END)
        print('[0] ' + str(round(getsize(inputFiles[0]) / (10 ** 6), 1)) + 'MiB')

        ###########################
        # --- Get the Variables ---
        ###########################

        # --- get the Magnetometer Data ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict_mag = loadDictFromFile(inputFiles[wFile], {})
        data_dict_mag['Epoch'][0] = array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_mag['Epoch'][0]])
        Done(start_time)

        # --- get the Attitude Data ---
        prgMsg(f'Loading data from {inputPath_modifier_attitude} Files')
        data_dict_attitude = loadDictFromFile(inputFiles_attitude[0], {})
        data_dict_attitude['Epoch'][0] = array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_attitude['Epoch'][0]])
        Done(start_time)

        # --- reduce dataset to science region ---
        prgMsg('Reducing data to science region')
        tTimes = [pycdf.lib.datetime_to_tt2000(targetTimes[wRocket - 4][0]),
                  pycdf.lib.datetime_to_tt2000(targetTimes[wRocket - 4][1])]

        lowCutoff, highCutoff = np.abs(data_dict_attitude['Epoch'][0] - tTimes[0]).argmin(), np.abs(data_dict_attitude['Epoch'][0] - tTimes[1]).argmin()
        for key, val in data_dict_attitude.items(): # apply reduction to attitude data
            if key == 'Alt':
                data_dict_attitude[key][0] = array(data_dict_attitude[key][0][lowCutoff:highCutoff])/1000 # convert to km
            else:
                data_dict_attitude[key][0] = array(data_dict_attitude[key][0][lowCutoff:highCutoff])

        lowCutoff, highCutoff = np.abs(data_dict_mag['Epoch'][0] - tTimes[0]).argmin(), np.abs(data_dict_mag['Epoch'][0] - tTimes[1]).argmin()
        for key, val in data_dict_mag.items(): # apply reduction to mag data
            data_dict_mag[key][0] = array(data_dict_mag[key][0][lowCutoff:highCutoff])

        Done(start_time)

        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ### INTERPOLATE ATTITUDE DATA ###
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        # --- [0] Collect Variables ---
        prgMsg('Collecting Variables')

        # define some storage variables
        dataKeys = ['Epoch', 'Alt', 'Latgd', 'Long', 'a11', 'a12', 'a13', 'a21', 'a22', 'a23', 'a31', 'a32', 'a33']
        dataKeysVal = [deepcopy(data_dict_mag['Epoch'][0]), [], [], [], [], [], [], [], [], [], [], [], []]
        attitudeData = [deepcopy(data_dict_attitude[key][0]) for key in dataKeys] # a list to contain the attitude only the data that I care about
        dataInterp_dict_attitude = {key : value for key, value in zip(dataKeys, dataKeysVal)}

        # Define the Rocket Magnetic Field Vector
        B_rkt = array([[data_dict_mag['Bx'][0][i], data_dict_mag['By'][0][i], data_dict_mag['Bz'][0][i]] for i in range(len(data_dict_mag['Epoch'][0]))])

        Done(start_time)

        if SECTION_determineOffsets:

            ############################
            # --- FIND SCALAR OFFSET ---
            ############################
            prgMsg('Determining best scalar offset')
            offsetsVals_b = np.linspace(offsetGridSearchRange_b[0], offsetGridSearchRange_b[1], Nn)
            offsetsVals_m = np.linspace(offsetGridSearchRange_m[0], offsetGridSearchRange_m[1], M)
            offsetsVals_A = np.linspace(offsetGridSearchRange_A[0], offsetGridSearchRange_A[1], K)
            offsetsVals_VecComp = np.linspace(offsetGridSearchRange_VecComp[0], offsetGridSearchRange_VecComp[1], L)
            bestOffset = [0, 1, 1E30] # format [offsetVal (in secs), ChiSquareVal]. Dummy values to be replaced

            if interactivePlot:

                #############################################
                # --- try a slider method ---
                Epoch_mag = array([ pycdf.lib.tt2000_to_datetime(tme) for tme in data_dict_mag['Epoch'][0]])

                # every time I call this, I need a new timebase and a new set of IGRF values
                def f(slope, intercept):
                    # --- [1] apply offset to attitude data ---
                    Epoch_attitude_loop = array([slope * attitudeData[0][i] + intercept for i in range(len(attitudeData[0]))])

                    # --- [2] Spline interpolate the time-adjusted attitude data onto the mag data ---
                    counter = 0
                    for key, newDataList in dataInterp_dict_attitude.items():
                        if key != 'Epoch':
                            # --- cubic interpolation ---
                            splCub = CubicSpline(Epoch_attitude_loop, attitudeData[counter])

                            # evaluate the interpolation at all the epoch_mag points
                            dataInterp_dict_attitude[key] = array(
                                [splCub(timeVal) for timeVal in dataInterp_dict_attitude['Epoch']])
                        counter += 1

                    # EVALUATE SPUN-UP IGRF MODEL

                    date = 2022 + 323 / 365  # Corresponds to 11/20/2022

                    # --- [3] Evaluate IGRF at new attitude coordinates ---
                    ### IGRF info ###
                    # [3] North Comp (+ N | - S)
                    # [4] East Comp (+ E | - W)
                    # [5] Vertical Comp (+ D | - U)
                    # [6] Total Field
                    IGRF = array([igrf_value(dataInterp_dict_attitude['Latgd'][i], dataInterp_dict_attitude['Long'][i],
                                             dataInterp_dict_attitude['Alt'][i], date) for i in
                                  range(len(dataInterp_dict_attitude['Epoch']))])
                    IGRF_ENU = array([[vec[4], vec[3], -1 * vec[5]] for vec in IGRF])

                    # --- [4] Despin the mag data ---
                    DCM = array([
                        [[dataInterp_dict_attitude['a11'][i], dataInterp_dict_attitude['a12'][i],
                          dataInterp_dict_attitude['a13'][i]],
                         [dataInterp_dict_attitude['a21'][i], dataInterp_dict_attitude['a22'][i],
                          dataInterp_dict_attitude['a23'][i]],
                         [dataInterp_dict_attitude['a31'][i], dataInterp_dict_attitude['a32'][i],
                          dataInterp_dict_attitude['a33'][i]]]
                        for i in range(len(dataInterp_dict_attitude['Epoch']))
                    ])  # construct the new, interpolated DCM matrix

                    # --- [5] get inverse of DCM, spin up IGRF ---
                    IGRF_spun = array(
                        [matmul(inv(DCM[i]), IGRF_ENU[i]) for i in range(len(dataInterp_dict_attitude['Epoch']))])
                    axes = array(['Bx','By','Bz'])
                    thisOne = np.where(axes==wAxis)[0][0]

                    return IGRF_spun[:,thisOne]

                index = [i for i in range(len(attitudeData[0]))]

                # Define initial parameters
                init_slope = (offsetGridSearchRange_m[1] + offsetGridSearchRange_m[0])/2
                init_intercept = (offsetGridSearchRange_b[1] + offsetGridSearchRange_b[0])/2
                init_Amplitude = (offsetGridSearchRange_A[1] + offsetGridSearchRange_A[0]) / 2

                # Create the figure and the line that we will manipulate
                fig, ax = plt.subplots(2)
                lineRKT, = ax[0].plot(Epoch_mag, data_dict_mag[wAxis][0], lw=2, color='black') # initialize rkt data
                line, = ax[0].plot(Epoch_mag, f(init_slope, init_intercept), lw=2,color='blue') # initialize the spun IGRF

                RKTdata = array(data_dict_mag[wAxis][0])
                init_sub = (init_Amplitude + RKTdata) - f(init_slope, init_intercept)
                lineSub, = ax[1].plot(Epoch_mag,init_sub,lw=2,color='orange')
                ax[0].set_ylabel(wAxis)
                ax[1].set_ylabel(f'RKT - {wAxis}_IGRF')
                ax[1].set_ylim(-200, 200)
                ax[1].set_xlabel('Epoch')

                # adjust the main plot to make room for the sliders
                fig.subplots_adjust(left=0.25, bottom=0.25)

                # Make a horizontal slider to control the slope.
                axSlope = fig.add_axes([0.25, 0.1, 0.65, 0.03])
                slope_slider = Slider(
                    ax=axSlope,
                    label='Slope',
                    valmin=offsetGridSearchRange_m[0],
                    valmax=offsetGridSearchRange_m[1],
                    valinit=init_slope,
                )

                # Make a horizontal slider to control the slope.
                axAmp = fig.add_axes([0.25, 0.15, 0.65, 0.03])
                amplitude_slider = Slider(
                    ax=axAmp,
                    label='Amplitude',
                    valmin=offsetGridSearchRange_A[0],
                    valmax=offsetGridSearchRange_A[1],
                    valinit=init_Amplitude,
                )

                # Make a vertically oriented slider to control the amplitude
                axIntercept = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
                intercept_slider = Slider(
                    ax=axIntercept,
                    label="Intercept",
                    valmin=offsetGridSearchRange_b[0],
                    valmax=offsetGridSearchRange_b[1],
                    valinit=init_intercept,
                    orientation="vertical"
                )
                fig.suptitle(f'Slope : {init_slope} \n'
                             f'Intercept : {init_intercept} \n'
                             f'Amplitude : {init_Amplitude} \n')

                # The function to be called anytime a slider's value changes
                def update(val):
                    newData = f(slope_slider.val, intercept_slider.val)
                    line.set_ydata(newData)
                    lineRKT.set_ydata(amplitude_slider.val + RKTdata)
                    lineSub.set_ydata((amplitude_slider.val + RKTdata)-newData)
                    fig.suptitle(f'Slope : {slope_slider.val} \n'
                             f'Intercept : {intercept_slider.val} \n'
                                 f'Amplitude : {amplitude_slider.val}')

                    fig.canvas.draw_idle()


                # register the update function with each slider
                intercept_slider.on_changed(update)
                slope_slider.on_changed(update)
                amplitude_slider.on_changed(update)

                plt.show()
            else:
                # --- Loop over the temporal offsets ---
                for loopIndex_b in tqdm(range(len(offsetsVals_b))):
                    for loopIndex_m in range(len(offsetsVals_m)):
                        for loopIndex_A in range(len(offsetsVals_A)):
                            for loopIndex_VecComp in range(len(offsetsVals_VecComp)):


                                # --- [1] apply offset to attitude data ---
                                Epoch_attitude_loop = array([offsetsVals_m[loopIndex_m]*attitudeData[0][i] + offsetsVals_b[loopIndex_b] for i in range(len(attitudeData[0]))])

                                # --- [2] Spline interpolate the time-adjusted attitude data onto the mag data ---
                                counter = 0
                                for key, newDataList in dataInterp_dict_attitude.items():
                                    if key != 'Epoch':
                                        # --- cubic interpolation ---
                                        splCub = CubicSpline(Epoch_attitude_loop, attitudeData[counter])

                                        # evaluate the interpolation at all the epoch_mag points
                                        dataInterp_dict_attitude[key] = array([splCub(timeVal) for timeVal in dataInterp_dict_attitude['Epoch']])
                                    counter += 1

                                # EVALUATE SPUN-UP IGRF MODEL

                                date = 2022 + 323 / 365  # Corresponds to 11/20/2022

                                # --- [3] Evaluate IGRF at new attitude coordinates ---
                                ### IGRF info ###
                                # [3] North Comp (+ N | - S)
                                # [4] East Comp (+ E | - W)
                                # [5] Vertical Comp (+ D | - U)
                                # [6] Total Field
                                IGRF = array([igrf_value(dataInterp_dict_attitude['Latgd'][i], dataInterp_dict_attitude['Long'][i], dataInterp_dict_attitude['Alt'][i], date) for i in range(len(dataInterp_dict_attitude['Epoch']))])
                                IGRF_ENU = array([ [vec[4], vec[3], -1 * vec[5]] for vec in IGRF])

                                # --- [4] Despin the mag data ---
                                DCM = array([
                                    [[dataInterp_dict_attitude['a11'][i], dataInterp_dict_attitude['a12'][i], dataInterp_dict_attitude['a13'][i]],
                                     [dataInterp_dict_attitude['a21'][i], dataInterp_dict_attitude['a22'][i], dataInterp_dict_attitude['a23'][i]],
                                     [dataInterp_dict_attitude['a31'][i], dataInterp_dict_attitude['a32'][i], dataInterp_dict_attitude['a33'][i]]]
                                    for i in range(len(dataInterp_dict_attitude['Epoch']))
                                ]) # construct the new, interpolated DCM matrix

                                # --- [5] get inverse of DCM, spin up IGRF ---
                                IGRF_spun = array([matmul(inv(DCM[i]), IGRF_ENU[i]) for i in range(len(dataInterp_dict_attitude['Epoch']))])

                                # --- [6] Analyze results of loop to see if chosen offset is good ---
                                vecOffset = [0,0,0]
                                vecOffset[wVecComp] = offsetsVals_VecComp[loopIndex_VecComp]

                                sub = (offsetsVals_A[loopIndex_A]*B_rkt + vecOffset) - IGRF_spun
                                ChiSquare = (1 / len(dataInterp_dict_attitude['Epoch'])) * (sum([dot(sub[i], sub[i]) for i in range(len(sub))]) ** 2)
                                print([offsetsVals_m[loopIndex_m], offsetsVals_b[loopIndex_b], offsetsVals_A[loopIndex_A], offsetsVals_VecComp[loopIndex_VecComp],vecOffset, ChiSquare], '\n')
                                if ChiSquare <= bestOffset[-1]:
                                    bestOffset = [offsetsVals_m[loopIndex_m], offsetsVals_b[loopIndex_b], offsetsVals_A[loopIndex_A], offsetsVals_VecComp[loopIndex_VecComp],vecOffset, ChiSquare]

            Done(start_time)

            print(f'The Best Offset was: {bestOffset}')
        else: # just do a normal de-spin assuming the DCM is perfect

            Epoch_attitude_loop = array([offsetGridSearchRange_m[wRocket-4]*attitudeData[0][i] + offsetResults_intercept[wRocket - 4] for i in range(len(attitudeData[0]))])

            prgMsg('Interpolating Attitude Data')
            spline_dict_attitude = {}
            dataInterp_dict_attitude = {'a11': [], 'a12': [], 'a13': [], 'a21': [], 'a22': [], 'a23': [], 'a31': [],
                                        'a32': [], 'a33': [],'Alt':[], 'Latgd':[], 'Long':[]}
            for key, newDataList in dataInterp_dict_attitude.items():
                if 'Epoch'.lower() not in key.lower():
                    # --- cubic interpolation ---
                    splCub = CubicSpline(Epoch_attitude_loop, data_dict_attitude[key][0])

                    # store the spline information for later
                    spline_dict_attitude = {**spline_dict_attitude, **{key: splCub}}

                    # evaluate the interpolation at all the epoch_mag points
                    dataInterp_dict_attitude[key] = np.array([splCub(timeVal) for timeVal in data_dict_mag['Epoch'][0]])
            Done(start_time)

            #############################
            # --- Despin the Mag data ---
            #############################
            prgMsg('De-spinning Mag Data')
            # form the DCM matrix
            DCM = np.array([
                [[dataInterp_dict_attitude['a11'][i], dataInterp_dict_attitude['a12'][i], dataInterp_dict_attitude['a13'][i]],
                 [dataInterp_dict_attitude['a21'][i], dataInterp_dict_attitude['a22'][i], dataInterp_dict_attitude['a23'][i]],
                 [dataInterp_dict_attitude['a31'][i], dataInterp_dict_attitude['a32'][i], dataInterp_dict_attitude['a33'][i]]]
                for i in range(len(data_dict_mag['Epoch'][0]))
            ])

            invDCM = array([inv(DCM[i]) for i in range(len(data_dict_mag['Epoch'][0]))])

            # apply the despin
            B_rkt_despun = np.array([np.matmul(DCM[i], B_rkt[i]) for i in range(len(data_dict_mag['Epoch'][0]))])

            Done(start_time)


            # --- create a spun-up IGRF Field ---
            date = 2022 + 323 / 365  # Corresponds to 11/20/2022
            IGRF = array([igrf_value(dataInterp_dict_attitude['Latgd'][i], dataInterp_dict_attitude['Long'][i], dataInterp_dict_attitude['Alt'][i], date) for i in range(len(data_dict_mag['Epoch'][0]))])
            IGRF_ENU = array([[IGRF[i][4], IGRF[i][3], -1 * IGRF[i][5]] for i in range(len(IGRF))])
            IGRF_spun = offsetResults_Amplitude[wRocket-4]*array([ matmul(invDCM[i],IGRF_ENU[i]) for i in range(len(IGRF))])


        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:

            prgMsg('Creating output file')

            # items I need for diangostics
            # [1] rocket data that's shifted by the scalar amount
            # [2] spun IGRF field
            # [3] Attempted RKT Despun Data

            # # construct integration_tad_files data B-vector and despin the integration_tad_files data
            # B_ENU = matmul(array([data_dict_mag['Bx'][0][i], data_dict_mag['By'][0][i], data_dict_mag['Bz'][0][i]]), IGRF_ENU)

            IGRF_east, IGRF_north, IGRF_up = array(IGRF_ENU[:, 0]), array(IGRF_ENU[:, 1]), array(IGRF_ENU[:, 2])
            IGRF_rkt_X, IGRF_rkt_Y, IGRF_rkt_Z = array(IGRF_spun[:, 0]), array(IGRF_spun[:, 1]), array(IGRF_spun[:, 2])
            IGRF_mag = array([norm(IGRF_ENU[i]) for i in range(len(IGRF_ENU))])


            ########################
            ### IGRF OUTPUT DATA ###
            ########################

            # --- Construct the data dict for IGRF ---
            labels = ['Epoch', 'IGRF_east', 'IGRF_north', 'IGRF_up', 'IGRF_X', 'IGRF_Y', 'IGRF_Z', 'IGRF_mag']
            dataToOutput = [array(Epoch_attitude_loop), IGRF_east, IGRF_north, IGRF_up, IGRF_rkt_X, IGRF_rkt_Y, IGRF_rkt_Z, IGRF_mag]
            data_dict_IGRF = buildDataDict(labels, dataToOutput)

            fileoutName = f'ACESII_{rocketID}_IGRF_interpolated_model.cdf'
            outputPath = f'{rocketFolderPath}{outputPath_modifier_IGRF}\{fliers[wflyer]}\\{fileoutName}'
            outputCDFdata(outputPath, data_dict_IGRF, ModelData, globalAttrsMod, "RingCore")

            # -- Construct the data dict for despun mag data ---
            labels = ['Epoch', 'B_east', 'B_north', 'B_up', 'B_x', 'B_y', 'B_z']
            dataToOutput = [data_dict_mag['Epoch'][0], B_rkt_despun[:, 0], B_rkt_despun[:, 1], B_rkt_despun[:, 2], B_rkt[:, 0], B_rkt[:, 1], B_rkt[:, 2]]
            data_dict_Despun = buildDataDict(labels, dataToOutput)

            fileoutName = f'ACESII_{rocketID}_l2_RingCore_Despun.cdf'
            outputPath = f'{rocketFolderPath}{outputPath_modifier_L2}\{fliers[wflyer]}\\{fileoutName}'
            outputCDFdata(outputPath, data_dict_Despun, ModelData, globalAttrsMod, "RingCore")

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

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*_RingCore*')) == 0:
    print(color.RED + 'There are no RingCore.cdf files in the specified directory' + color.END)
else:
    for filesNo in wFiles:
        RingCore_Despin(wRocket, filesNo, rocketFolderPath, justPrintFileNames, wflyer)
