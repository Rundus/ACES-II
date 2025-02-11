#--- L1_to_L2_langmuir_fixed.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert the engineering Langmuir data to scientifically useful units. Also renames
# "Boom_Monitor_1/2" --> "Fixed_Boom_Monitor_1/2" etc

# it was discovered that ni_swept for the low AND high Flyer start with an extra value that should not be there.
# This is fixed in L0_to_L1.py

# It was discovered that the gain for the Swept LPs is probably too high. This means the instrument enters saturation
# very quickly. So both the real data and calibration data saturate too quickly


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

# --- --- --- --- ---

from src.my_imports import *

start_time = time.time()
# --- --- --- --- ---

#################
# --- TOGGLES ---
#################

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
useNanoAmps = True

###############################
# --- FIXED VOLTAGE TOGGLES ---
###############################
unit_conversion = 1e9 # 1 for Amps 10**9 for nano amps, etc
SECTION_fixedProbeCal = True
applyFixedCalCurve = False
plotFixedCalCurve = True

##################################
# --- SWEPT TO VOLTAGE TOGGLES ---
##################################
# Does a little analysis on the "step" variable to determine the average values of the steps.
# Needed to calculate from "Step" variable to voltage
# MUST BE ==TRUE WHEN SWEPT PROBE ==TRUE
SECTION_stepToVoltage = False

# the capacitive effect between the probe and plasma cause an RC decay on the data.
downSample_RCeffect = True # This toggle only uses the 10th (the last) datapoint of each voltage setpoint to eliminate this effect
keepThisManyPoints = 1 # how many points to keep when you downsample RC effect

#########################
# --- SWEPT PROBE CAL ---
#########################
SECTION_sweptProbe = False # use the swept Probe calibration data to convert analog values to current
removeADCValsOutofCalRange = True
applySweptCalCurve = True # False - data is in analog values, True - data is in current.
inputCalFile_Iowa = f'{ACES_data_folder}\calibration\LP_calibration\high\ACESII_LP_Cals_swept_Iowa.cdf'


# --- break into curves toggles ---
breakIntoCurves = True
targetVoltage_min = -1 # only care about voltage sweeps above this voltage value. Nominally -1
digitalVariance = 5 # how much the digitized step point can vary when looking for the top and bottom of curves. nominally = 5
indvEpochThresh = 15000000 # Value, in tt2000, that determines the time diff needed between epoch points to identify particular sweeps

#####################
# --- DATA OUTPUT ---
#####################
outputData = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---

import scipy
from warnings import filterwarnings # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.
filterwarnings("ignore")
from collections import Counter

def L1_to_Langmuir(wRocket, rocketFolderPath, justPrintFileNames, wflyer):

    outputFolderPath = rocketFolderPath + r'L2\\'

    # --- Get ACESII rocket Attributes ---
    rocketAttrs,b,c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'Langmuir'
    outputModelData = L2_TRICE_Quick(wflyer)

    # Set the paths for the file names
    L1Files = glob(f'{rocketFolderPath}L1\{fliers[wflyer]}\*_lp_*')
    LangmuirFiles = glob(f'{outputFolderPath}\{fliers[wflyer]}\*_langmuir_*')
    LangmuirSweptCalFiles = glob(f'{rocketFolderPath}\calibration\LP_calibration\{fliers[wflyer]}\*_345deg_*')
    L1_names = [ifile.replace(f'{rocketFolderPath}L1\{fliers[wflyer]}\\', '') for ifile in L1Files]
    L1_names_searchable = [ifile.replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace('l1_', '').replace('_v00', '') for ifile in L1_names]
    dataFile_name = L1_names_searchable[0].replace(f'{rocketFolderPath}L1\{fliers[wflyer]}\\', '').replace('lp_','').replace('ni_','').replace('ne_swept_','').replace('step_','').replace('ni_swept_','').replace('deltaNdivN_','')
    fileoutName = rf'ACESII_{rocketAttrs.rocketID[wRocket-4]}_l2_langmuir' + '.cdf'
    wInstr = [0, 'LangmuirProbe']

    if justPrintFileNames:
            for i, file in enumerate(L1Files):
                print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, L1_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Creating langmuir data for {fliers[wRocket-4]} flyer' + color.END)

        # --- get the data from the tmCDF file ---
        prgMsg('Loading data from L1Files')

        # Collect the LP data except deltaNdivN into a data dict
        data_dict = {}
        for file in L1Files:
            if 'deltaNdivN' not in file:
                with pycdf.CDF(file) as L1DataFile:
                    for key, val in L1DataFile.items():
                        if key not in data_dict:
                            data_dict = {**data_dict, **{key : [L1DataFile[key][...] , {key:val for key,val in L1DataFile[key].attrs.items()  }  ]  }  }

        # Collect the LP data except deltaNdivN into a data dict
        data_dict_cal = {}
        for file in LangmuirSweptCalFiles:
            if 'deltaNdivN' not in file:
                with pycdf.CDF(file) as LangmuirSweptCalFiles:
                    for key, val in LangmuirSweptCalFiles.items():
                        if key not in data_dict_cal:
                            data_dict_cal = {**data_dict_cal, **{key: [LangmuirSweptCalFiles[key][...], {key: val for key, val in LangmuirSweptCalFiles[key].attrs.items()}]}}
        Done(start_time)

        #####################
        # --- Fixed Probe ---
        #####################
        # description: Loads in the fixed calibrated data (provided by scott, found in missionAttributes.py).

        if SECTION_fixedProbeCal:
            prgMsg('Converting fixed probe to Voltage')
            def calFunction_fixed(x, A, B):
                y = A*x + B
                return y

            fixedCalResistances = rocketAttrs.LPFixed_calResistances[wRocket - 4]
            probeBias = rocketAttrs.LPFixedProbeBias[wRocket - 4]
            calibrationCurrents = []
            analog_vals = []

            # convert calibration data to current (in nanoamps)
            # NOTE: WE don't care about the "open case"
            for key, val in fixedCalResistances.items():
                if key != 'Open':
                    analog_vals.append(val)

                    if useNanoAmps:
                        calibrationCurrents.append(unit_conversion * probeBias / key)
                    else:
                        calibrationCurrents.append(probeBias / key)

            analog_vals, calibrationCurrents = np.array(analog_vals), np.array(calibrationCurrents)

            # apply a natural log scale to the data in order to prepare it for fitting
            calibrationCurrents = np.array([np.log(-1*cur) for cur in calibrationCurrents])

            # Fit a linear line to the log'd data
            parameters, covariance = scipy.optimize.curve_fit(calFunction_fixed, analog_vals, calibrationCurrents, maxfev=10000)

            if plotFixedCalCurve:
                figWidth =10
                figHeight =10
                Label_Fontsize = 15
                Title_Fontsize = 15

                import matplotlib
                matplotlib.rc('figure', figsize=(5, 5))
                plt.rcParams["figure.figsize"] = (figHeight, figWidth)
                xDataFit = np.array([i for i in range(1, 4096)])
                yDataFit = [calFunction_fixed(val, *parameters) for val in xDataFit]
                plt.plot(xDataFit, yDataFit, color='red')
                plt.scatter(analog_vals, calibrationCurrents)
                plt.xlabel('ADC Value',fontsize=Label_Fontsize)
                plt.ylabel(r'Ln($I_{cal}$) [nA]',fontsize=Label_Fontsize)
                plt.suptitle(f'FIXED LP - {rocketAttrs.rocketID[wRocket-4]}\n'
                             'Calculated Calibration Current vs Analog Value',fontsize=Title_Fontsize)
                plt.legend(['ln(y) = mx + b\n'
                           f'm: {parameters[0]}\n'
                           f'b: {parameters[1]}'])
                plt.savefig(rf'C:\Users\cfelt\Desktop\rockets\ACES-II\Instruments\fixedLangmuirProbe\Calibration\ACESII_{fliers[wflyer]}_CalCurve.png')

            # Apply the calibration function curve
            index = np.abs(data_dict['Epoch_ni'][0] - dt.datetime(2022,11,20,17,25,44,500000)).argmin()

            if applyFixedCalCurve:
                caldCurrent = np.array([
                    np.exp(calFunction_fixed(data_dict['ni'][0][i], parameters[0],parameters[1]))
                    for i in range(len(data_dict['ni'][0]))])
                fixedLPunits = 'nA'
            else:
                caldCurrent = np.array(data_dict['ni'][0])
                fixedLPunits = 'ADC'

            # apply a quality assurance step:
            if applyFixedCalCurve:
                if wRocket == 4:
                    for i in range(len(caldCurrent)):
                        if np.abs(caldCurrent[i]) > 2000:
                            caldCurrent[i] = rocketAttrs.epoch_fillVal
                elif wRocket == 5:
                    for i in range(len(caldCurrent)):
                        if np.abs(caldCurrent[i]) > 2000:
                            caldCurrent[i] = rocketAttrs.epoch_fillVal

            # --- --- --- --- --- --- ----
            # --- FIXED ERROR ANALYSIS ---
            # --- --- --- --- --- --- ----

            # (1) The error in the plasma density comes from the error in the fit coefficents for the conversion between I_probe and ADC
            # % error n_i = sqrt(a^2 * delta a^2 + delta b^2)

            # (2) to get delta a, delta b we can do a reduced chi-square analysis on the fit:
            # chi^2 = 1/nu SUM (f(x_i) - y_i)^2 / (sigmaX_i^2 + sigmaY_i^2)
            # here:
            # (i) f: a*ADC+ b
            # (ii) y_i: ith value of ln(I_probe)
            # (iii): sigmaX_i: error in ADC value
            # (iv): SigmaY_i: error in ln(I_probe)
            #
            # we can work out what delta a, delta b are from optimization calculus to get:
            # (i): delta a = Sxx/DELTA
            # (ii): delta b = S/DELTA
            # for S = SUM(1/sqrt(sigmaX_i^2 +sigmaY_i^2)), Sxx = SUM(x_i^2/sqrt(sigmaX_i^2 +sigmaY_i^2))

            # (3) we need delta ADC and delta Ln(I_probe).
            # (a) The best we can do is delta ADC = +/- 0.5 ADC
            # (b) for Ln(I_probe), we assume voltage error of +/- 0.01V and known resistance error or 1% we have"
            IprobePercentError = np.sqrt((0.01/5.05)**2 + 0.01**2) # = about 0.01
            # (c) we convert ^^^ to delta Ln(I_probe) by delta ln(I_probe)= \delta I_probe /Iprobe BUT delta I_probe ~ 0.01 * Iprobe / Iprobe ===> 0.01
            deltaLnProbe = IprobePercentError


            Sxx = sum([ analog_vals[i]**2 / np.sqrt((0.5**2) + deltaLnProbe**2) for i in range(len(analog_vals))])
            Sx = sum([analog_vals[i] / np.sqrt((0.5 ** 2) + deltaLnProbe**2) for i in range(len(analog_vals))])
            S = sum([1 / np.sqrt((0.5 ** 2) + deltaLnProbe**2) for i in range(len(analog_vals))])
            DELTA = S*Sxx - (Sx)**2
            delta_a = Sxx/DELTA
            delta_b = S/DELTA
            delta_n = np.array([np.sqrt(  (parameters[0]*delta_a)**2 + (delta_b)**2  )])

            data_dict = {**data_dict, **{'fixed_current': [caldCurrent, {'LABLAXIS': 'current',
                                                            'DEPEND_0': 'fixed_Epoch',
                                                            'DEPEND_1': None,
                                                            'DEPEND_2': None,
                                                            'FILLVAL': rocketAttrs.epoch_fillVal,
                                                            'FORMAT': 'E12.2',
                                                            'UNITS': fixedLPunits,
                                                            'VALIDMIN': caldCurrent.min(),
                                                            'VALIDMAX': caldCurrent.max(),
                                                            'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
            data_dict = {**data_dict, **{'ni_percent_error': [delta_n, {'LABLAXIS': 'ni_percent_error',
                                                                         'DEPEND_0': None,
                                                                         'DEPEND_1': None,
                                                                         'DEPEND_2': None,
                                                                         'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                         'FORMAT': 'E12.2',
                                                                         'UNITS': 'cm^-3',
                                                                         'VALIDMIN': 0,
                                                                         'VALIDMAX': 1,
                                                                         'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]}}
            data_dict['fixed_Epoch'] = data_dict.pop('Epoch_ni') # rename to fixed
            Done(start_time)


        #########################
        # --- STEP to Voltage ---
        #########################
        # description: The swept LP "step" variable contains analog values that of 100 steps (101 values in total) for values up and down
        # determining a nominal voltage to assign each of these values will indicate what voltage was applied
        # to the swept langmuir probes.

        if SECTION_stepToVoltage:

            prgMsg('Calculating swept step voltage')

            # calculate the epoch index of beginning/end of sample range
            Epoch_step = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_step'][0][i]) for i in range(len(data_dict['Epoch_step'][0]))])
            sampleStart = np.abs(np.array(Epoch_step - pycdf.lib.datetime_to_tt2000(rocketAttrs.Epoch_range_to_determine_stepDAC[wRocket-4][0]))).argmin()
            sampleEnd = np.abs(np.array(Epoch_step - pycdf.lib.datetime_to_tt2000(rocketAttrs.Epoch_range_to_determine_stepDAC[wRocket-4][1]))).argmin()

            adjustments = [[2, 2], [10, 11]]

            # determine the values of step for each step
            sampleData = data_dict['step'][0][sampleStart+adjustments[wRocket-4][0]:sampleEnd+adjustments[wRocket-4][1]]
            stepsDigitalVals = [round(sum(sampleData[i*10:(i+1)*10])/10) for i in range(round(len(sampleData)/10))]

            counted = Counter(stepsDigitalVals)
            counted_reduced = dict(sorted(counted.items()))

            # filter out values of the step so that only 101 remain
            loopDict = [counted_reduced for i in range(4)]
            for j in range(len(loopDict)):
                keys = []
                removekeys = []

                for key, val in loopDict[j].items():
                    keys.append(key)

                for i in range(len(keys)-1):
                    if np.abs(keys[i+1] - keys[i]) < 9:
                        if loopDict[j][keys[i+1]] > loopDict[j][keys[i]]:
                            removekeys.append(keys[i])
                        elif loopDict[j][keys[i+1]] < loopDict[j][keys[i]]:
                            removekeys.append(keys[i+1])
                        elif loopDict[j][keys[i]] == loopDict[j][keys[i]]:
                            removekeys.append(keys[i+1])
                    if loopDict[j][keys[i]] < 5:
                        removekeys.append(keys[i])

                for thing in removekeys:
                    counted_reduced.pop(thing,None)

            ###############################################
            # --- REMOVE RC EFFECT BY /10  DOWNSAMPLING ---
            ###############################################
            # description: There is an RC decay that appears on the lower applied voltage values of the LP char. curves
            # downsample the data to only the "end" of the RC decay so as to remove it's appearance in the data
            if downSample_RCeffect:

                # Find where the bottom of the first sweep occurs:
                stepStart = np.abs(data_dict['step'][0] - list(counted_reduced.keys())[0]).argmin()
                n_i_swept = np.array(data_dict['ni_swept'][0][stepStart::])
                n_e_swept = np.array(data_dict['ne_swept'][0][stepStart::])
                sweptStep = np.array(data_dict['step'][0][stepStart::])
                Epoch_sweptCurrent = np.array(data_dict['Epoch_step'][0][stepStart::])
                n_i_swept_Div10, n_e_Swept_Div10, sweptStep_Div10, Epoch_sweptCurrent_Div10 = [], [], [], []

                # Downsample the data
                for i in range(0, len(sweptStep)-10, 10): #each step contains 10 points

                    for j in range(keepThisManyPoints):
                        n_i_swept_Div10.append(n_i_swept[i + (9 - (keepThisManyPoints-1) + j)])
                        n_e_Swept_Div10.append(n_e_swept[i + 9 - (keepThisManyPoints-1) + j])
                        sweptStep_Div10.append(sweptStep[i + 9 - (keepThisManyPoints-1) + j])
                        Epoch_sweptCurrent_Div10.append(Epoch_sweptCurrent[i + 9 - (keepThisManyPoints-1) + j])

                data_dict['ni_swept'][0] = np.array(n_i_swept_Div10)
                data_dict['ne_swept'][0] = np.array(n_e_Swept_Div10)
                data_dict['step'][0] = np.array(sweptStep_Div10)
                data_dict['Epoch_step'][0] = np.array(Epoch_sweptCurrent_Div10)

            # --- linear conversion ---
            # description: linearly convert the "Step" variable to voltage using a linear fit (since we know the min/max) of the applied LP
            voltageRange = rocketAttrs.LPswept_voltage_range[wRocket - 4]
            finalStep_digital = max(counted_reduced)
            initialStep_digital = min(counted_reduced)
            slope = (voltageRange[1] - voltageRange[0]) / (finalStep_digital - initialStep_digital)
            intercept = voltageRange[0]

            stepVoltage = np.array([slope*data_dict['step'][0][i] + intercept if data_dict['step'][0][i] not in [-1,65535] else data_dict['step'][0][i] for i in range(len(data_dict['step'][0]))])
            data_dict = {**data_dict, **{'step_Voltage': [stepVoltage, {'LABLAXIS': 'step Voltage', 'DEPEND_0': 'Epoch_step', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -1e30, 'FORMAT': 'E12.2', 'UNITS': 'Volts', 'VALIDMIN': stepVoltage.min(), 'VALIDMAX': stepVoltage.max(), 'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}

            Done(start_time)

        #####################
        # --- Swept Probe ---
        #####################

        if SECTION_sweptProbe:

            # Calculate the swept current
            n_i_swept = list(data_dict['ni_swept'][0])
            n_e_swept = list(data_dict['ne_swept'][0])
            sweptStep_temp = np.array(data_dict['step_Voltage'][0])
            Epoch_sweptCurrent_temp = np.array(data_dict['Epoch_step'][0])

            # --- prepare the calibration data ---
            if applySweptCalCurve:
                prgMsg('Collecting Swept Calibration Fit Data')
                def sweptCal_Analog_to_Current(analogVal, fitParams):
                    return np.exp(fitParams[0]*analogVal + fitParams[1])

                # get the fit parameters from "csv_to_cdf_LPcal.py"
                data_dict_iowaCal = loadDictFromFile(inputCalFile_Iowa, {})

                # get the calibration data
                fit_params = data_dict_iowaCal['fit_params'][0] # format: [[#1], [#2]]
                Done(start_time)

            prgMsg('Calculating swept LP current')
            if removeADCValsOutofCalRange:
                fitRegions = data_dict_iowaCal['fitRegions'][0]  # format: [[#1], [#2]]

                # remove data out of range: n_i
                for i in range(len(n_i_swept)):
                    if (fitRegions[1][0] > n_i_swept[i]) or (fitRegions[1][1] < n_i_swept[i]):
                        n_i_swept[i] = rocketAttrs.epoch_fillVal

                # remove data out of range: n_e
                for i in range(len(n_e_swept)):
                    if (fitRegions[0][0] > n_e_swept[i]) or (fitRegions[0][1] < n_e_swept[i]):
                        n_e_swept[i] = rocketAttrs.epoch_fillVal

            if breakIntoCurves:
                Epoch_sweptCurrent = []
                sweptCurrent_ne = []
                sweptCurrent_ni = []
                step_sweptVoltage = []

                # Reduce the sweptCurrent data to only include data with stepVoltage >= targetVoltage_min
                for i in range(len(data_dict['step_Voltage'][0])):
                    if data_dict['step_Voltage'][0][i] >= targetVoltage_min:
                        sweptCurrent_ne.append(n_e_swept[i])
                        sweptCurrent_ni.append(n_i_swept[i])
                        Epoch_sweptCurrent.append(Epoch_sweptCurrent_temp[i])
                        step_sweptVoltage.append(sweptStep_temp[i])

                # Reduce data to only specified epoch range that contains good data
                epochIndexLow = np.abs(np.array(Epoch_sweptCurrent) - rocketAttrs.startEndLangmuirBreakIntoCurves[wRocket - 4][0]).argmin()
                epochIndexHigh = np.abs(np.array(Epoch_sweptCurrent) - rocketAttrs.startEndLangmuirBreakIntoCurves[wRocket - 4][1]).argmin()
                sweptCurrent_ne = sweptCurrent_ne[epochIndexLow:epochIndexHigh + 1]
                sweptCurrent_ni = sweptCurrent_ni[epochIndexLow:epochIndexHigh + 1]
                step_sweptVoltage = step_sweptVoltage[epochIndexLow:epochIndexHigh + 1]
                Epoch_sweptCurrent = np.array([pycdf.lib.datetime_to_tt2000(Epoch_sweptCurrent[i]) for i in range(epochIndexLow,epochIndexHigh+1)])

                # Identify Individual sweeps first by the large gap in Epoch value
                indvSweepIndicies = []
                sweepNo = 0

                for i in range(len(Epoch_sweptCurrent)-1):

                    if (Epoch_sweptCurrent[i+1] - Epoch_sweptCurrent[i]) > indvEpochThresh:

                        if sweepNo == 0:
                            indvSweepIndicies.append([0, i])
                        else:
                            indvSweepIndicies.append([indvSweepIndicies[-1][1]+1, i])

                        sweepNo += 1

                # the Above algorithim doesn't include the final sweep, so I add it in here
                indvSweepIndicies.append([indvSweepIndicies[-1][1]+1, len(Epoch_sweptCurrent)])

                # Separate each individual sweep into upleg and downleg
                sweptCurrent_ne_New = []
                sweptCurrent_ni_New = []
                Epoch_sweptCurrent_New = []
                sweptStep_New = []
                for index in range(len(indvSweepIndicies)):
                    indvCurrent_ne = sweptCurrent_ne[indvSweepIndicies[index][0]:indvSweepIndicies[index][1] + 1]
                    indvCurrent_ni = sweptCurrent_ni[indvSweepIndicies[index][0]:indvSweepIndicies[index][1] + 1]
                    indvEpoch = Epoch_sweptCurrent[indvSweepIndicies[index][0]:indvSweepIndicies[index][1] + 1]
                    indvStep = step_sweptVoltage[indvSweepIndicies[index][0]:indvSweepIndicies[index][1] + 1]

                    # # create a temporary n_e - n_i variable that's useful for sorting
                    # ADCcurrent_temp

                    if downSample_RCeffect:
                        # Take the top value of "step" and split the curve based on that
                        middleIndex = np.array(indvStep).argmax()
                    else:
                        # take the top 9 values of 'indvCurrent_ne' and split the curve at the middle value
                        countedIndicies = sorted(range(len(indvCurrent_ne)), key=lambda i: indvCurrent_ne[i])[-9:]
                        countedIndicies.sort()
                        middleIndex = countedIndicies[int((len(countedIndicies) - 1) / 2)]

                    # Break up the curve
                    currentUpLeg_ne = indvCurrent_ne[0:middleIndex]
                    currentUpLeg_ni = indvCurrent_ni[0:middleIndex]
                    epochUpLeg = indvEpoch[0:middleIndex]
                    stepUpLeg = indvStep[0:middleIndex]

                    currentDownLeg_ne = indvCurrent_ne[middleIndex:]
                    currentDownLeg_ni = indvCurrent_ni[middleIndex:]
                    epochDownLeg = indvEpoch[middleIndex:]
                    stepDownLeg = indvStep[middleIndex:]

                    # Store the broken up curve data
                    sweptCurrent_ne_New.append(currentUpLeg_ne)
                    sweptCurrent_ne_New.append(currentDownLeg_ne[::-1])

                    sweptCurrent_ni_New.append(currentUpLeg_ni)
                    sweptCurrent_ni_New.append(currentDownLeg_ni[::-1])

                    Epoch_sweptCurrent_New.append(epochUpLeg)
                    Epoch_sweptCurrent_New.append(epochDownLeg)
                    sweptStep_New.append(stepUpLeg)
                    sweptStep_New.append(stepDownLeg[::-1])

                # prepare final outputs of the breakCurves algorithm
                sweptCurrent_ne = []
                sweptCurrent_ni = []
                Epoch_sweptCurrent = []
                step_sweptVoltage = []

                # Flatten the data. np.flatten does not work for some reason
                for i in range(len(Epoch_sweptCurrent_New)):
                    for thing in sweptCurrent_ne_New[i]:
                        sweptCurrent_ne.append(thing)
                    for thing in sweptCurrent_ni_New[i]:
                        sweptCurrent_ni.append(thing)
                    for thing in Epoch_sweptCurrent_New[i]:
                        Epoch_sweptCurrent.append(pycdf.lib.tt2000_to_datetime(thing))
                    for thing in sweptStep_New[i]:
                        step_sweptVoltage.append(thing)

                ###########################################################
                # --- APPLY THE CALIBRATION CURVES FROM SWEPT_PROBE_CAL ---
                ###########################################################
                if applySweptCalCurve:
                    ne_current = np.array([sweptCal_Analog_to_Current(analogVal, fit_params[0]) if (analogVal != rocketAttrs.epoch_fillVal) else rocketAttrs.epoch_fillVal for analogVal in sweptCurrent_ne])
                    ni_current = np.array([sweptCal_Analog_to_Current(analogVal, fit_params[1]) if (analogVal != rocketAttrs.epoch_fillVal) else rocketAttrs.epoch_fillVal for analogVal in sweptCurrent_ni])
                else:
                    ne_current = np.array(sweptCurrent_ne)
                    ni_current = np.array(sweptCurrent_ni)
                sweptCurrent = np.array(ne_current - ni_current)
                Epoch_sweptCurrent = np.array(Epoch_sweptCurrent)
                step_sweptVoltage = np.array(step_sweptVoltage)
            else:
                if applySweptCalCurve:
                    ne_current = np.array([sweptCal_Analog_to_Current(analogVal, fit_params[0]) if (analogVal != rocketAttrs.epoch_fillVal) else rocketAttrs.epoch_fillVal for analogVal in n_e_swept])
                    ni_current = np.array([sweptCal_Analog_to_Current(analogVal, fit_params[1]) if (analogVal != rocketAttrs.epoch_fillVal) else rocketAttrs.epoch_fillVal for analogVal in n_i_swept])
                    sweptCurrent = np.array([ne_current[i] - ni_current[i] if (ne_current[i] != rocketAttrs.epoch_fillVal) or (ni_current[i] != rocketAttrs.epoch_fillVal) else rocketAttrs.epoch_fillVal for i in range(len(Epoch_sweptCurrent_temp))])

                else:
                    sweptCurrent = np.array([n_e_swept[i] - n_i_swept[i] if (n_e_swept[i] != rocketAttrs.epoch_fillVal) or (n_i_swept[i] != rocketAttrs.epoch_fillVal) else rocketAttrs.epoch_fillVal for i in range(len(Epoch_sweptCurrent_temp))])

                Epoch_sweptCurrent = np.array(Epoch_sweptCurrent_temp)
                step_sweptVoltage = np.array(sweptStep_temp)

            # do a quality assurance check
            for i, val in enumerate(sweptCurrent):
                if np.abs(val) > 1E10:
                    sweptCurrent[i] = rocketAttrs.epoch_fillVal

            units = 'nA' if applySweptCalCurve else 'Analog'
            if applySweptCalCurve:
                if useNanoAmps:
                    units = 'nA'
                    sweptCurrent = np.array(sweptCurrent)
                else:
                    units = 'A'
                    sweptCurrent = np.array(sweptCurrent) / unit_conversion
            else:
                units = 'Analog'


            data_dict = {**data_dict, **{'swept_Current': [sweptCurrent, {'LABLAXIS': 'swept_Current',
                                                            'DEPEND_0': 'Epoch_swept_Current',
                                                            'DEPEND_1': None,
                                                            'DEPEND_2': None,
                                                            'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                            'UNITS': units,
                                                            'VALIDMIN': sweptCurrent.min(),
                                                            'VALIDMAX': sweptCurrent.max(),
                                                            'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}

            data_dict = {**data_dict, **{'Epoch_swept_Current': [Epoch_sweptCurrent, {'LABLAXIS':'Epoch_swept_Current',
                                                                    'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None,
                                                                    'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2', 'UNITS': 'ns',
                                                                    'VALIDMIN': Epoch_sweptCurrent.min(), 'VALIDMAX': Epoch_sweptCurrent.max(), 'VAR_TYPE': 'support_data',
                                                                    'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000',
                                                                    'TIME_SCALE': 'Terrestrial Time','REFERENCE_POSITION': 'Rotating Earth Geoid',
                                                                    'SCALETYP': 'linear'}]}}
            data_dict = {**data_dict, **{'swept_Voltage': [step_sweptVoltage, {'LABLAXIS': 'swept_Voltage',
                                                                          'DEPEND_0': 'Epoch_swept_Current',
                                                                          'DEPEND_1': None,
                                                                          'DEPEND_2': None,
                                                                          'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                                          'UNITS': 'Volts',
                                                                          'VALIDMIN': step_sweptVoltage.min(),
                                                                          'VALIDMAX': step_sweptVoltage.max(),
                                                                          'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
            Done(start_time)

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        if outputData:

            # remove uneeded data.
            removeThese = ['ni', 'ni_swept', 'ne_swept', "Epoch_ni_swept", 'Epoch_ne_swept', 'step', 'EXP_Current', '28V_Monitor', 'Epoch_monitors']
            for thing in removeThese:
                data_dict.pop(thing)

            # rename Boom_Monitors
            data_dict['Swept_Boom_Monitor'] = data_dict.pop('Boom_Monitor_1')
            data_dict['Swept_Boom_Monitor'][1]['DEPEND_0'] = "Epoch_Swept_Boom_Monitor"
            data_dict['Epoch_Swept_Boom_Monitor'] = data_dict.pop('Epoch_monitor_1')

            data_dict['Fixed_Boom_Monitor'] = data_dict.pop('Boom_Monitor_2')
            data_dict['Fixed_Boom_Monitor'][1]['DEPEND_0'] = "Epoch_Fixed_Boom_Monitor"
            data_dict['Epoch_Fixed_Boom_Monitor'] = data_dict.pop('Epoch_monitor_2')

            prgMsg('Creating output file')

            outputPath = f'{outputFolderPath}\\' + f'{fliers[wRocket - 4]}\{fileoutName}'

            globalAttrsMod['Descriptor'] = wInstr[1]
            outputCDFdata(outputPath, data_dict,instrNam='Langmuir')

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

if len(glob(f'{rocketFolderPath}L1\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        L1_to_Langmuir(wRocket, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        L1_to_Langmuir(wRocket, rocketFolderPath, justPrintFileNames,wflyer)