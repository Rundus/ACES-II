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
import time
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
wRocket = 5

# --- DATA OUTPUT ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.ACESII.Processing.Langmuir.toggles import L1toL2_SweptLPToggles as fToggles
from warnings import filterwarnings # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.
filterwarnings("ignore")
from collections import Counter

def L1_to_L2_langmuir_swept(wflyer, justPrintFileNames):

    # --- Get ACESII rocket Attributes ---
    rocket_folder_path = DataPaths.ACES_data_folder
    rocketID = ACESII.payload_IDs[wflyer]
    global_attrs_mod = deepcopy(ACESII.global_attributes[wflyer])
    global_attrs_mod['Logical_source'] = global_attrs_mod['Logical_source'] + 'Langmuir'
    global_attrs_mod['Descriptor'] = 'Langmuir'

    # Set the paths for the file names
    data_repository = f'{rocket_folder_path}{fToggles.inputPath_modifier}\{ACESII.fliers[wflyer]}{fToggles.modifier}\\'

    input_files_ne = glob(data_repository + '*_ne_swept_*')
    input_names_ne = [ifile.replace(data_repository, '') for ifile in input_files_ne]

    input_files_ni = glob(data_repository + '*_ni_swept_*')
    input_names_ni = [ifile.replace(data_repository, '') for ifile in input_files_ni]

    input_files_step = glob(data_repository + '*_lp_step_*')
    input_names_step = [ifile.replace(data_repository, '') for ifile in input_files_step]

    input_names_searchable = [ifile.replace(fToggles.inputPath_modifier.lower() + '_', '').replace('_v00', '') for ifile in input_names_ne+input_names_ni+input_names_step]

    if justPrintFileNames:
        for i, file in enumerate(input_files_ne+input_files_ni+input_files_step):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return

    print('\n')
    print(stl.color.UNDERLINE + f'Creating swept langmuir data for {ACESII.fliers[wflyer]} flyer' + stl.color.END)

    # --- get the data from the tmCDF file ---
    stl.prgMsg('Loading data from L1Files')
    data_dict_step = stl.loadDictFromFile(input_files_step[0])
    data_dict_ne = stl.loadDictFromFile(input_files_ne[0])
    data_dict_ni = stl.loadDictFromFile(input_files_ni[0])

    # --- prepare the output data ---
    data_dict_output = {}
    stl.Done(start_time)

    ###################################
    # --- Calibrate STEP to Voltage ---
    ###################################
    # description: The swept LP "step" variable contains analog values that of 100 steps (101 values in total) for values up and down
    # determining a nominal voltage to assign each of these values will indicate what voltage was applied
    # to the swept langmuir probes.

    stl.prgMsg('Calculating swept step voltage')

    # calculate the epoch index of beginning/end of sample range
    step_start = np.abs(np.array(data_dict_step['Epoch_step'][0] - ACESII.LP_epoch_range_to_determine_step_DAC[wflyer][0])).argmin()
    step_end = np.abs(np.array(data_dict_step['Epoch_step'][0] - ACESII.LP_epoch_range_to_determine_step_DAC[wflyer][1])).argmin()

    # determine the values of step for each step
    adjustments = [[2, 2], [10, 11]]
    step_data = deepcopy(data_dict_step['step'][0][step_start+adjustments[wflyer][0]:step_end+adjustments[wflyer][1]])
    stepsDigitalVals = [round(sum(step_data[i*10:(i+1)*10])/10) for i in range(round(len(step_data)/10))]

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


    #####################################
    # --- CONVERT STEP DAC TO VOLTAGE ---
    #####################################
    # --- linear conversion ---
    # description: linearly convert the "Step" variable to voltage using a linear fit (since we know the min/max) of the applied LP
    probe_sweep_voltage_range = ACESII.LP_swept_voltage_range[wflyer]
    slope = (probe_sweep_voltage_range[1] - probe_sweep_voltage_range[0]) / (max(counted_reduced) - min(counted_reduced))
    intercept = probe_sweep_voltage_range[0]
    step_voltage = slope * deepcopy(data_dict_step['step'][0]) + intercept
    data_dict_output = {**data_dict_output, **{'step_voltage': [step_voltage, {'LABLAXIS': 'Probe Voltage', 'DEPEND_0': 'Epoch',
                                                                'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': ACESII.epoch_fillVal,
                                                                'FORMAT': 'E12.2', 'UNITS': 'Volts',
                                                                'VALIDMIN': step_voltage.min(),
                                                                'VALIDMAX': step_voltage.max(),
                                                                'VAR_TYPE': 'data',
                                                                'SCALETYP': 'linear'}]}}

    #######################################
    # --- TRIM DATA TO FIRST FULL SWEEP ---
    #######################################
    # Find where the bottom of the first sweep occurs:
    steps_start_idx = np.abs(data_dict_step['step'][0] - list(counted_reduced.keys())[0]).argmin()
    n_i_swept = np.array(data_dict_ni['ni_swept'][0][steps_start_idx::])
    n_e_swept = np.array(data_dict_ne['ne_swept'][0][steps_start_idx::])
    Epoch = np.array(data_dict_step['Epoch_step'][0][steps_start_idx::])

    data_dict_output = {**data_dict_output, **{'ion_swept_current': [n_i_swept, {'LABLAXIS': 'ion current',
                                                                   'DEPEND_0': 'Epoch',
                                                                    'FILLVAL': ACESII.epoch_fillVal,
                                                                    'FORMAT': 'E12.2', 'UNITS': 'nA',
                                                                    'VALIDMIN': n_i_swept.min(),
                                                                    'VALIDMAX': n_i_swept.max(),
                                                                    'VAR_TYPE': 'data',
                                                                    'SCALETYP': 'linear'}]}}

    data_dict_output = {**data_dict_output, **{'electron_swept_current': [n_e_swept, {'LABLAXIS': 'electron current',
                                                                   'DEPEND_0': 'Epoch',
                                                                   'FILLVAL': ACESII.epoch_fillVal,
                                                                   'FORMAT': 'E12.2', 'UNITS': 'nA',
                                                                   'VALIDMIN': n_e_swept.min(),
                                                                   'VALIDMAX': n_e_swept.max(),
                                                                   'VAR_TYPE': 'data',
                                                                   'SCALETYP': 'linear'}]}}

    data_dict_output = {**data_dict_output, **{'Epoch': [Epoch, {'LABLAXIS': 'Epoch', 'DEPEND_0': None,
                                                                    'FILLVAL': ACESII.epoch_fillVal,
                                                                   'FORMAT': 'E12.2', 'UNITS': 'UTC',
                                                                   'VALIDMIN': Epoch.min(),
                                                                   'VALIDMAX': Epoch.max(),
                                                                   'VAR_TYPE': 'data',
                                                                   'SCALETYP': 'linear'}]}}

    data_dict_output['step_voltage'][0] = data_dict_output['step_voltage'][0][steps_start_idx::]
    stl.Done(start_time)


    ###############################################
    # --- REMOVE RC EFFECT BY /10  DOWNSAMPLING ---
    ###############################################
    if fToggles.down_sample_RC_effect_bool:

        stl.prgMsg('Removing RC effect')

        # description: There is an RC decay that appears on the lower applied voltage values of the LP char. curves
        # downsample the data to only the "end" of the RC decay so as to remove it's appearance in the data
        n_i_swept_Div10, n_e_Swept_Div10, swept_voltage_Div10, Epoch_Div10 = [], [], [], []

        # Downsample the data
        for i in range(0, len(step_voltage)-10, 10): #each step contains 10 points

            for j in range(fToggles.keep_this_many_points):
                n_i_swept_Div10.append(data_dict_output['ion_swept_current'][0][i + (9 - (fToggles.keep_this_many_points-1) + j)])
                n_e_Swept_Div10.append(data_dict_output['electron_swept_current'][0][i + 9 - (fToggles.keep_this_many_points-1) + j])
                swept_voltage_Div10.append(data_dict_output['step_voltage'][0][i + 9 - (fToggles.keep_this_many_points-1) + j])
                Epoch_Div10.append(data_dict_output['Epoch'][0][i + 9 - (fToggles.keep_this_many_points-1) + j])

        data_dict_output['ion_swept_current'][0] = np.array(n_i_swept_Div10)
        data_dict_output['electron_swept_current'][0] = np.array(n_e_Swept_Div10)
        data_dict_output['step_voltage'][0] = np.array(swept_voltage_Div10)
        data_dict_output['Epoch'][0] = np.array(Epoch_Div10)
        stl.Done(start_time)

    # #####################
    # # --- Swept Probe ---
    # #####################
    #
    # if SECTION_sweptProbe:
    #
    #     # Calculate the swept current
    #     n_i_swept = list(data_dict['ni_swept'][0])
    #     n_e_swept = list(data_dict['ne_swept'][0])
    #     sweptStep_temp = np.array(data_dict['step_Voltage'][0])
    #     Epoch_sweptCurrent_temp = np.array(data_dict['Epoch_step'][0])
    #
    #     # --- prepare the calibration data ---
    #     if applySweptCalCurve:
    #         stl.prgMsg('Collecting Swept calibration Fit Data')
    #         def sweptCal_Analog_to_Current(analogVal, fitParams):
    #             return np.exp(fitParams[0]*analogVal + fitParams[1])
    #
    #         # get the fit parameters from "csv_to_cdf_LPcal.py"
    #         data_dict_iowaCal = loadDictFromFile(inputCalFile_Iowa, {})
    #
    #         # get the calibration data
    #         fit_params = data_dict_iowaCal['fit_params'][0] # format: [[#1], [#2]]
    #         stl.Done(start_time)
    #
    #     stl.prgMsg('Calculating swept LP current')
    #     if removeADCValsOutofCalRange:
    #         fitRegions = data_dict_iowaCal['fitRegions'][0]  # format: [[#1], [#2]]
    #
    #         # remove data out of range: n_i
    #         for i in range(len(n_i_swept)):
    #             if (fitRegions[1][0] > n_i_swept[i]) or (fitRegions[1][1] < n_i_swept[i]):
    #                 n_i_swept[i] = rocketAttrs.epoch_fillVal
    #
    #         # remove data out of range: n_e
    #         for i in range(len(n_e_swept)):
    #             if (fitRegions[0][0] > n_e_swept[i]) or (fitRegions[0][1] < n_e_swept[i]):
    #                 n_e_swept[i] = rocketAttrs.epoch_fillVal
    #
    #     if breakIntoCurves:
    #         Epoch_sweptCurrent = []
    #         sweptCurrent_ne = []
    #         sweptCurrent_ni = []
    #         step_sweptVoltage = []
    #
    #         # Reduce the sweptCurrent data to only include data with stepVoltage >= targetVoltage_min
    #         for i in range(len(data_dict['step_Voltage'][0])):
    #             if data_dict['step_Voltage'][0][i] >= targetVoltage_min:
    #                 sweptCurrent_ne.append(n_e_swept[i])
    #                 sweptCurrent_ni.append(n_i_swept[i])
    #                 Epoch_sweptCurrent.append(Epoch_sweptCurrent_temp[i])
    #                 step_sweptVoltage.append(sweptStep_temp[i])
    #
    #         # Reduce data to only specified epoch range that contains good data
    #         epochIndexLow = np.abs(np.array(Epoch_sweptCurrent) - rocketAttrs.startEndLangmuirBreakIntoCurves[wRocket - 4][0]).argmin()
    #         epochIndexHigh = np.abs(np.array(Epoch_sweptCurrent) - rocketAttrs.startEndLangmuirBreakIntoCurves[wRocket - 4][1]).argmin()
    #         sweptCurrent_ne = sweptCurrent_ne[epochIndexLow:epochIndexHigh + 1]
    #         sweptCurrent_ni = sweptCurrent_ni[epochIndexLow:epochIndexHigh + 1]
    #         step_sweptVoltage = step_sweptVoltage[epochIndexLow:epochIndexHigh + 1]
    #         Epoch_sweptCurrent = np.array([pycdf.lib.datetime_to_tt2000(Epoch_sweptCurrent[i]) for i in range(epochIndexLow,epochIndexHigh+1)])
    #
    #         # Identify Individual sweeps first by the large gap in Epoch value
    #         indvSweepIndicies = []
    #         sweepNo = 0
    #
    #         for i in range(len(Epoch_sweptCurrent)-1):
    #
    #             if (Epoch_sweptCurrent[i+1] - Epoch_sweptCurrent[i]) > indvEpochThresh:
    #
    #                 if sweepNo == 0:
    #                     indvSweepIndicies.append([0, i])
    #                 else:
    #                     indvSweepIndicies.append([indvSweepIndicies[-1][1]+1, i])
    #
    #                 sweepNo += 1
    #
    #         # the Above algorithim doesn't include the final sweep, so I add it in here
    #         indvSweepIndicies.append([indvSweepIndicies[-1][1]+1, len(Epoch_sweptCurrent)])
    #
    #         # Separate each individual sweep into upleg and downleg
    #         sweptCurrent_ne_New = []
    #         sweptCurrent_ni_New = []
    #         Epoch_sweptCurrent_New = []
    #         sweptStep_New = []
    #         for index in range(len(indvSweepIndicies)):
    #             indvCurrent_ne = sweptCurrent_ne[indvSweepIndicies[index][0]:indvSweepIndicies[index][1] + 1]
    #             indvCurrent_ni = sweptCurrent_ni[indvSweepIndicies[index][0]:indvSweepIndicies[index][1] + 1]
    #             indvEpoch = Epoch_sweptCurrent[indvSweepIndicies[index][0]:indvSweepIndicies[index][1] + 1]
    #             indvStep = step_sweptVoltage[indvSweepIndicies[index][0]:indvSweepIndicies[index][1] + 1]
    #
    #             # # create a temporary n_e - n_i variable that's useful for sorting
    #             # ADCcurrent_temp
    #
    #             if downSample_RCeffect:
    #                 # Take the top value of "step" and split the curve based on that
    #                 middleIndex = np.array(indvStep).argmax()
    #             else:
    #                 # take the top 9 values of 'indvCurrent_ne' and split the curve at the middle value
    #                 countedIndicies = sorted(range(len(indvCurrent_ne)), key=lambda i: indvCurrent_ne[i])[-9:]
    #                 countedIndicies.sort()
    #                 middleIndex = countedIndicies[int((len(countedIndicies) - 1) / 2)]
    #
    #             # Break up the curve
    #             currentUpLeg_ne = indvCurrent_ne[0:middleIndex]
    #             currentUpLeg_ni = indvCurrent_ni[0:middleIndex]
    #             epochUpLeg = indvEpoch[0:middleIndex]
    #             stepUpLeg = indvStep[0:middleIndex]
    #
    #             currentDownLeg_ne = indvCurrent_ne[middleIndex:]
    #             currentDownLeg_ni = indvCurrent_ni[middleIndex:]
    #             epochDownLeg = indvEpoch[middleIndex:]
    #             stepDownLeg = indvStep[middleIndex:]
    #
    #             # Store the broken up curve data
    #             sweptCurrent_ne_New.append(currentUpLeg_ne)
    #             sweptCurrent_ne_New.append(currentDownLeg_ne[::-1])
    #
    #             sweptCurrent_ni_New.append(currentUpLeg_ni)
    #             sweptCurrent_ni_New.append(currentDownLeg_ni[::-1])
    #
    #             Epoch_sweptCurrent_New.append(epochUpLeg)
    #             Epoch_sweptCurrent_New.append(epochDownLeg)
    #             sweptStep_New.append(stepUpLeg)
    #             sweptStep_New.append(stepDownLeg[::-1])
    #
    #         # prepare final outputs of the breakCurves algorithm
    #         sweptCurrent_ne = []
    #         sweptCurrent_ni = []
    #         Epoch_sweptCurrent = []
    #         step_sweptVoltage = []
    #
    #         # Flatten the data. np.flatten does not work for some reason
    #         for i in range(len(Epoch_sweptCurrent_New)):
    #             for thing in sweptCurrent_ne_New[i]:
    #                 sweptCurrent_ne.append(thing)
    #             for thing in sweptCurrent_ni_New[i]:
    #                 sweptCurrent_ni.append(thing)
    #             for thing in Epoch_sweptCurrent_New[i]:
    #                 Epoch_sweptCurrent.append(pycdf.lib.tt2000_to_datetime(thing))
    #             for thing in sweptStep_New[i]:
    #                 step_sweptVoltage.append(thing)
    #
    #         ###########################################################
    #         # --- APPLY THE CALIBRATION CURVES FROM SWEPT_PROBE_CAL ---
    #         ###########################################################
    #         if applySweptCalCurve:
    #             ne_current = np.array([sweptCal_Analog_to_Current(analogVal, fit_params[0]) if (analogVal != rocketAttrs.epoch_fillVal) else rocketAttrs.epoch_fillVal for analogVal in sweptCurrent_ne])
    #             ni_current = np.array([sweptCal_Analog_to_Current(analogVal, fit_params[1]) if (analogVal != rocketAttrs.epoch_fillVal) else rocketAttrs.epoch_fillVal for analogVal in sweptCurrent_ni])
    #         else:
    #             ne_current = np.array(sweptCurrent_ne)
    #             ni_current = np.array(sweptCurrent_ni)
    #         sweptCurrent = np.array(ne_current - ni_current)
    #         Epoch_sweptCurrent = np.array(Epoch_sweptCurrent)
    #         step_sweptVoltage = np.array(step_sweptVoltage)
    #     else:
    #         if applySweptCalCurve:
    #             ne_current = np.array([sweptCal_Analog_to_Current(analogVal, fit_params[0]) if (analogVal != rocketAttrs.epoch_fillVal) else rocketAttrs.epoch_fillVal for analogVal in n_e_swept])
    #             ni_current = np.array([sweptCal_Analog_to_Current(analogVal, fit_params[1]) if (analogVal != rocketAttrs.epoch_fillVal) else rocketAttrs.epoch_fillVal for analogVal in n_i_swept])
    #             sweptCurrent = np.array([ne_current[i] - ni_current[i] if (ne_current[i] != rocketAttrs.epoch_fillVal) or (ni_current[i] != rocketAttrs.epoch_fillVal) else rocketAttrs.epoch_fillVal for i in range(len(Epoch_sweptCurrent_temp))])
    #
    #         else:
    #             sweptCurrent = np.array([n_e_swept[i] - n_i_swept[i] if (n_e_swept[i] != rocketAttrs.epoch_fillVal) or (n_i_swept[i] != rocketAttrs.epoch_fillVal) else rocketAttrs.epoch_fillVal for i in range(len(Epoch_sweptCurrent_temp))])
    #
    #         Epoch_sweptCurrent = np.array(Epoch_sweptCurrent_temp)
    #         step_sweptVoltage = np.array(sweptStep_temp)
    #
    #     # do a quality assurance check
    #     for i, val in enumerate(sweptCurrent):
    #         if np.abs(val) > 1E10:
    #             sweptCurrent[i] = rocketAttrs.epoch_fillVal
    #
    #     units = 'nA' if applySweptCalCurve else 'Analog'
    #     if applySweptCalCurve:
    #         if useNanoAmps:
    #             units = 'nA'
    #             sweptCurrent = np.array(sweptCurrent)
    #         else:
    #             units = 'A'
    #             sweptCurrent = np.array(sweptCurrent) / unit_conversion
    #     else:
    #         units = 'Analog'
    #
    #
    #     data_dict = {**data_dict, **{'swept_Current': [sweptCurrent, {'LABLAXIS': 'swept_Current',
    #                                                     'DEPEND_0': 'Epoch_swept_Current',
    #                                                     'DEPEND_1': None,
    #                                                     'DEPEND_2': None,
    #                                                     'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
    #                                                     'UNITS': units,
    #                                                     'VALIDMIN': sweptCurrent.min(),
    #                                                     'VALIDMAX': sweptCurrent.max(),
    #                                                     'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
    #
    #     data_dict = {**data_dict, **{'Epoch_swept_Current': [Epoch_sweptCurrent, {'LABLAXIS':'Epoch_swept_Current',
    #                                                             'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None,
    #                                                             'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2', 'UNITS': 'ns',
    #                                                             'VALIDMIN': Epoch_sweptCurrent.min(), 'VALIDMAX': Epoch_sweptCurrent.max(), 'VAR_TYPE': 'support_data',
    #                                                             'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000',
    #                                                             'TIME_SCALE': 'Terrestrial Time','REFERENCE_POSITION': 'Rotating Earth Geoid',
    #                                                             'SCALETYP': 'linear'}]}}
    #     data_dict = {**data_dict, **{'swept_Voltage': [step_sweptVoltage, {'LABLAXIS': 'swept_Voltage',
    #                                                                   'DEPEND_0': 'Epoch_swept_Current',
    #                                                                   'DEPEND_1': None,
    #                                                                   'DEPEND_2': None,
    #                                                                   'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
    #                                                                   'UNITS': 'Volts',
    #                                                                   'VALIDMIN': step_sweptVoltage.min(),
    #                                                                   'VALIDMAX': step_sweptVoltage.max(),
    #                                                                   'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
    #     stl.Done(start_time)
    #


    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        stl.prgMsg('Creating output file')

        # remove uneeded data.
        keep_keys = ['EXP_Current',
                     '28V_Monitor',
                     'Boom_Monitor_1',
                     'Epoch_monitor_1']

        # update the output data_dict
        for key in keep_keys:
            data_dict_output = {**data_dict_output,
                                **{key: deepcopy(data_dict_ne[key])}
                                }

        # rename Boom_Monitors - BOOM MONITOR #1 is swept, #2 is Fixed
        data_dict_output['swept_boom_monitor'] = data_dict_output.pop('Boom_Monitor_1')
        data_dict_output['swept_boom_monitor'][1]['DEPEND_0'] = "Epoch_boom_monitor"
        data_dict_output['Epoch_boom_monitor'] = data_dict_output.pop('Epoch_monitor_1')
        data_dict_output['28V_Monitor'][1]['DEPEND_0'] = "Epoch_boom_monitor"


        output_file_name = rf'ACESII_{rocketID}_l2_langmuir_swept.cdf'
        output_path = f'{rocket_folder_path}\\{fToggles.outputPath_modifier}\\{ACESII.fliers[wflyer]}\{output_file_name}'
        stl.outputCDFdata(output_path, data_dict_output, instrNam='Langmuir', globalAttrsMod=global_attrs_mod)
        stl.Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
target_files = f'{DataPaths.ACES_data_folder}{fToggles.inputPath_modifier}\{ACESII.fliers[wRocket-4]}{fToggles.modifier}\*_ni_swept_*'

if len(glob(target_files)) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    L1_to_L2_langmuir_swept(wRocket-4, justPrintFileNames)