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

#################
# --- TOGGLES ---
#################
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import scipy
from src.ACESII.Processing.Langmuir.toggles import L1toL2_FixedLPToggles as fToggles

#######################
# --- MAIN FUNCTION ---
#######################

def L1_to_L2_langmuir_fixed(wflyer, justPrintFileNames):

    # --- Get ACESII rocket Attributes ---
    rocket_folder_path = DataPaths.ACES_data_folder
    rocketID = ACESII.payload_IDs[wflyer]
    global_attrs_mod = deepcopy(ACESII.global_attributes[wflyer])
    global_attrs_mod['Logical_source'] = global_attrs_mod['Logical_source'] + 'Langmuir'
    global_attrs_mod['Descriptor'] = 'Langmuir'

    # Set the paths for the file names
    data_repository = f'{rocket_folder_path}{fToggles.inputPath_modifier}\{ACESII.fliers[wflyer]}{fToggles.modifier}\\'
    input_files = glob(data_repository+'*_ni_playback*')
    input_names = [ifile.replace(data_repository, '') for ifile in input_files]
    input_names_searchable = [ifile.replace(fToggles.inputPath_modifier.lower() + '_', '').replace('_v00', '') for ifile in input_names]

    if justPrintFileNames:
        for i, file in enumerate(input_files):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return

    print(stl.color.UNDERLINE + f'Converting L1 to L2 langmuir data for {ACESII.fliers[wflyer]} flyer' + stl.color.END)

    # --- get the data from the CDF file ---
    data_dict = stl.loadDictFromFile(inputFilePath=input_files[0])

    # --- prepare the output ---
    data_dict_output = {}

    #################################
    # --- Fixed Probe calibration ---
    #################################
    # description: Loads in the fixed calibrated data (provided by scott, found in missionAttributes_OUTDATED.py).

    def calFunction_fixed(x, A, B):
        y = A*x + B
        return y

    cal_currents = [] # cal currents in NANOAMPs
    analog_vals = []

    # convert calibration data to current (in nanoamps)
    # NOTE: WE don't care about the "open case"
    for key, val in ACESII.LP_fixed_cal_resistances[wflyer].items():
        if key != 'Open':
            analog_vals.append(val)
            cal_currents.append(1E9*ACESII.LP_fixed_probe_bias[wflyer] / key)

    analog_vals, cal_currents = np.array(analog_vals), np.array(cal_currents)

    # apply a natural log scale to the data in order to prepare it for fitting
    cal_currents = np.log(-1*cal_currents)

    # Fit a linear line to the log'd data
    params, cov = scipy.optimize.curve_fit(calFunction_fixed, analog_vals, cal_currents, maxfev=10000)

    # output the fit parameters
    data_dict_output = {**data_dict_output, **{'cal_curve_fit_params': [np.array(params), {}]}}

    # calculate the probe current from the calibration curve
    LP_fixed_current = np.exp(calFunction_fixed(data_dict['ni'][0], params[0], params[1]))
    LP_units = 'nA'

    data_dict_output = {**data_dict_output, **{'fixed_current': [LP_fixed_current, {'LABLAXIS': 'current',
                                                                 'DEPEND_0': 'Epoch',
                                                                 'DEPEND_1': None,
                                                                 'DEPEND_2': None,
                                                                 'FILLVAL': ACESII.epoch_fillVal,
                                                                 'FORMAT': 'E12.2',
                                                                 'UNITS': LP_units,
                                                                 'VALIDMIN': LP_fixed_current.min(),
                                                                 'VALIDMAX': LP_fixed_current.max(),
                                                                 'VAR_TYPE': 'data',
                                                                 'SCALETYP': 'linear'}]}}


    # --- quality assurance ---
    # if wflyer == 0:
    #     for i in range(len(caldCurrent)):
    #         if np.abs(caldCurrent[i]) > 2000:
    #             caldCurrent[i] = rocketAttrs.epoch_fillVal
    # elif wflyer == 1:
    #     for i in range(len(caldCurrent)):
    #         if np.abs(caldCurrent[i]) > 2000:
    #             caldCurrent[i] = rocketAttrs.epoch_fillVal


    # # --- --- --- --- --- --- ----
    # # --- FIXED ERROR ANALYSIS ---
    # # --- --- --- --- --- --- ----
    #
    # # (1) The error in the plasma density comes from the error in the fit coefficents for the conversion between I_probe and ADC
    # # % error n_i = sqrt(a^2 * delta a^2 + delta b^2)
    #
    # # (2) to get delta a, delta b we can do a reduced chi-square analysis on the fit:
    # # chi^2 = 1/nu SUM (f(x_i) - y_i)^2 / (sigmaX_i^2 + sigmaY_i^2)
    # # here:
    # # (i) f: a*ADC+ b
    # # (ii) y_i: ith value of ln(I_probe)
    # # (iii): sigmaX_i: error in ADC value
    # # (iv): SigmaY_i: error in ln(I_probe)
    # #
    # # we can work out what delta a, delta b are from optimization calculus to get:
    # # (i): delta a = Sxx/DELTA
    # # (ii): delta b = S/DELTA
    # # for S = SUM(1/sqrt(sigmaX_i^2 +sigmaY_i^2)), Sxx = SUM(x_i^2/sqrt(sigmaX_i^2 +sigmaY_i^2))
    #
    # # (3) we need delta ADC and delta Ln(I_probe).
    # # (a) The best we can do is delta ADC = +/- 0.5 ADC
    # # (b) for Ln(I_probe), we assume voltage error of +/- 0.01V and known resistance error or 1% we have"
    # IprobePercentError = np.sqrt((0.01/5.05)**2 + 0.01**2) # = about 0.01
    # # (c) we convert ^^^ to delta Ln(I_probe) by delta ln(I_probe)= \delta I_probe /Iprobe BUT delta I_probe ~ 0.01 * Iprobe / Iprobe ===> 0.01
    # deltaLnProbe = IprobePercentError
    #
    #
    # Sxx = sum([ analog_vals[i]**2 / np.sqrt((0.5**2) + deltaLnProbe**2) for i in range(len(analog_vals))])
    # Sx = sum([analog_vals[i] / np.sqrt((0.5 ** 2) + deltaLnProbe**2) for i in range(len(analog_vals))])
    # S = sum([1 / np.sqrt((0.5 ** 2) + deltaLnProbe**2) for i in range(len(analog_vals))])
    # DELTA = S*Sxx - (Sx)**2
    # delta_a = Sxx/DELTA
    # delta_b = S/DELTA
    # delta_n = np.array([np.sqrt(  (parameters[0]*delta_a)**2 + (delta_b)**2  )])
    #
    # data_dict = {**data_dict, **{'fixed_current': [caldCurrent, {'LABLAXIS': 'current',
    #                                                 'DEPEND_0': 'fixed_Epoch',
    #                                                 'DEPEND_1': None,
    #                                                 'DEPEND_2': None,
    #                                                 'FILLVAL': rocketAttrs.epoch_fillVal,
    #                                                 'FORMAT': 'E12.2',
    #                                                 'UNITS': fixedLPunits,
    #                                                 'VALIDMIN': caldCurrent.min(),
    #                                                 'VALIDMAX': caldCurrent.max(),
    #                                                 'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
    # data_dict = {**data_dict, **{'ni_percent_error': [delta_n, {'LABLAXIS': 'ni_percent_error',
    #                                                              'DEPEND_0': None,
    #                                                              'DEPEND_1': None,
    #                                                              'DEPEND_2': None,
    #                                                              'FILLVAL': rocketAttrs.epoch_fillVal,
    #                                                              'FORMAT': 'E12.2',
    #                                                              'UNITS': 'cm^-3',
    #                                                              'VALIDMIN': 0,
    #                                                              'VALIDMAX': 1,
    #                                                              'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]}}
    # data_dict['fixed_Epoch'] = data_dict.pop('Epoch_ni') # rename to fixed
    # stl.Done(start_time)
    #
    #
    #
    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        # remove uneeded data.
        keep_keys = ['EXP_Current',
                     '28V_Monitor',
                     'Boom_Monitor_2',
                     'Epoch_monitor_2',
                     'Epoch_ni']

        # update the output data_dict
        for key in keep_keys:
            data_dict_output = {**data_dict_output,
                                **{key: deepcopy(data_dict[key])}
                                }

        # rename Boom_Monitors - BOOM MONITOR #1 is swept, #2 is Fixed
        data_dict_output['fixed_boom_monitor'] = data_dict_output.pop('Boom_Monitor_2')
        data_dict_output['fixed_boom_monitor'][1]['DEPEND_0'] = "Epoch_boom_monitor"
        data_dict_output['Epoch_boom_monitor'] = data_dict_output.pop('Epoch_monitor_2')
        data_dict_output['Epoch'] = data_dict_output.pop('Epoch_ni')
        data_dict_output['28V_Monitor'][1]['DEPEND_0'] = "Epoch_boom_monitor"

        stl.prgMsg('Creating output file')
        output_file_name = rf'ACESII_{rocketID}_l2_langmuir_fixed.cdf'
        output_path = f'{rocket_folder_path}\\{fToggles.outputPath_modifier}\\{ACESII.fliers[wflyer]}\{output_file_name}'
        stl.outputCDFdata(output_path, data_dict_output, instrNam='Langmuir', globalAttrsMod=global_attrs_mod)
        stl.Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
target_files = f'{DataPaths.ACES_data_folder}{fToggles.inputPath_modifier}\{ACESII.fliers[wRocket-4]}{fToggles.modifier}\*.cdf'

if len(glob(target_files)) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    L1_to_L2_langmuir_fixed(wRocket-4, justPrintFileNames)