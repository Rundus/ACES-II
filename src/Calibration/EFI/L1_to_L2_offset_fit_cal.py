# --- L1_to_L2_offset_fit_cal.py ---
# Description: For the auroral coordinate vxB calibration values,
# [1] determine the GAIN value (y=mx) needed to bring the Tangent component of the
# E-Field into alignment with the vxB term.
# [2] THEN apply the gain correction term to ALL components
# [3] THEN subtract the vxB term to get a fully calibrated E-Field
# [4] Output the EFI in ENU coordinates as L2 data



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
import numpy as np
from src.my_imports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
from scipy.interpolate import CubicSpline

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 4 -> ACES-II High Flier
# 5 -> ACES-II Low Flier
wRocket = 5
wFiles = [0]
outputData = True

# --- TOGGLES ---
Plot_vxB_rawEField_data = True
Plot_corrected_data = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

def L1_to_L2_offset_fit_cal(wRocket, justPrintFileNames):

    # --- Load the Data ---
    stl.prgMsg(f'Loading data')
    data_dict_vxB = stl.loadDictFromFile(glob('C:\Data\ACESII\calibration\EFI_rkt_convection_calibration\low\*.cdf*')[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}

    # --- --- --- --- --- --- ---
    # --- COLLECTION FIT DATA ---
    # --- --- --- --- --- --- ---
    # Description: Collect a subset of data to calibrate the E-Field data to
    target_time_low = np.abs(data_dict_vxB['Epoch'][0] - dt.datetime(2022, 11, 20, 17, 23, 46)).argmin()
    target_time_high = np.abs(data_dict_vxB['Epoch'][0] - dt.datetime(2022, 11, 20, 17, 29, 00)).argmin()

    # define a working data dictonary
    data_dict_fit_data = {}
    for key, val in data_dict_vxB.items():
        data_dict_fit_data = {**data_dict_fit_data,
                              f'{key}':[deepcopy(val[0][target_time_low:target_time_high]),deepcopy(val[1])]}

    if Plot_vxB_rawEField_data:

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3)
        keys = ['N', 'T', 'p']

        for idx,key in enumerate(keys):
            if key == 'N':
                ax[idx].plot(data_dict_fit_data['Epoch'][0], data_dict_fit_data[f'E_{key}_raw'][0], label=f'E_{key}')
                ax[idx].plot(data_dict_fit_data['Epoch'][0], data_dict_fit_data[f'vxB_{key}'][0], label=f'vxB_{key}',
                             color='red')
            elif key == 'T':
                ax[idx].plot(data_dict_fit_data['Epoch'][0], data_dict_fit_data[f'E_{key}_raw'][0], label=f'E_{key}')
                ax[idx].plot(data_dict_fit_data['Epoch'][0], data_dict_fit_data[f'vxB_{key}'][0], label=f'vxB_{key}',
                             color='red')
            else:
                ax[idx].plot(data_dict_fit_data['Epoch'][0], data_dict_fit_data[f'E_{key}_raw'][0], label=f'E_{key}')
                ax[idx].plot(data_dict_fit_data['Epoch'][0], data_dict_fit_data[f'vxB_{key}'][0], label=f'vxB_{key}',
                             color='red')

            ax[idx].legend()

        plt.show()

    # E_N_correction = 3.31
    E_N_correction = np.sqrt(2)
    E_T_correction = np.sqrt(2)
    E_p_correction = np.sqrt(2)

    if Plot_corrected_data:

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3)
        keys = ['N', 'T', 'p']

        for idx, key in enumerate(keys):
            if key == 'N':
                ax[idx].plot(data_dict_fit_data['Epoch'][0], (E_N_correction) * data_dict_fit_data[f'E_{key}_raw'][0], label=f'E_{key} (Gain Applied)')
            elif key == 'T':
                ax[idx].plot(data_dict_fit_data['Epoch'][0], (E_T_correction) * data_dict_fit_data[f'E_{key}_raw'][0], label=f'E_{key} (Gain Applied)')
            else:
                ax[idx].plot(data_dict_fit_data['Epoch'][0], (E_p_correction) * data_dict_fit_data[f'E_{key}_raw'][0], label=f'E_{key} (Gain Applied)')

            ax[idx].plot(data_dict_fit_data['Epoch'][0], data_dict_fit_data[f'vxB_{key}'][0], label=f'vxB_{key}',color='red')
            ax[idx].legend()

        plt.show()

    # Apply the gain-correction term
    gain_corrections = [E_N_correction, E_T_correction, E_T_correction]
    keys = ['N', 'T', 'p']
    for idx, key in enumerate(keys):
        correct_data = deepcopy(data_dict_vxB[f'E_{key}_raw'][0]) * gain_corrections[idx] - deepcopy( data_dict_vxB[f'vxB_{key}'][0])

        data_dict_output = {**data_dict_output,
                            **{f'E_{key}': [correct_data, deepcopy(data_dict_vxB[f'E_{key}_raw'][1])]}}

    data_dict_output = {**data_dict_output,
                        **{'Epoch': deepcopy(data_dict_vxB['Epoch'])}}

    data_dict_output['E_N'][1]['LABLAXIS'] = 'E_Normal'
    data_dict_output['E_T'][1]['LABLAXIS'] = 'E_Tangent'
    data_dict_output['E_p'][1]['LABLAXIS'] = 'E_Field_Aligned'

    # Add in the electric field magnitude
    E_Field = np.array([deepcopy(data_dict_output['E_N'][0]), deepcopy(data_dict_output['E_T'][0]),
                        deepcopy(data_dict_output['E_p'][0])]).T
    Emag = np.array([np.linalg.norm(val) for val in E_Field])
    data_dict_output = {**data_dict_output,
                        **{'Emag': [Emag, deepcopy(data_dict_output['E_N'][1])]}}
    data_dict_output['Emag'][1]['LABLAXIS'] = 'Emag'

    # --- --- --- --- --- --- --- --- --- --- --- ---
    # --- Estimate neutral wind speed Requirements ---
    # --- --- --- --- --- --- --- --- --- --- --- ---

    # Assume u_p = 0
    u_T = deepcopy(data_dict_output['E_N'][0]) /deepcopy(data_dict_vxB['B_p'][0])
    u_N = -1*deepcopy(data_dict_output['E_T'][0]) / deepcopy(data_dict_vxB['B_p'][0])

    data_dict_output = {**data_dict_output, **{'u_T': [u_T, {'LABLAXIS': 'u_T', 'DEPEND_0': 'Epoch',
                                                      'DEPEND_1': None, 'DEPEND_2': None,
                                                      'FILLVAL': ACESII.epoch_fillVal,
                                                      'FORMAT': 'E12.2', 'UNITS': 'm/s',
                                                      'VAR_TYPE': 'data',
                                                      'SCALETYP': 'linear'}],
                        'u_N': [u_N, {'LABLAXIS': 'u_N', 'DEPEND_0': 'Epoch',
                                      'DEPEND_1': None, 'DEPEND_2': None,
                                      'FILLVAL': ACESII.epoch_fillVal,
                                      'FORMAT': 'E12.2', 'UNITS': 'm/s',
                                      'VAR_TYPE': 'data',
                                      'SCALETYP': 'linear'}]
                        }}

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')
        fileoutName = rf'ACESII_{ACESII.payload_IDs[wRocket-4]}_l2_E_Field_auroral_fullCal.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output, instrNam='EFI')
        stl.Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
L1_to_L2_offset_fit_cal(wRocket, justPrintFileNames)