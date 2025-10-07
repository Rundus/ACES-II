# --- L2_estimate_floating_potential.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Use the moment the swept n_e current turns on in the L2 swept LP data
# to get the voltage. This corresponds to ~0 net current and thus the floating potential

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal

from src.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False # Just print the names of files

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

# --- OutputData ---
outputData = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.Science.LangmuirProbes.toggles import FloatingPotentialToggles as fToggles

#######################
# --- MAIN FUNCTION ---
#######################
def L2_estimate_floating_potential(wflyer, justPrintFileNames):

    # --- FILE I/O ---
    rocket_folder_path = DataPaths.ACES_data_folder
    rocketID = ACESII.payload_IDs[wflyer]
    globalAttrsMod = deepcopy(ACESII.global_attributes[wflyer])
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'LangmuirData'

    # Set the paths for the file names
    data_repository = f'{rocket_folder_path}{fToggles.inputPath_modifier}\\{ACESII.fliers[wflyer]}{fToggles.modifier}\\'
    input_files = glob(data_repository+'*_swept.cdf*')
    input_names = [ifile.replace(data_repository, '') for ifile in input_files]
    input_names_searchable = [ifile.replace(fToggles.inputPath_modifier.lower() + '_', '').replace('_v00', '') for ifile in input_names]
    dataFile_name = input_files[0].replace(data_repository,'')

    if justPrintFileNames:
        for i, file in enumerate(input_files):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(os.path.getsize(file) / (10 ** 6), 1)))
        return

    print('\n')
    print(stl.color.UNDERLINE + f'Converting to {fToggles.outputPath_modifier} data for {dataFile_name}' + stl.color.END)
    print('[' + str(0) + ']   ' + str(round(os.path.getsize(input_files[0]) / (10 ** 6), 1)) + 'MiB')

    # --- get the data from the L2 file ---
    stl.prgMsg(f'Loading data from {fToggles.inputPath_modifier} Files')
    data_dict_swept = stl.loadDictFromFile(input_files[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}

    ##################################################################
    # --- Calculate the plasma density from Ion Saturation Current ---
    ##################################################################
    stl.prgMsg('Calculating Fixed ni')

    electron_current = data_dict_swept['electron_swept_current'][0]
    step_voltage = data_dict_swept['step_voltage'][0]
    Epoch = data_dict_swept['Epoch'][0]

    varAttrs = {'LABLAXIS': 'plasma density',
                'DEPEND_0': 'Epoch',
               'DEPEND_1': None,
               'DEPEND_2': None,
               'FILLVAL': ACESII.epoch_fillVal,
               'FORMAT': 'E12.2',
               'UNITS': '!Ncm!A-3!N',
               'VALIDMIN': 0,
               'VALIDMAX': 0,
               'VAR_TYPE': 'data', 'SCALETYP': 'linear'}

    # --- find where the "step-up" point occurs ---
    diff = np.diff(electron_current)
    diff_peaks_idx = np.where(diff > 300)[0]

    floating_voltages = step_voltage[diff_peaks_idx]
    floating_epochs = Epoch[diff_peaks_idx]

    floating_voltages_diff = np.diff(floating_voltages)
    floating_voltages_bad_indicies = np.where(np.abs(floating_voltages_diff)>1)[0]

    floating_voltages = np.delete(floating_voltages,floating_voltages_bad_indicies)
    floating_epochs = np.delete(floating_epochs,floating_voltages_bad_indicies)

    floating_voltages_smooth = scipy.signal.savgol_filter(floating_voltages,
                                                          window_length=20,
                                                          polyorder=3)

    fig, ax = plt.subplots()

    if wflyer == 0:
        ax.axvline(dt.datetime(2022,11,20,17,21,15),color='red',label='E-Region')
    elif wflyer == 1:
        ax.axvline(dt.datetime(2022, 11, 20, 17, 21, 30), color='red',label='E-Region')

    ax.scatter(floating_epochs, -1*floating_voltages, color='tab:blue')
    ax.plot(floating_epochs, -1*floating_voltages_smooth, color='black', linewidth=2)
    # ax.plot(floating_voltages_diff)
    ax.set_title(f'ACESII {rocketID}')
    ax.set_ylabel('Payload Floating Voltage', fontsize=15)
    ax.set_xlabel('Epoch (UTC)',fontsize=15)
    ax.set_ylim(-3,3)
    ax.legend()
    plt.show()



    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:
        stl.prgMsg('Creating output file')
        fileoutName_fixed = f'ACESII_{rocketID}_langmuir_fixed.cdf'
        outputPath = f'{rocket_folder_path}{fToggles.outputPath_modifier}\{ACESII.fliers[wflyer]}\\{fileoutName_fixed}'
        stl.outputCDFdata(outputPath, data_dict_output, instrNam= 'Langmuir Probe', globalAttrsMod=globalAttrsMod)
        stl.Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
target_files = f'{DataPaths.ACES_data_folder}{fToggles.inputPath_modifier}\{ACESII.fliers[wRocket-4]}{fToggles.modifier}\*.cdf'

if len(glob(target_files)) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    L2_estimate_floating_potential(wRocket-4, justPrintFileNames)

