#--- L2_to_floating_potential_swept.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Load the LP swept data. Isolate individual sweeps.
# Interpolate those sweeps and find the voltage where the currents = 0


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"



# --- --- --- --- ---
from src.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

#################
# --- TOGGLES ---
#################

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

# --- DATA OUTPUT ---
outputData = False

plot_peaks = False
plot_individual_sweeps = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def L2_to_floating_potential_swept(wRocket):

    # --- Get ACES-II rocket Attributes ---
    stl.prgMsg('Loading data')
    input_files = glob(rf'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket-4]}\\*langmuir_swept.cdf*')
    data_dict_LP_currents = stl.loadDictFromFile(input_files[0])

    # --- prepare the output data ---
    data_dict_output = {}
    stl.Done(start_time)

    ###########################################
    # --- BREAK DATA INTO INDIVIDUAL CURVES ---
    ###########################################

    # create the current sweeps
    LP_current_swept = deepcopy(data_dict_LP_currents['electron_swept_current'][0] - data_dict_LP_currents['ion_swept_current'][0])
    step_voltage = deepcopy(data_dict_LP_currents['step_voltage'][0])

    # Find the peaks in the step voltage - these will isolate our sweeps
    peaks, _ = find_peaks(step_voltage, height=2.05)
    minimums,_ = find_peaks(-1*step_voltage, height=4.66)

    if plot_peaks:
        fig, ax = plt.subplots()
        ax.plot(step_voltage)
        ax.scatter(peaks, step_voltage[peaks],s=50,color='red')
        ax.scatter(minimums, step_voltage[minimums], s=50, color='purple')
        plt.show()

    # separate sweeps into individual datasets
    sweeps_current = []
    sweeps_voltage = []
    for i in range(len(peaks)):
        low_idx = minimums[i]
        high_idx = peaks[i]
        sweeps_current.append(LP_current_swept[low_idx:high_idx])
        sweeps_voltage.append(step_voltage[low_idx:high_idx])

    if plot_individual_sweeps:
        if wSweeps == []:
        elif:
        for i in range(len(sweeps_current)):
            fig, ax = plt.subplots()
            fig.suptitle(f'Sweep No {i+1}')
            ax.plot(sweeps_voltage[i],sweeps_current[i])
            ax.set_ylabel('Current [nA]')
            ax.set_xlabel('Voltage [V]')
            plt.show()


    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        stl.prgMsg('Creating output file')
        output_file_name = rf'ACESII_{ACESII.payload_IDs[wRocket-4]}_FloatingPotential.cdf'
        output_path = rf'C:\Data\ACESII\science\payload_potential\{ACESII.fliers[wRocket-4]}\\' + output_file_name
        stl.outputCDFdata(output_path, data_dict_output)
        stl.Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
L2_to_floating_potential_swept(wRocket)