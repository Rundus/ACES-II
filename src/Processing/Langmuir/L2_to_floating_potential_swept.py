#--- L2_to_floating_potential_swept.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Load the LP swept data. Isolate individual sweeps.
# Interpolate those sweeps and find the voltage where the currents = 0


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

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
wRocket = 4
wSweeps = [22, 23]

# --- DATA OUTPUT ---
plot_peaks = True
plot_individual_sweeps = False
outputData = True

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

    # remove parts of the data I don't need
    if wRocket-4 == 0:
        target_time = [
                        dt.datetime(2022, 11, 20, 17, 21),
                        dt.datetime(2022, 11, 20, 17, 29, 58),
                       ]

    elif wRocket-4 == 1:
        target_time = [
                        dt.datetime(2022, 11, 20, 17, 23, 4),
                        dt.datetime(2022, 11, 20, 17, 28, 18),
                       ]

    for key in ['step_voltage','electron_swept_current','ion_swept_current','Epoch']:
        low_idx = np.abs(data_dict_LP_currents['Epoch'][0]-target_time[0]).argmin()
        high_idx = np.abs(data_dict_LP_currents['Epoch'][0] - target_time[1]).argmin()
        data_dict_LP_currents[key][0] = data_dict_LP_currents[key][0][low_idx:high_idx+1]

    # --- prepare the output data ---
    data_dict_output = {}
    stl.Done(start_time)

    ###########################################
    # --- BREAK DATA INTO INDIVIDUAL CURVES ---
    ###########################################

    # create the current sweeps
    LP_current_swept = deepcopy(data_dict_LP_currents['electron_swept_current'][0] - data_dict_LP_currents['ion_swept_current'][0])
    step_voltage = deepcopy(data_dict_LP_currents['step_voltage'][0])

    # === Use the flipping polarity at the tops/bottoms of the step_voltges to define where the sweeps are ===
    diff = np.diff(step_voltage)
    transistion_idxs = np.where(np.abs(diff) > 0.03)[0]
    minimum_idxs = []
    maximum_idxs = []
    for i in range(1,len(transistion_idxs)-1):
        prev = step_voltage[transistion_idxs[i-1]]
        cur = step_voltage[transistion_idxs[i]]
        nex = step_voltage[transistion_idxs[i+1]]

        # The HighFlight has some outliers in the step_voltages, handle these
        if np.abs(transistion_idxs[i-1] - transistion_idxs[i])<9 or np.abs(transistion_idxs[i+1] - transistion_idxs[i])<9:
            continue
        else:
            # If everything is good, continue to use the polarity of the transition point to determine the sweep
            if (cur - prev) < 0 and (nex - cur) > 0:
                minimum_idxs.append(transistion_idxs[i])
            elif (cur - prev) > 0 and (nex - cur) < 0:
                maximum_idxs.append(transistion_idxs[i])

    # require the distance between minimums be at least 2000 points
    bad_idxs = [i for i in range(1, len(maximum_idxs)-1) if maximum_idxs[i + 1] - maximum_idxs[i-1] == 2000]
    maximum_idxs = np.delete(maximum_idxs, bad_idxs)
    bad_idxs = [i for i in range(1, len(minimum_idxs)-1) if minimum_idxs[i + 1] - minimum_idxs[i-1] == 2000]
    minimum_idxs = np.delete(minimum_idxs, bad_idxs)

    if plot_peaks:
        fig, ax = plt.subplots()
        ax.scatter([i for i in range(len(step_voltage))],step_voltage, s=20, color='tab:blue')

        ax.scatter(minimum_idxs, step_voltage[minimum_idxs], s=20, color='purple')
        ax.scatter(maximum_idxs, step_voltage[maximum_idxs], s=20, color='red')
        plt.show()



    # --- separate sweeps into individual datasets ---
    # For each diff_peak_idx, find the nearest peak "above" this index and nearest minimum "below" this value
    sweeps_current = []
    sweeps_voltage = []
    sweeps_epoch = []
    for i in range(len(maximum_idxs)-1):

        # === rising sweep ===
        sweeps_current.append(LP_current_swept[minimum_idxs[i]:maximum_idxs[i]+1])
        sweeps_voltage.append(step_voltage[minimum_idxs[i]:maximum_idxs[i]+1])
        sweeps_epoch.append(data_dict_LP_currents['Epoch'][0][minimum_idxs[i]])

        # === falling sweep ===
        sweeps_current.append(LP_current_swept[maximum_idxs[i]:minimum_idxs[i+1] + 1])
        sweeps_voltage.append(step_voltage[maximum_idxs[i]:minimum_idxs[i+1] + 1])
        sweeps_epoch.append(data_dict_LP_currents['Epoch'][0][maximum_idxs[i]])

    # --- Interpolate the sweeps to get the zeropoint ---
    floating_potential = []
    for i in range(len(sweeps_current)):
        # Linear interpolate where y=0 to get the x-point
        y_offset = sweeps_current[i]
        y = sweeps_current[i]
        x = sweeps_voltage[i]
        pos = np.where((y_offset[1:] * y_offset[:-1]) <= 0)[0]
        x1 = x[pos]
        x2 = x[pos + 1]
        y1 = y[pos]
        y2 = y[pos + 1]
        xp = (0 - y1) / (y2 - y1) * (x2 - x1) + x1
        floating_potential.append(np.nanmean(xp))

    if plot_individual_sweeps:
        if wSweeps == []:
            for i in range(len(sweeps_current)):
                fig, ax = plt.subplots()
                fig.suptitle(f'Sweep No {i + 1}')
                ax.plot(sweeps_voltage[i], sweeps_current[i])
                ax.axvline(x=floating_potential[i],color='red')
                ax.set_ylabel('Current [nA]')
                ax.set_xlabel('Voltage [V]')
                print(floating_potential[i])
                plt.show()
        else:
            for i in wSweeps:
                fig, ax = plt.subplots()
                fig.suptitle(f'Sweep No {i + 1}')
                ax.plot(sweeps_voltage[i], sweeps_current[i])
                ax.axvline(x=floating_potential[i])
                ax.set_ylabel('Current [nA]')
                ax.set_xlabel('Voltage [V]')
                print(floating_potential[i])
                plt.show()


    # remove any nans from the floating potential data
    from scipy.signal import savgol_filter
    good_idxs = np.where(np.isnan(floating_potential) == False)
    floating_potential = np.array(floating_potential)[good_idxs]
    sweeps_epoch = np.array(sweeps_epoch)[good_idxs]

    bad_idxs = np.where(floating_potential<0)
    floating_potential = np.delete(floating_potential,bad_idxs)
    sweeps_epoch = np.delete(sweeps_epoch,bad_idxs)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:


        data_dict_output = {
                            'Epoch':[np.array(sweeps_epoch),deepcopy(data_dict_LP_currents['Epoch'][1])],
                            'floating_potential':[np.array(floating_potential), {'DEPEND_0':'Epoch','UNITS':'Volts','LABLAXIS':'|Floating Potential|','VAR_TYPE':'data'}],
                            'sweeps_voltage' : [np.array(sweeps_voltage),{'UNITS':'V','LABL_AXIS':'Voltage','VAR_TYPE':'data'}],
                            'sweeps_current': [np.array(sweeps_current),{'UNITS':'nA','LABL_AXIS':'Current','VAR_TYPE':'data'}],
                            }

        stl.prgMsg('Creating output file')
        output_file_name = rf'ACESII_{ACESII.payload_IDs[wRocket-4]}_FloatingPotential.cdf'
        output_path = rf'C:\Data\ACESII\science\payload_potential\{ACESII.fliers[wRocket-4]}\\' + output_file_name
        stl.outputCDFdata(output_path, data_dict_output)
        stl.Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
L2_to_floating_potential_swept(wRocket)