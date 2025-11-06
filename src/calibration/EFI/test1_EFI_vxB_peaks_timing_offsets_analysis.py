# --- test1_EFI_vxB_peaks_timing_offsets_analysis.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: determine the timing between the EFI peaks and vxB. Use this to
# determine what the time-offset between the vxB and DCM should be.

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"

import matplotlib.pyplot as plt
import numpy as np
import spaceToolsLib

from src.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
time_scale = 1000 # converts to ms
E_field_scale = 1

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

# --- OutputData ---
outputData = True
filter_data = True
plot_filtered_data = False
plot_fits = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.my_imports import *
import spaceToolsLib as stl
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline

#######################
# --- MAIN FUNCTION ---
#######################
def cal2_EFI_vxB_offsets_analysis(wRocket):

    # --- FILE I/O ---
    data_folder_path = rf'{DataPaths.ACES_data_folder}\\calibration\\\EFI_cal1_spin_vxB_into_rktFrm_calibration\\{ACESII.fliers[wRocket-4]}\\'

    # get the EFI files
    input_files = glob(data_folder_path + '*vxB_rktFrm.cdf*')

    # --- get the data from the L2 file ---
    stl.prgMsg(f'Loading data Files')
    data_dict_EFI = stl.loadDictFromFile(input_files[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}

    ################################################################
    # --- [1] Find the deltaT between the peaks of vxB and EFI ---
    ################################################################

    # low-pass butterworth filter the EFI data before finding peaks


    if filter_data:
        stl.prgMsg('Filtering Data')
        filterType = 'bandpass'
        fs = 1E9/(pycdf.lib.datetime_to_tt2000(data_dict_EFI['Epoch'][0][200000+1]) - pycdf.lib.datetime_to_tt2000(data_dict_EFI['Epoch'][0][200000]) )
        order = 2
        fs_cutoff_low = 0.4
        fs_cutoff_high = 1.5

        E_filtered = [[], [], []]
        for idx, key in enumerate(['X', 'Y', 'Z']):
            E_filtered[idx] = stl.butterFilter().butter_filter(
                data=data_dict_EFI[f'E_{key}'][0],
                order=order,
                lowcutoff=fs_cutoff_low,
                highcutoff=fs_cutoff_high,
                filtertype=filterType,
                fs=fs
                             )

        if plot_filtered_data:

            fig, ax = plt.subplots(3)

            for idx, key in enumerate(['X','Y','Z']):
                ax[idx].plot(data_dict_EFI['Epoch'][0],data_dict_EFI[f'E_{key}'][0], color='tab:blue')
                ax[idx].plot(data_dict_EFI['Epoch'][0], E_filtered[idx], color='tab:red')
                ax[idx].set_ylim(-0.15, 0.15)
                ax[idx].set_xlim(dt.datetime(2022, 11, 20, 17, 24), dt.datetime(2022, 11, 20, 17, 25))

            plt.show()

        data_dict_EFI['E_X'][0] = deepcopy(np.array(E_filtered[0]))
        data_dict_EFI['E_Y'][0] = deepcopy(np.array(E_filtered[1]))
        data_dict_EFI['E_Z'][0] = deepcopy(np.array(E_filtered[2]))
        stl.Done(start_time)

    # Find the peaks in the EFI data
    stl.prgMsg('Finding Peaks')
    peaks_EFI = [[[], []], [[], []], [[], []]] # stores the high and low peaks of the EFI
    peaks_vxB = [[[], []], [[], []], [[], []]] # stores the high and low peaks of the EFI

    peak_params = {'height': 0.04,
                   'distance': 1500,
                   'width': 500}

    for i in range(2):
        # Find the peaks in the EFI data
        peaks_EFI[0][i], _ = find_peaks(x=np.power(-1, (i + 2))*data_dict_EFI['E_X'][0], **peak_params)
        peaks_EFI[1][i], _ = find_peaks(x=np.power(-1, (i + 2))*data_dict_EFI['E_Y'][0], **peak_params)
        peaks_EFI[2][i], _ = find_peaks(x=np.power(-1, (i + 2))*data_dict_EFI['E_Z'][0], **peak_params)

        # Find the peaks in the vxB data
        peaks_vxB[0][i], _ = find_peaks(x=np.power(-1, (i + 2))*(data_dict_EFI['vxB_X'][0]), **peak_params)
        peaks_vxB[1][i], _ = find_peaks(x=np.power(-1, (i + 2))*(data_dict_EFI['vxB_Y'][0]), **peak_params)
        peaks_vxB[2][i], _ = find_peaks(x=np.power(-1, (i + 2))*(data_dict_EFI['vxB_Z'][0]), **peak_params)

    stl.Done(start_time)

    # [1] get the times of the EFI and vxB peaks in seconds since T0.
    stl.prgMsg('Finding Time Delay')
    T0 = dt.datetime(2022, 11, 20, 17, 20)
    peak_times_EFI = [[[], []], [[], []], [[], []]]  # stores the high and low peaks of the EFI
    peak_times_vxB = [[[], []], [[], []], [[], []]]  # stores the high and low peaks of the EFI
    for i in range(2):
        for idx, key in enumerate(['X', 'Y', 'Z']):
            if key !='X':
                # EFI
                times = data_dict_EFI['Epoch'][0][peaks_EFI[idx][i]]
                peak_times_EFI[idx][i] = stl.EpochTo_T0_Rocket(times, T0=T0)

                # vxB
                times = data_dict_EFI['Epoch'][0][peaks_vxB[idx][i]]
                peak_times_vxB[idx][i] = stl.EpochTo_T0_Rocket(times, T0=T0)

    # [2] Limit peak data starting from T_start
    start_datetime = dt.datetime(2022, 11, 20, 17, 24, 00, 500000)
    T_start_idx = np.abs(data_dict_EFI['Epoch'][0] - start_datetime).argmin()

    for i in range(2): # Shorten the peaks/data to only indices above T_start_idx
        for idx, key in enumerate(['X', 'Y', 'Z']):

            if key !='X':

                # EFI
                finder_idxs =np.where(peaks_EFI[idx][i] > T_start_idx)
                peaks_EFI[idx][i] = peaks_EFI[idx][i][finder_idxs]
                peak_times_EFI[idx][i] = peak_times_EFI[idx][i][finder_idxs]

                # vxB
                finder_idxs = np.where(peaks_vxB[idx][i] > T_start_idx)
                peaks_vxB[idx][i] = peaks_vxB[idx][i][finder_idxs]
                peak_times_vxB[idx][i] = peak_times_vxB[idx][i][finder_idxs]

    # [3] Find the time delay between subsequent peaks in the X-Y data
    # Note: we need the vxB peak times since our linear trend needs to be correct the vxB time, i.e. (time correction) = (slope)*vxB_time + intercept

    deltaT = [[], [], []]
    deltaT_epoch = [[], [], []]
    for idx, key in enumerate(['X', 'Y', 'Z']):
        for i in range(len(peaks_EFI[idx][0])-1):
            if idx != 0: # Ignore the X-axis since it has no clear peaks

                if key == 'Y': # the Y-axis starts on minimums
                    # calculate the deltaT for minimums
                    deltaT[idx].append(peak_times_vxB[idx][1][i] - peak_times_EFI[idx][1][i])
                    deltaT_epoch[idx].append(data_dict_EFI['Epoch'][0][peaks_EFI[idx][1][i]])

                    # calculate the deltaT for maximums
                    deltaT[idx].append(peak_times_vxB[idx][0][i] - peak_times_EFI[idx][0][i])
                    deltaT_epoch[idx].append(data_dict_EFI['Epoch'][0][peaks_EFI[idx][0][i]])
                elif key == 'Z': # the Z-axis startson maximums
                    # calculate the deltaT for maximums
                    deltaT[idx].append(peak_times_vxB[idx][0][i] - peak_times_EFI[idx][0][i])
                    deltaT_epoch[idx].append(data_dict_EFI['Epoch'][0][peaks_EFI[idx][0][i]])

                    # calculate the deltaT for minimums
                    deltaT[idx].append(peak_times_vxB[idx][1][i] - peak_times_EFI[idx][1][i])
                    deltaT_epoch[idx].append(data_dict_EFI['Epoch'][0][peaks_EFI[idx][1][i]])
            else:
                continue

    stl.Done(start_time)

    # [4] Fit the DeltaT
    stl.prgMsg('Fitting DeltaT')
    def fitFunc(x,A,B):
        return A*x+B
    params_1, cov = curve_fit(fitFunc, stl.EpochTo_T0_Rocket(deltaT_epoch[1],T0=T0), np.array(deltaT[1]))
    params_2, cov = curve_fit(fitFunc, stl.EpochTo_T0_Rocket(deltaT_epoch[2],T0=T0), np.array(deltaT[2]))
    rogers_values = [0.3842, 0.37736]
    stl.Done(start_time)

    # ---------------
    # PLOT EVERYTHING
    # ---------------
    if plot_fits:
        stl.prgMsg('Plotting everything')
        fig, ax = plt.subplots(3)
        xData = data_dict_EFI['Epoch'][0]

        for idx, key in enumerate(['Y', 'Z']):
            dat_idx = idx +1
            # EFI
            yData = data_dict_EFI[f'E_{key}'][0]
            ax[idx].plot(xData, yData*E_field_scale, color='tab:blue', label=f'E_{key}')
            ax[idx].scatter(xData[peaks_EFI[dat_idx][0]], yData[peaks_EFI[dat_idx][0]]*E_field_scale, s=40, color='green')
            ax[idx].scatter(xData[peaks_EFI[dat_idx][1]], yData[peaks_EFI[dat_idx][1]]*E_field_scale, s=40, color='purple')

            # vxB
            yData = data_dict_EFI[f'vxB_{key}'][0]
            ax[idx].plot(xData, yData, color='tab:red', label=f'vxB_{key}')
            ax[idx].scatter(xData[peaks_vxB[dat_idx][0]], yData[peaks_vxB[dat_idx][0]], s=40, color='green')
            ax[idx].scatter(xData[peaks_vxB[dat_idx][1]], yData[peaks_vxB[dat_idx][1]], s=40, color='purple')

            # figure configuration
            ax[idx].set_ylabel(f'{key}-Axis [V/m]')
            ax[idx].set_ylim(-0.15, 0.15)
            ax[idx].set_xlim(start_datetime - dt.timedelta(seconds=4), start_datetime+dt.timedelta(seconds=30))
            ax[idx].legend()

        ax[2].plot(deltaT_epoch[1], np.array(deltaT[1])*time_scale, label='Y', color='tab:blue')
        ax[2].plot(deltaT_epoch[1], fitFunc(stl.EpochTo_T0_Rocket(deltaT_epoch[1], T0=T0), *params_1)*time_scale, label=f'Y [slope:{round(params_1[0],6)}, b:{round(params_1[1],6)}]', color='tab:blue', alpha=0.5, linestyle='--')

        ax[2].plot(deltaT_epoch[2], np.array(deltaT[2])*time_scale, label='Z', color='tab:red')
        ax[2].plot(deltaT_epoch[2], fitFunc(stl.EpochTo_T0_Rocket(deltaT_epoch[2], T0=T0), *params_2)*time_scale, label=f'Z [slope:{round(params_2[0],6)}, b:{round(params_2[1],6)}]', color='tab:red', alpha=0.5, linestyle='--')

        ax[2].grid(which='both')
        ax[2].set_ylabel('Peaks $\Delta$T [ms]\nT_vxB - T_E')
        ax[2].legend()
        plt.tight_layout()
        plt.show()
        stl.Done(start_time)

    stl.prgMsg('Shifting vxB by fit vals')

    # Apply the fit functions to the vxB data - shift the vxB data backward in time to match the EFI data
    Epoch_T0 = np.array(stl.EpochTo_T0_Rocket(data_dict_EFI['Epoch'][0], T0=T0))

    # --- Interpolate the (Epoch_new, vxB_component) data on the OLD timebase ---

    # Y-Axis
    Epoch_new = Epoch_T0 + fitFunc(Epoch_T0, *params_1)
    E_Y_new = np.interp(Epoch_T0,Epoch_new,data_dict_EFI['E_Y'][0])

    # Z-Axis
    Epoch_new = Epoch_T0 + fitFunc(Epoch_T0, *params_2)
    E_Z_new = np.interp(Epoch_T0, Epoch_new, data_dict_EFI['E_Z'][0])

    # X-Axis
    vxB_X_new = np.interp(Epoch_T0, Epoch_new, data_dict_EFI['vxB_X'][0])

    # Interpolate the NEW attitude DCM matrix elements into OLD timebase
    DCM_new = np.zeros(shape=(len(data_dict_EFI['DCM'][0]), 3, 3))
    for idx1 in range(3):
        for idx2 in range(3):
            DCM_new[:, idx1, idx2] = np.interp(Epoch_T0, Epoch_new, data_dict_EFI['DCM'][0][:,idx1,idx2])
    stl.Done(start_time)


    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')
        data_dict_output = {
                            'Epoch': deepcopy(data_dict_EFI['Epoch']),
                            'Epoch_EFI_T0': [Epoch_T0, deepcopy(data_dict_EFI['Epoch'][1])],
                            'Y_fitParams': [np.array([params_1]), {}],
                            'Z_fitParams': [np.array([params_2]), {}],
                            'deltaT_Y': [deltaT_epoch[1], {}],
                            'deltaT_Z': [deltaT_epoch[2], {}],
                            'Y_fit': [np.array(fitFunc(stl.EpochTo_T0_Rocket(deltaT_epoch[1], T0=T0), *params_1)), {'DEPEND_0':'deltaT_Y'}],
                            'Z_fit': [np.array(fitFunc(stl.EpochTo_T0_Rocket(deltaT_epoch[2], T0=T0), *params_2)), {'DEPEND_0':'deltaT_Z'}],
                            'DCM': [np.array(data_dict_EFI['DCM'][0]), deepcopy(data_dict_EFI['DCM'][1])],
                            'E_Y':[np.array(E_Z_new), deepcopy(data_dict_EFI['E_Y'][1])],
                            'E_Z':[np.array(E_Z_new), deepcopy(data_dict_EFI['E_Z'][1])],
                            'E_X':deepcopy(data_dict_EFI['E_X'])
                            }

        fileoutName = f'ACESII_{ACESII.payload_IDs[wRocket-4]}_l1_EFI_vxB_rktFrm.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\calibration\EFI_cal2_timing_offset_calibration\\{ACESII.fliers[wRocket-4]}\\'+ fileoutName
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
cal2_EFI_vxB_offsets_analysis(wRocket)

