# --- test2_EFI_timeOffset_minimum_variance_analysis.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: minimize the variance between the
# EFI and vxB profiles

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"

import matplotlib.pyplot as plt
import numpy as np

from src.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
E_field_scale = 1/(np.sqrt(2))

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
from src.my_imports import *
import spaceToolsLib as stl
from scipy.interpolate import CubicSpline

#######################
# --- MAIN FUNCTION ---
#######################
def cal2_EFI_timeOffset_minimum_variance_analysis(wRocket):
    # get the EFI data
    data_folder_path = rf'C:\Data\ACESII\calibration\EFI_cal1_spin_vxB_into_rktFrm_calibration\\{ACESII.fliers[wRocket - 4]}\\'
    input_files = glob(data_folder_path + '*vxB_rktFrm.cdf*')
    data_dict_EFI = stl.loadDictFromFile(input_files[0])

    # get the attitude data
    input_files_attitude = glob(rf'C:\Data\ACESII\attitude\{ACESII.fliers[wRocket - 4]}\*.cdf*')
    data_dict_attitude = stl.loadDictFromFile(input_files_attitude[0])

    ##############################################
    # --- Reduce data to only relevant regions ---
    ##############################################

    # Entire Data set
    # target_times = [dt.datetime(2022, 11, 20, 17, 24, 00),
    #                 dt.datetime(2022, 11, 20, 17, 27, 40)]

    # Beginning Region - Quiet
    target_times = [dt.datetime(2022, 11, 20, 17, 24, 00),
                    dt.datetime(2022, 11, 20, 17, 25, 00)]

    # Middle region - Very Active
    # target_times = [dt.datetime(2022, 11, 20, 17, 25, 00),
    #                 dt.datetime(2022, 11, 20, 17, 26, 40)]

    # # End region - Active
    # target_times = [dt.datetime(2022, 11, 20, 17, 26, 40),
    #                 dt.datetime(2022, 11, 20, 17, 27, 40)]


    low_idx = np.abs(data_dict_EFI['Epoch'][0] - target_times[0]).argmin()
    high_idx = np.abs(data_dict_EFI['Epoch'][0] - target_times[1]).argmin()

    ################################
    # --- Loop Over time-offsets ---
    ################################

    # Define some setup variables
    Epoch_EFI = deepcopy(data_dict_EFI['Epoch'][0][low_idx:high_idx])
    T0_EFI = stl.EpochTo_T0_Rocket(Epoch_EFI, T0=dt.datetime(2022, 11, 20, 17, 20))
    vxB_ENU = np.array([data_dict_EFI['vxB_E_RingCore'][0][low_idx:high_idx],
                        data_dict_EFI['vxB_N_RingCore'][0][low_idx:high_idx],
                        data_dict_EFI['vxB_Up_RingCore'][0][low_idx:high_idx]
                        ]).T
    E_rkt = np.array([data_dict_EFI['E_X'][0][low_idx:high_idx],
                      data_dict_EFI['E_Y'][0][low_idx:high_idx],
                      data_dict_EFI['E_Z'][0][low_idx:high_idx]]).T

    T0_attitude = stl.EpochTo_T0_Rocket(deepcopy(data_dict_attitude['Epoch'][0]),
                                        T0=dt.datetime(2022, 11, 20, 17, 20))

    # loop through time-offsets and calculate the variances
    ChiVal = []
    time_offset = []

    N = 20
    deltaT = np.linspace(-0.2, -0.06, N)
    for tme in tqdm(deltaT):
        DCM = np.zeros(shape=(len(Epoch_EFI), 3, 3))
        T0_interp = deepcopy(T0_attitude) + tme
        for i in range(3):
            for j in range(3):
                cs = CubicSpline(T0_interp, data_dict_attitude[f'a{i + 1}{j + 1}'][0])
                DCM[:, i, j] = cs(T0_EFI)

        # --- Apply the DCM to the EFI ---

        E_ENU = (1 / np.sqrt(2)) * np.array([np.matmul(DCM[i], E_rkt[i]) for i in range(len(DCM))])

        # --- Calculate the variance ---
        ChiVal.append(np.sum(np.power(vxB_ENU-E_ENU,2)/np.abs(vxB_ENU)))
        time_offset.append(tme)

    fig, ax = plt.subplots()
    ax.scatter(time_offset,ChiVal, label='$T_{min}=$' + f'{round(time_offset[np.array(ChiVal).argmin()],6)}')
    ax.legend()
    ax.set_ylabel('$|\chi|$ Values')
    ax.set_xlabel('Time Offsets [s]')
    ax.grid(True)
    plt.show()




# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
cal2_EFI_timeOffset_minimum_variance_analysis(wRocket)

