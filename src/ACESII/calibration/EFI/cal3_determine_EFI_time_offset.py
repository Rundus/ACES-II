# --- cal3_determine_EFI_time_offset.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: minimize the variance between the
# IDEAL EFI from cal2 and the time-adjusted
# EFI data using small deltaT's on the attitue DCM

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"

import matplotlib.pyplot as plt

from src.ACESII.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
E_field_scale = 1/(np.sqrt(2))

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

# --- OutputData ---
outputData = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.ACESII.my_imports import *
import spaceToolsLib as stl
from scipy.interpolate import CubicSpline

#######################
# --- MAIN FUNCTION ---
#######################
def cal3_determine_EFI_time_offset(wRocket):

    # get the EFI files with vxB included
    data_folder_path = rf'C:\Data\ACESII\calibration\EFI\cal1_spin_vxB_into_rktFrm_calibration\\{ACESII.fliers[wRocket - 4]}\\'
    input_files = glob(data_folder_path + '*vxB_rktFrm*')
    data_dict_EFI = stl.loadDictFromFile(input_files[0])

    # get the IDEAL EFI files
    input_files = glob(rf'C:\Data\ACESII\calibration\EFI\cal2_form_ideal_EFI_dataset\\{ACESII.fliers[wRocket - 4]}\\' + '*.cdf*')
    data_dict_EFI_ideal = stl.loadDictFromFile(input_files[0])

    # get the attitude data
    data_dict_attitude = stl.loadDictFromFile(glob(rf"C:\Data\ACESII\\attitude\\{ACESII.fliers[wRocket - 4]}\\*.cdf*")[0])

    ##############################################
    # --- Reduce data to only relevant regions ---
    ##############################################

    # Entire Data set
    target_times = [dt.datetime(2022, 11, 20, 17, 24, 00),
                    dt.datetime(2022, 11, 20, 17, 27, 40)]

    # Beginning Region - Quiet
    # target_times = [dt.datetime(2022, 11, 20, 17, 24, 00),
    #                 dt.datetime(2022, 11, 20, 17, 25, 00)]

    # Middle region - Very Active
    # target_times = [dt.datetime(2022, 11, 20, 17, 25, 00),
    #                 dt.datetime(2022, 11, 20, 17, 26, 40)]

    # # End region - Active
    # target_times = [dt.datetime(2022, 11, 20, 17, 26, 40),
    #                 dt.datetime(2022, 11, 20, 17, 27, 40)]


    low_idx = np.abs(data_dict_EFI['Epoch'][0] - target_times[0]).argmin()
    high_idx = np.abs(data_dict_EFI['Epoch'][0] - target_times[1]).argmin()

    #########################
    # --- SETUP VARIABLES ---
    #########################

    # Define some setup variables
    Epoch_EFI = deepcopy(data_dict_EFI['Epoch'][0][low_idx:high_idx])
    T0_EFI = stl.EpochTo_T0_Rocket(Epoch_EFI, T0=dt.datetime(2022, 11, 20, 17, 20))
    vxB_ENU = np.array([data_dict_EFI['vxB_E'][0][low_idx:high_idx],
                        data_dict_EFI['vxB_N'][0][low_idx:high_idx],
                        data_dict_EFI['vxB_Up'][0][low_idx:high_idx]
                        ]).T
    E_rkt = np.array([data_dict_EFI['E_X'][0][low_idx:high_idx],
                      data_dict_EFI['E_Y'][0][low_idx:high_idx],
                      data_dict_EFI['E_Z'][0][low_idx:high_idx]]).T

    T0_attitude = stl.EpochTo_T0_Rocket(deepcopy(data_dict_attitude['Epoch'][0]), T0=dt.datetime(2022, 11, 20, 17, 20))

    E_ENU_ideal = np.array([
        data_dict_EFI_ideal['E_E'][0][low_idx:high_idx],
        data_dict_EFI_ideal['E_N'][0][low_idx:high_idx],
        data_dict_EFI_ideal['E_Up'][0][low_idx:high_idx],
    ]).T

    ################################
    # --- Loop Over time-offsets ---
    ################################
    # Loop Toggles
    ChiVal = []
    offsets = []
    slopes = []
    N = 2
    M = 2
    deltaT_offset = np.linspace(-0.15, -0.09, N)
    deltaT_slope = np.linspace(0, 0.00006, M)

    # Perform the Loop
    for idx1 in tqdm(range(N)):
        for idx2 in range(M):
            DCM = np.zeros(shape=(len(Epoch_EFI), 3, 3))
            T0_interp = deepcopy(T0_attitude) + (deepcopy(T0_attitude)*deltaT_slope[idx2] + deltaT_offset[idx1])
            for i in range(3):
                for j in range(3):
                    cs = CubicSpline(T0_interp, data_dict_attitude[f'a{i + 1}{j + 1}'][0])
                    DCM[:, i, j] = cs(T0_EFI)

            # --- Apply the DCM to the EFI and subtract vxB---
            E_ENU = (1 / np.sqrt(2)) * np.array([np.matmul(DCM[i], E_rkt[i]) for i in range(len(DCM))]) - vxB_ENU

            # --- Calculate the variance ---
            ChiVal.append(np.abs(np.sum(np.power(E_ENU_ideal-E_ENU,2)/E_ENU_ideal)))
            offsets.append(deltaT_offset[idx1])
            slopes.append(deltaT_slope[idx2])

    fig, ax = plt.subplots()
    cmap_ax = ax.pcolormesh(np.reshape(offsets,(N,M)),
                            np.reshape(slopes,(N,M)),
                            np.reshape(ChiVal,(N,M)),
                            cmap='viridis',
                            norm='log')

    ax.set_ylabel('Slopes')
    ax.set_xlabel('Time Offsets [s]')
    ax.grid(True)
    fig.colorbar(cmap_ax)
    min_idx =np.array(ChiVal).argmin()
    ax.legend(['$\DeltaT=$' + f'{round(offsets[min_idx], 6)}\n' + f'm={round(slopes[min_idx], 10)}'])
    # print('$\DeltaT=$' + f'{round(offsets[min_idx],6)}\n')
    # print(f'm={round(slopes[min_idx],10)}')
    plt.show()


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
cal3_determine_EFI_time_offset(wRocket)

