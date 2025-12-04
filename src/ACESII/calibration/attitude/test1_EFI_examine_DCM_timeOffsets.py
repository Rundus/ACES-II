# --- test2_EFI_examine_DCM_timeOFFsets.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: take a small portion of rocket EFI data and attitude data. Despin
# them through direct applicaiton/interpolation then compare to a time-shifted DCM
# to see what the DCM should be shifted by

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
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

target_times = [dt.datetime(2022,11,20,17,23,55),
                dt.datetime(2022,11,20,17,24,20)]
# target_times = [dt.datetime(2022, 11, 20, 17, 26, 30),
#                     dt.datetime(2022, 11, 20, 17, 27, 40)]

init_offset = -0.0693
init_slope = -7E-6

# BEST PARAMETEr FIT SO FAR
# time_offset = -0.058 # in seconds. Determined from cal3
# time_slope = -0.000074 # in seconds. Determined from cal3

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.ACESII.my_imports import *
from scipy.interpolate import CubicSpline
import spaceToolsLib as stl
from matplotlib.widgets import Slider

#######################
# --- MAIN FUNCTION ---
#######################
def test2_EFI_examine_DCM_timeOFFsets(wRocket):

    # get the EFI data
    data_folder_path = rf'C:\Data\ACESII\calibration\EFI\\cal1_spin_vxB_into_rktFrm_calibration\\{ACESII.fliers[wRocket-4]}\\'
    input_files = glob(data_folder_path + '*vxB_rktFrm.cdf*')
    data_dict_EFI = stl.loadDictFromFile(input_files[0])

    # get the Ideal EFI data
    data_folder_path = rf'C:\Data\ACESII\calibration\EFI\\cal2_form_ideal_EFI_dataset\\{ACESII.fliers[wRocket - 4]}\\'
    input_files = glob(data_folder_path + '*.cdf*')
    data_dict_EFI_ideal = stl.loadDictFromFile(input_files[0])

    # get the attitude data
    input_files_attitude = glob(rf'C:\Data\ACESII\attitude\{ACESII.fliers[wRocket-4]}\*.cdf*')
    data_dict_attitude = stl.loadDictFromFile(input_files_attitude[0])

    ############################
    # --- Reduce EFI Dataset ---
    ############################
    low_idx = np.abs(data_dict_EFI['Epoch'][0]-target_times[0]).argmin()
    high_idx = np.abs(data_dict_EFI['Epoch'][0] - target_times[1]).argmin()
    E_rkt = np.array([data_dict_EFI['E_X'][0][low_idx:high_idx],
                  data_dict_EFI['E_Y'][0][low_idx:high_idx],
                  data_dict_EFI['E_Z'][0][low_idx:high_idx]]).T

    Epoch_EFI = deepcopy(data_dict_EFI['Epoch'][0][low_idx:high_idx])
    T0_EFI = stl.EpochTo_T0_Rocket(Epoch_EFI,T0=dt.datetime(2022,11,20,17,20))

    # --- Interpolate Attitude DCM onto EFI timebase ---
    DCM = np.zeros(shape=(len(Epoch_EFI),3,3))
    T0_attitude = stl.EpochTo_T0_Rocket(deepcopy(data_dict_attitude['Epoch'][0]),T0=dt.datetime(2022,11,20,17,20))
    for i in range(3):
        for j in range(3):
            cs = CubicSpline(T0_attitude,data_dict_attitude[f'a{i+1}{j+1}'][0])
            DCM[:, i, j] = cs(T0_EFI)

    # --- Apply the DCM to the EFI ---
    E_ENU = (1/np.sqrt(2))*np.array([np.matmul(DCM[i],E_rkt[i]) for i in range(len(E_rkt))])
    E_ENU_ideal = np.array([data_dict_EFI_ideal['E_E'][0][low_idx:high_idx],
                            data_dict_EFI_ideal['E_N'][0][low_idx:high_idx],
                            data_dict_EFI_ideal['E_Up'][0][low_idx:high_idx]
                            ]).T
    vxB_ENU = np.array([data_dict_EFI['vxB_E'][0][low_idx:high_idx],
                        data_dict_EFI['vxB_N'][0][low_idx:high_idx],
                        data_dict_EFI['vxB_Up'][0][low_idx:high_idx]
                        ]).T

    #########################
    # --- Plot Everything ---
    #########################

    fig, ax = plt.subplots(3)
    comp_names = ['E','N','Up']

    # East Component
    axNo = 0
    ax[axNo].plot(Epoch_EFI, E_ENU[:, axNo]- vxB_ENU[:, axNo],label=f'E_{comp_names[axNo]}'+'/$\sqrt{2}$', color='tab:blue')
    ax[axNo].plot(Epoch_EFI, E_ENU_ideal[:, axNo], label=f'E_{comp_names[axNo]}' + ' (Ideal)', color='tab:green', alpha=0.5)
    line1_adjust, = ax[axNo].plot(Epoch_EFI, E_ENU[:, axNo] - vxB_ENU[:, axNo], label=f'E_{comp_names[axNo]}' + '/$\sqrt{2} - vxB$ (time-adjusted)', color='tab:purple')
    ax[axNo].set_ylim(-0.03, 0.03)
    ax[axNo].axhline(0, color='red')
    ax[axNo].legend()

    # North Component
    axNo = 1
    ax[axNo].plot(Epoch_EFI,E_ENU[:, axNo]- vxB_ENU[:, axNo], label=f'E_{comp_names[axNo]}' + '/$\sqrt{2}$', color='tab:blue')
    ax[axNo].plot(Epoch_EFI, E_ENU_ideal[:, axNo], label=f'E_{comp_names[axNo]}' + ' (Ideal)', color='tab:green', alpha=0.5)
    line2_adjust, = ax[axNo].plot(Epoch_EFI, E_ENU[:, axNo] - vxB_ENU[:, axNo], label=f'E_{comp_names[axNo]}' + '/$\sqrt{2} - vxB$ (time-adjusted)', color='tab:purple')
    ax[axNo].set_ylim(-0.03, 0.03)
    ax[axNo].axhline(0,color='red')
    ax[axNo].legend()

    # Up Component
    axNo = 2
    ax[axNo].plot(Epoch_EFI, E_ENU[:, axNo]- vxB_ENU[:, axNo], label=f'E_{comp_names[axNo]}' + '/$\sqrt{2}$', color='tab:blue')
    ax[axNo].plot(Epoch_EFI, E_ENU_ideal[:, axNo], label=f'E_{comp_names[axNo]}' + ' (Ideal)', color='tab:green', alpha=0.5)
    line3_adjust, = ax[axNo].plot(Epoch_EFI, E_ENU[:, axNo] - vxB_ENU[:, axNo], label=f'E_{comp_names[axNo]}' + '/$\sqrt{2} - vxB$ (time-adjusted)', color='tab:purple')
    ax[axNo].set_ylim(-0.03, 0.03)
    ax[axNo].axhline(0, color='red')
    ax[axNo].legend()

    # Add the time slider
    axOffset = fig.add_axes([0.25,0.07,0.65,0.02])

    offset_slider = Slider(
        ax=axOffset,
        label='Offset [s]',
        valmin=-0.20,
        valmax=-0.01,
        valinit=init_offset,
    )

    axSlope = fig.add_axes([0.25, 0.025, 0.65, 0.02])

    slope_slider = Slider(
        ax=axSlope,
        label='Slope ',
        valmin=-0.0002,
        valmax=0.0002,
        valinit=init_slope,
    )

    def update(val):

        # calculate the new DCM
        DCM_adjust = np.zeros(shape=(len(Epoch_EFI), 3, 3))
        T0_attitude = stl.EpochTo_T0_Rocket(deepcopy(data_dict_attitude['Epoch'][0]), T0=dt.datetime(2022, 11, 20, 17, 20))

        T0_attitude_new = T0_attitude*(1 + slope_slider.val)+offset_slider.val
        for i in range(3):
            for j in range(3):
                cs = CubicSpline(T0_attitude_new, data_dict_attitude[f'a{i + 1}{j + 1}'][0])
                DCM_adjust[:, i, j] = cs(T0_EFI)

        # Apply the time-adjusted DCM to the EFI data
        E_ENU_adjust = (1 / np.sqrt(2)) *np.array([np.matmul(DCM_adjust[i], E_rkt[i]) for i in range(len(E_rkt))]) - vxB_ENU

        # Update the EFI Data
        line1_adjust.set_ydata(E_ENU_adjust[:, 0])
        line2_adjust.set_ydata(E_ENU_adjust[:, 1])
        line3_adjust.set_ydata(E_ENU_adjust[:, 2])
        fig.canvas.draw_idle()

    slope_slider.on_changed(update)
    offset_slider.on_changed(update)
    plt.show()






# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
test2_EFI_examine_DCM_timeOFFsets(wRocket)

