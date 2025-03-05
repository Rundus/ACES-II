# --- Langmuir_ni_upleg_downleg---
# Use the ACES-II LP data to analyze upleg/downleg n_i data from the
# payloads and try to estimate a nominal n_e vs altitude for the whole flight

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np


from src.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---


#################
# --- TOGGLES ---
#################
figure_width, figure_height = 10, 13
plot_LineWidth = 3
plot_textFontSize = 20
plot_MarkerSize = 20
plot_Colors = ['tab:blue', 'tab:red', 'tab:orange', 'tab:green', 'tab:purple', 'tab:olive', 'black']
text_FontSize = 25
title_FontSize = 14
label_fontsize = 22
labels_subplot_fontsize = 19
legend_FontSize = 15
legend_SubAxes_FontSize = 15
tick_LabelSize = 15
tick_SubplotLabelSize = 15
tick_Width = 2
tick_Length = 4
cbar_FontSize = 15
dpi = 100

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4
bool_make_plot_of_upleg_downleg = True

# --- OutputData ---
outputData = False
output_folder = r'C:\Users\cfelt\Desktop\Research\ACESII\Feltman2025_ACESII_JouleHeating\PLOTS\SCIENCE\density_upleg_downleg\\'

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.mission_attributes import *
from src.data_paths import DataPaths



def langmuir_upleg_downleg(wflyer):


    # --- load in the payload data ---
    stl.prgMsg('Loading Data')
    rocket_ID = ACESII.payload_IDs[wflyer]
    data_dict_LP = stl.loadDictFromFile(DataPaths.ACES_data_folder + f'\\L3\\Langmuir\\{ACESII.fliers[wflyer]}\\ACESII_{rocket_ID}_l3_langmuir_fixed.cdf')
    data_dict_eepaa = stl.loadDictFromFile(DataPaths.ACES_data_folder + f'\\L2\\{ACESII.fliers[wflyer]}\\ACESII_{rocket_ID}_l2_eepaa_fullCal.cdf')
    data_dict_attitude = stl.loadDictFromFile(DataPaths.ACES_data_folder + f'\\attitude\\{ACESII.fliers[wflyer]}\\ACESII_{rocket_ID}_Attitude_Solution.cdf')

    # --- determine the point of apogee  ---
    apogee_time = data_dict_attitude['Epoch'][0][np.argmax(data_dict_attitude['Alt'][0])]
    stl.Done(start_time)

    # --- split the data into upleg and downleg ---
    split_idx = np.abs(data_dict_LP['Epoch'][0] - apogee_time).argmin()
    ni_upleg = data_dict_LP['ni'][0][0:split_idx]
    ni_downleg = data_dict_LP['ni'][0][split_idx:]
    data_dict_attitude_interp = stl.InterpolateDataDict(InputDataDict=data_dict_attitude,
                                                        InputEpochArray=data_dict_attitude['Epoch'][0],
                                                        wKeys=['Alt', 'Lat'],
                                                        targetEpochArray=data_dict_LP['Epoch'][0])

    lat_upleg = data_dict_attitude_interp['Lat'][0][0:split_idx]
    alt_upleg = data_dict_attitude_interp['Alt'][0][0:split_idx]

    lat_downleg = data_dict_attitude_interp['Lat'][0][split_idx:]
    alt_downleg = data_dict_attitude_interp['Alt'][0][split_idx:]


    # --- indicate the regions of auroral activity vs. quiet ---
    auroral_start = [dt.datetime(2022,11,20,17,22,4), dt.datetime(2022,11,20,17,28,52,612000)]



    if bool_make_plot_of_upleg_downleg:

        stl.prgMsg('Making Plot of Upleg/Downleg')

        # --- make the plot ---
        fig, ax = plt.subplots()
        fig.set_size_inches(figure_width, figure_height)
        ax.plot(ni_upleg, alt_upleg/stl.m_to_km, label='UpLeg', color='tab:blue')
        ax.plot(ni_downleg, alt_downleg / stl.m_to_km, label='DownLeg', color='tab:red')
        ax.set_ylabel('Altitude [km]', fontsize=label_fontsize)
        ax.set_ylim(0, 430)
        ax.set_xlim(1E3, 1E7)
        ax.set_xscale('log')
        ax.set_xlabel('$n_{i}$ [cm$^{-3}$]', fontsize=label_fontsize)
        ax.tick_params(axis='y', which='major', labelsize=tick_LabelSize, width=tick_Width, length=tick_Length)
        ax.tick_params(axis='y', which='minor', labelsize=tick_LabelSize, width=tick_Width, length=tick_Length / 2)
        ax.tick_params(axis='x', which='major', labelsize=tick_LabelSize + 5, width=tick_Width, length=tick_Length)
        ax.tick_params(axis='x', which='minor', labelsize=tick_LabelSize + 5, width=tick_Width, length=tick_Length / 2)
        ax.legend(fontsize=legend_FontSize)
        plt.tight_layout()
        plt.savefig(output_folder + f'ACESII_LP_{rocket_ID}_upLeg_downLeg.png',dpi=dpi)

        stl.Done(start_time)

# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
langmuir_upleg_downleg(wRocket-4)