# --- Plot2_rkt_data_stackplot.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Stack Plot detailing the derived quantities from the rocket instruments


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np
import spaceToolsLib as stl
from src.my_imports import *

start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from src.Papers.ACESII_Joule_Heating.file_toggles import *
print(stl.color.UNDERLINE + f'Plot2_rkt_derived_data' + stl.color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
plt.rcParams["font.family"] = "Arial"
dpi = 800
Jpara_scale = 1E-6
LP_scale = 1

# --- Cbar ---
cbarMin, cbarMax = 1E8, 1E11
cbar_TickLabelSize = 14
# LP_limit = [0,2.75]
LP_limit = [1E4, 3E5]
my_cmap = stl.apl_rainbow_black0_cmap()
my_cmap.set_bad(color=(0,0,0))

# --- Plot toggles ---
Figure_width = 13 # in inches
Figure_height =15# in inches
Text_FontSize = 20
Label_FontSize = 20
Tick_FontSize = 20
Tick_Length = 5
Tick_Width = 2
Plot_LineWidth = 1.5
Label_Padding = 15
Tick_Padding = 10
Legend_FontSize = 15
cbar_FontSize = 25

# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---
stl.prgMsg('Loading Data')


# EEPAA/IEPAA Parallel Current -
data_dict_EEPAA_current_high = stl.loadDictFromFile(glob('C:\Data\ACESII\science\ESA_currents\high\*eepaa.cdf*')[0])
data_dict_EEPAA_current_low = stl.loadDictFromFile(glob('C:\Data\ACESII\science\ESA_currents\low\*eepaa.cdf*')[0])
data_dict_IEPAA_current_high = stl.loadDictFromFile(glob('C:\Data\ACESII\science\ESA_currents\high\*iepaa.cdf*')[0])
data_dict_IEPAA_current_low = stl.loadDictFromFile(glob('C:\Data\ACESII\science\ESA_currents\low\*iepaa.cdf*')[0])

# L-Shell Data
data_dict_Lshell_low = stl.loadDictFromFile('C:\Data\ACESII\coordinates\Lshell\low\ACESII_36364_Lshell.cdf')
data_dict_Lshell_high = stl.loadDictFromFile('C:\Data\ACESII\coordinates\Lshell\high\ACESII_36359_Lshell.cdf')

# DERPA Temperature Data
data_dict_DERPA1_high = stl.loadDictFromFile('C:\Data\ACESII\L3\DERPA\high\ACESII_36359_l3_ERPA1_fullCal.cdf')
data_dict_DERPA2_high = stl.loadDictFromFile('C:\Data\ACESII\L3\DERPA\high\ACESII_36359_l3_ERPA2_fullCal.cdf')
data_dict_DERPA1_low = stl.loadDictFromFile('C:\Data\ACESII\L3\DERPA\low\ACESII_36364_l3_ERPA1_fullCal.cdf')
data_dict_DERPA2_low = stl.loadDictFromFile('C:\Data\ACESII\L3\DERPA\low\ACESII_36364_l3_ERPA2_fullCal.cdf')

# EISCAT calibration data
data_dict_EISCAT_cal_high = stl.loadDictFromFile('C:\Data\ACESII\calibration\LP_postFlight_calibration\high\ACESII_36359_postFlight_cal.cdf')
data_dict_EISCAT_cal_low = stl.loadDictFromFile('C:\Data\ACESII\calibration\LP_postFlight_calibration\low\ACESII_36364_postFlight_cal.cdf')

# Langmuir Probe Temperature Data
data_dict_LP_low = stl.loadDictFromFile(glob('C:\Data\ACESII\L3\Langmuir\low\*langmuir_fixed_fullCal*')[0])
data_dict_LP_high = stl.loadDictFromFile(glob('C:\Data\ACESII\L3\Langmuir\high\*langmuir_fixed_fullCal*')[0])

# Load MPI Data
# data_dict_MPI_low = stl.loadDictFromFile(glob('C:\Data\ACESII\L3\Langmuir\high\*langmuir_fixed*')[0])

stl.Done(start_time)


##########################
# --- --- --- --- --- ---
# --- PREPARE THE DATA ---
# --- --- --- --- --- ---
##########################
# [1] --- Downsample data to only the specified L-Shell range ---
stl.prgMsg('Down-sampling Data')

# Get the simulation range from the spatial environment data
data_dict_sim_spatial = stl.loadDictFromFile('C:\Data\physicsModels\ionosphere\spatial_environment\spatial_environment.cdf')
simLShell_min = deepcopy(data_dict_sim_spatial['simLShell'][0][0])
simLShell_max = deepcopy(data_dict_sim_spatial['simLShell'][0][-1])


####################
# --- HIGH FLYER ---
####################

# TODO: Add the swept Langmuir data
data_dicts_high = [data_dict_EEPAA_current_high, data_dict_IEPAA_current_high, data_dict_DERPA1_high, data_dict_DERPA2_high]

# --- Interpolate L-Shell into data ---
Epoch_Lshell_high_T0 = stl.EpochTo_T0_Rocket(data_dict_Lshell_high['Epoch'][0], T0=dt.datetime(2022, 11, 20, 17, 20, 00))
cs = CubicSpline(Epoch_Lshell_high_T0, data_dict_Lshell_high['L-Shell'][0])

for idx, ddict in enumerate(data_dicts_high):
    Epoch_T0 = stl.EpochTo_T0_Rocket(ddict['Epoch'][0], T0=dt.datetime(2022, 11, 20, 17, 20, 00))
    data_dicts_high[idx] = {**data_dicts_high[idx], **{'L-Shell': [cs(Epoch_T0), {}]}}

# SPECIAL: Interpolate Ion J_para data onto EEPAA Epoch
Epoch_T0_iepaa = stl.EpochTo_T0_Rocket(deepcopy(data_dicts_high[1]['Epoch'][0]), T0=dt.datetime(2022, 11, 20, 17, 20, 00))
Epoch_T0_eepaa = stl.EpochTo_T0_Rocket(deepcopy(data_dicts_high[0]['Epoch'][0]), T0=dt.datetime(2022, 11, 20, 17, 20, 00))
keys = ['j_para', 'L-Shell']
for key in keys:
    cs = CubicSpline(Epoch_T0_iepaa, deepcopy(data_dicts_high[1][f'{key}'][0]))
    data_dicts_high[1][f'{key}'][0] = cs(Epoch_T0_eepaa)
data_dicts_high[1]['Epoch'][0] = deepcopy(data_dicts_high[1]['Epoch'][0])

# --- Reduce the Data to the Simulation Size ---
for idx, ddict in enumerate(data_dicts_high):
    if idx == 0:
        low_idx = np.abs(data_dict_EEPAA_current_high['L-Shell'][0] - simLShell_min).argmin()
        high_idx = np.abs(data_dict_EEPAA_current_high['L-Shell'][0] - simLShell_max).argmin()
        keys = ['Epoch', 'j_para', 'L-Shell']
    elif idx == 1:
        low_idx = np.abs(data_dict_IEPAA_current_high['L-Shell'][0] - simLShell_min).argmin()
        high_idx = np.abs(data_dict_IEPAA_current_high['L-Shell'][0] - simLShell_max).argmin()
        keys = ['Epoch', 'j_para', 'L-Shell']
    elif idx == 2:
        low_idx = np.abs(data_dict_DERPA1_high['L-Shell'][0] - simLShell_min).argmin()
        high_idx = np.abs(data_dict_DERPA1_high['L-Shell'][0] - simLShell_max).argmin()
        keys = ['Epoch', 'temperature','L-Shell']
    elif idx == 3:
        low_idx = np.abs(data_dict_DERPA2_high['L-Shell'][0] - simLShell_min).argmin()
        high_idx = np.abs(data_dict_DERPA2_high['L-Shell'][0] - simLShell_max).argmin()
        keys = ['Epoch', 'temperature','L-Shell']
    # elif idx == 3:
    #     continue
    for key in keys:
        ddict[key][0] = ddict[key][0][low_idx:high_idx+1]
data_dict_EEPAA_current_high, data_dict_IEPAA_current_high, data_dict_DERPA1_high,data_dict_DERPA2_high = deepcopy(data_dicts_high[0]),deepcopy(data_dicts_high[1]),deepcopy(data_dicts_high[2]),deepcopy(data_dicts_high[3])

###################
# --- LOW FLYER ---
###################

data_dicts_low = [data_dict_EEPAA_current_low, data_dict_IEPAA_current_low, data_dict_DERPA1_low, data_dict_DERPA2_low]

# --- Interpolate L-Shell into data ---
Epoch_Lshell_low_T0 = stl.EpochTo_T0_Rocket(data_dict_Lshell_low['Epoch'][0], T0=dt.datetime(2022, 11, 20, 17, 20, 00))
cs = CubicSpline(Epoch_Lshell_low_T0, data_dict_Lshell_low['L-Shell'][0])

for idx,ddict in enumerate(data_dicts_low):
    Epoch_T0 = stl.EpochTo_T0_Rocket(data_dicts_low[idx]['Epoch'][0], T0=dt.datetime(2022, 11, 20, 17, 20, 00))
    data_dicts_low[idx] = {**data_dicts_low[idx], **{'L-Shell': [cs(Epoch_T0), {}]}}

# SPECIAL: Interpolate Ion J_para data onto EEPAA Epoch
Epoch_T0_iepaa = stl.EpochTo_T0_Rocket(deepcopy(data_dicts_low[1]['Epoch'][0]), T0=dt.datetime(2022, 11, 20, 17, 20, 00))
Epoch_T0_eepaa = stl.EpochTo_T0_Rocket(deepcopy(data_dicts_low[0]['Epoch'][0]), T0=dt.datetime(2022, 11, 20, 17, 20, 00))
keys = ['j_para', 'L-Shell']
for key in keys:
    cs = CubicSpline(Epoch_T0_iepaa, deepcopy(data_dicts_low[1][f'{key}'][0]))
    data_dicts_low[1][f'{key}'][0] = cs(Epoch_T0_eepaa)
data_dicts_low[1]['Epoch'][0] = deepcopy(data_dicts_low[0]['Epoch'][0])


# --- Reduce the Data to the Simulation Size ---
for idx, ddict in enumerate(data_dicts_low):
    if idx == 0:
        low_idx = np.abs(data_dict_EEPAA_current_low['L-Shell'][0] - simLShell_min).argmin()
        high_idx = np.abs(data_dict_EEPAA_current_low['L-Shell'][0] - simLShell_max).argmin()
        keys = ['Epoch', 'j_para', 'L-Shell']
    elif idx == 1:
        low_idx = np.abs(data_dict_IEPAA_current_low['L-Shell'][0] - simLShell_min).argmin()
        high_idx = np.abs(data_dict_IEPAA_current_low['L-Shell'][0] - simLShell_max).argmin()
        keys = ['Epoch', 'j_para', 'L-Shell']
    elif idx == 2:
        low_idx = np.abs(data_dict_DERPA1_low['L-Shell'][0] - simLShell_min).argmin()
        high_idx = np.abs(data_dict_DERPA1_low['L-Shell'][0] - simLShell_max).argmin()
        keys = ['Epoch', 'temperature','L-Shell']
    elif idx == 2:
        low_idx = np.abs(data_dict_DERPA2_low['L-Shell'][0] - simLShell_min).argmin()
        high_idx = np.abs(data_dict_DERPA2_low['L-Shell'][0] - simLShell_max).argmin()
        keys = ['Epoch', 'temperature','L-Shell']
    # elif idx == 3:
    #     continue
    for key in keys:
        ddict[key][0] = ddict[key][0][low_idx:high_idx+1]
data_dict_EEPAA_current_low, data_dict_IEPAA_current_low, data_dict_DERPA1_low, data_dict_DERPA2_low = deepcopy(data_dicts_low[0]),deepcopy(data_dicts_low[1]),deepcopy(data_dicts_low[2]),deepcopy(data_dicts_low[3])
stl.Done(start_time)




############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

stl.prgMsg('Plotting Data')

fig, ax = plt.subplots(6, height_ratios=[1, 1, 0.35, 1, 1, 1], sharex=False)
fig.set_figwidth(Figure_width)
fig.set_figheight(Figure_height)

### HIGH FLYER DATA ###

# --- HF EEPAA/IEPAA PARALLEL CURRENT---
axNo = 0
ax[axNo].plot(data_dict_EEPAA_current_high['L-Shell'][0], data_dict_EEPAA_current_high['j_para'][0]/Jpara_scale, label='EEPAA $J_{\parallel}$',color='tab:blue')
ax[axNo].plot(data_dict_IEPAA_current_high['L-Shell'][0], data_dict_IEPAA_current_high['j_para'][0]/Jpara_scale, label='IEPAA $J_{\parallel}$',color='tab:red')
ax[axNo].legend(fontsize=Legend_FontSize)
# ax[axNo].plot(data_dict_EEPAA_current_high['L-Shell'][0], (data_dict_IEPAA_current_high['j_para'][0]+data_dict_EEPAA_current_high['j_para'][0])/Jpara_scale, label='EEPAA $J_{\parallel}$',color='black')
ax[axNo].set_ylabel('Parallel Current\n[$\mu$A]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
ax[axNo].set_ylim(-7, 12)


# --- HF Electron Temperatue (DERPA/LP) ---
axNo +=1
ax[axNo].plot(data_dict_EISCAT_cal_high['L-Shell'][0], data_dict_EISCAT_cal_high['Te'][0],color='tab:purple', linewidth=Plot_LineWidth, label='EISCAT Background $T_{e}$')
ax[axNo].plot(data_dict_DERPA1_high['L-Shell'][0],data_dict_DERPA1_high['temperature'][0], color='tab:red', linewidth=Plot_LineWidth, label='ERPA1 $T_{e}$')
ax[axNo].plot(data_dict_DERPA2_high['L-Shell'][0],data_dict_DERPA2_high['temperature'][0], color='tab:blue', linewidth=Plot_LineWidth, label='ERPA2 $T_{e}$')
ax[axNo].set_ylabel('Temperature\n[eV]', fontsize=Label_FontSize,labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='both', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
leg=ax[axNo].legend(fontsize=Legend_FontSize,loc='upper right')
ax[axNo].set_ylim(0.05,0.5)



# --- BREAK AXIS ---
axNo +=1
ax[axNo].axis('off')


### LOW FLYER DATA ###

# --- LF EEPAA/IEPAA PARALLEL CURRENT---
axNo +=1
ax[axNo].plot(data_dict_EEPAA_current_low['L-Shell'][0], data_dict_EEPAA_current_low['j_para'][0]/Jpara_scale, label='EEPAA $J_{\parallel}$',color='tab:blue')
ax[axNo].plot(data_dict_IEPAA_current_low['L-Shell'][0], data_dict_IEPAA_current_low['j_para'][0]/Jpara_scale, label='IEPAA $J_{\parallel}$',color='tab:red')
ax[axNo].legend(fontsize=Legend_FontSize)
# ax[axNo].plot(data_dict_EEPAA_current_low['L-Shell'][0], (data_dict_IEPAA_current_low['j_para'][0]+data_dict_EEPAA_current_low['j_para'][0])/Jpara_scale, label='EEPAA $J_{\parallel}$',color='black')
ax[axNo].set_ylabel('Parallel Current\n[$\mu$A]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
ax[axNo].set_ylim(-7, 7)

# --- LF Electron Temperatue (DERPA/LP) ---
axNo +=1
ax[axNo].plot(data_dict_EISCAT_cal_low['L-Shell'][0], data_dict_EISCAT_cal_low['Te'][0],color='tab:purple', linewidth=Plot_LineWidth, label='EISCAT Background $T_{e}$')
ax[axNo].plot(data_dict_DERPA1_low['L-Shell'][0],data_dict_DERPA1_low['temperature'][0], color='tab:red', linewidth=Plot_LineWidth, label='ERPA1 $T_{e}$')
ax[axNo].plot(data_dict_DERPA2_low['L-Shell'][0],data_dict_DERPA2_low['temperature'][0], color='tab:blue', linewidth=Plot_LineWidth, label='ERPA2 $T_{e}$')
ax[axNo].set_ylabel('Temperature\n[eV]', fontsize=Label_FontSize,labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='both', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
leg=ax[axNo].legend(fontsize=Legend_FontSize,loc='upper right')
ax[axNo].set_ylim(0.05, 0.5)

# ---Low Flyer MPI ---
axNo +=1
ax[axNo].text(x=0.5,y=0.5,s='MPI Data', fontsize=Text_FontSize*2, transform=ax[axNo].transAxes)
# ax[axNo].plot(data_dict_DERPA_low['L-Shell'][0],data_dict_DERPA_low['temperature'][0], color='tab:blue', linewidth=Plot_LineWidth+2, label='ERPA2 $T_{e}$')
# ax[axNo].set_ylabel('Temp', fontsize=Label_FontSize,labelpad=Label_Padding)

# # --- get L-Shell labels and Alt Labels together ---
xTickLabels = ax[1].axes.get_xticklabels()
xTick_Locations = [tickVal.get_text() for tickVal in xTickLabels]
xTick_Locations = [float(tickVal.get_text()) for tickVal in xTickLabels]
xTick_newLabels_high = [f'{LshellVal}\n{round(data_dict_Lshell_high["Alt"][0][np.abs(data_dict_Lshell_high["L-Shell"][0] - LshellVal).argmin()]/stl.m_to_km)}' for LshellVal in xTick_Locations]
xTick_newLabels_low = [f'{LshellVal}\n{round(data_dict_Lshell_low["Alt"][0][np.abs(data_dict_Lshell_low["L-Shell"][0] - LshellVal).argmin()]/stl.m_to_km)}' for LshellVal in xTick_Locations]

ax[1].set_xticks(xTick_Locations, labels=xTick_newLabels_high)
ax[1].set_xlabel('L-Shell \n Alt [km]', fontsize=Tick_FontSize, weight='bold')
ax[1].xaxis.set_label_coords(-0.085, -0.03)


ax[5].set_xticks(xTick_Locations,labels=xTick_newLabels_low)
ax[5].set_xlabel('L-Shell \n Alt [km]', fontsize=Tick_FontSize, weight='bold')
ax[5].xaxis.set_label_coords(-0.085, -0.03)

for num in [1, 5]:
    ax[num].tick_params(axis='x', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width, pad=Tick_Padding)

# -- Do some minor adjustments to labels/margins/limits ---
for i in range(6):
    ax[i].margins(0)
    ax[i].set_xlim(simLShell_min, simLShell_max)
fig.align_ylabels(ax[:])

for i in [0, 1, 3, 4, 5]:
    ax[i].grid(alpha=0.9, which='both')

fig.subplots_adjust(left=0.13, bottom=0.06,right=0.99, top=0.98,hspace=0)  # remove the space between plots
plt.savefig(rf'{path_to_plots}\Plot2\Plot2_rkt_derived_data.png', dpi=dpi)
stl.Done(start_time)






