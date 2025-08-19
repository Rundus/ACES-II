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
Jpara_scale = 1E-9
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
data_dict_EEPAA_current_high = stl.loadDictFromFile(glob('C:\Data\ACESII\science\ESA_currents\high\*.cdf*')[0])
data_dict_EEPAA_current_low = stl.loadDictFromFile(glob('C:\Data\ACESII\science\ESA_currents\high\*.cdf*')[0])

# data_dict_IEPAA_current_high = stl.loadDictFromFile(glob('C:\Data\ACESII\science\ESA_currents\high\*.cdf*')[0])
# data_dict_IEPAA_current_low = stl.loadDictFromFile(glob('C:\Data\ACESII\science\ESA_currents\high\*.cdf*')[0])

# L-Shell Data
data_dict_Lshell_low = stl.loadDictFromFile('C:\Data\ACESII\coordinates\Lshell\low\ACESII_36364_Lshell.cdf')
data_dict_Lshell_high = stl.loadDictFromFile('C:\Data\ACESII\coordinates\Lshell\high\ACESII_36359_Lshell.cdf')

# DERPA Temperature Data
data_dict_DERPA_high = stl.loadDictFromFile('C:\Data\ACESII\L2\high\ACESII_36359_l2_ERPA.cdf')
data_dict_DERPA_low = stl.loadDictFromFile('C:\Data\ACESII\L2\low\ACESII_36364_l2_ERPA.cdf')

# Langmuir Probe Temperature Data
data_dict_LP_low = stl.loadDictFromFile(glob('C:\Data\ACESII\L3\Langmuir\low\*langmuir_fixed*')[0])
data_dict_LP_high = stl.loadDictFromFile(glob('C:\Data\ACESII\L3\Langmuir\high\*langmuir_fixed*')[0])

# Load MPI Data
# data_dict_MPI_low = stl.loadDictFromFile(glob('C:\Data\ACESII\L3\Langmuir\high\*langmuir_fixed*')[0])

stl.Done(start_time)


##########################
# --- --- --- --- --- ---
# --- PREPARE THE DATA ---
# --- --- --- --- --- ---
##########################
stl.prgMsg('Down-sampling Data')
# [1] --- Downsample data to only the specified L-Shell range ---

# Get the simulation range from the spatial environment data
data_dict_sim_spatial = stl.loadDictFromFile('C:\Data\physicsModels\ionosphere\spatial_environment\spatial_environment.cdf')
simLShell_min = deepcopy(data_dict_sim_spatial['simLShell'][0][0])
simLShell_max = deepcopy(data_dict_sim_spatial['simLShell'][0][-1])

# HIGH FLYER: reduce the Parallel Current data to the simulation size
low_idx = np.abs(data_dict_Lshell_high['L-Shell'][0] - simLShell_min).argmin()
high_idx = np.abs(data_dict_Lshell_high['L-Shell'][0] - simLShell_max).argmin()

for key in data_dict_EEPAA_current_high.keys():
    if key in ['Epoch', 'J_parallel']:
        data_dict_EEPAA_current_high[key][0] = data_dict_EEPAA_current_high[key][0][low_idx:high_idx+1]

# HIGH FLYER: reduce the swept Langmuir Probe Data
# TODO:


# HIGH FLYER: reduce the DERPA Data to the simulation size
low_idx = np.abs(data_dict_DERPA_high['Epoch1_temp'][0] - simLShell_min).argmin()
high_idx = np.abs(data_dict_DERPA_high['Epoch1_temp'][0] - simLShell_max).argmin()
for key in data_dict_Lshell_high.keys():
    if key in ['Epoch1_temp','']
# TODO:

# HIGH FLYER: reduce the L-Shell dataset itself
low_idx = np.abs(data_dict_Lshell_high['L-Shell'][0] - simLShell_min).argmin()
high_idx = np.abs(data_dict_Lshell_high['L-Shell'][0] - simLShell_max).argmin()
for key in data_dict_Lshell_high.keys():
    if key in ['Epoch', 'L-Shell', 'Alt']:
        data_dict_Lshell_high[key][0] = data_dict_Lshell_high[key][0][low_idx:high_idx+1]




# HIGH FLYER: interpolate L-Shell into Parallel Current data
Epoch_Jparallel_high_T0 = stl.EpochTo_T0_Rocket(data_dict_EEPAA_current_high['Epoch'][0], T0=dt.datetime(2022,11,20,17,20,00))
Epoch_Lshell_high_T0 = stl.EpochTo_T0_Rocket(data_dict_Lshell_high['Epoch'][0], T0=dt.datetime(2022,11,20,17,20,00))
cs = CubicSpline(Epoch_Lshell_high_T0, data_dict_Lshell_high['L-Shell'][0])
data_dict_EEPAA_current_high = {**data_dict_EEPAA_current_high, **{'L-Shell':[cs(Epoch_Jparallel_high_T0),{}]}}

# HIGH FLYER: Swept Langmuir Probe
# TODO:

# HIGH FLYER: interpolate L-Shell into DERPA data
Epoch_DERPA_high_T0 = stl.EpochTo_T0_Rocket(data_dict_DERPA_high['Epoch1_temp'][0], T0=dt.datetime(2022,11,20,17,20,00))
data_dict_DERPA_high = {**data_dict_DERPA_high, **{'L-Shell':[cs(Epoch_DERPA_high_T0), {}]}}

fig, ax = plt.subplots()
plt.plot(data_dict_DERPA_high['L-Shell'][0],data_dict_DERPA_high['ERPA1_temp'][0])
plt.show()







# LOW FLYER: reduce the EEPAA data
low_idx = np.abs(data_dict_Lshell_low['L-Shell'][0] - simLShell_min).argmin()
high_idx = np.abs(data_dict_Lshell_low['L-Shell'][0] - simLShell_max).argmin()
for key in data_dict_EEPAA_current_low.keys():
    if key in ['Epoch', 'J_parallel']:
        data_dict_EEPAA_current_low[key][0] = data_dict_EEPAA_current_low[key][0][low_idx:high_idx+1]

# LOW FLYER: reduce the Swept Langmuir Probe Data
# TODO:

# LOW FLYER: reduce the L-Shell data
low_idx = np.abs(data_dict_Lshell_low['L-Shell'][0] - simLShell_min).argmin()
high_idx = np.abs(data_dict_Lshell_low['L-Shell'][0] - simLShell_max).argmin()
for key in data_dict_Lshell_low.keys():
    if key in ['Epoch', 'L-Shell','Alt']:
        data_dict_Lshell_low[key][0] = data_dict_Lshell_low[key][0][low_idx:high_idx+1]
stl.Done(start_time)

# LOW FLYER: interpolate L-Shell into EEPAA/IEPAA currents data
Epoch_Jparallel_low_T0 = stl.EpochTo_T0_Rocket(data_dict_EEPAA_current_low['Epoch'][0], T0=dt.datetime(2022,11,20,17,20,00))
Epoch_Lshell_low_T0 = stl.EpochTo_T0_Rocket(data_dict_Lshell_low['Epoch'][0], T0=dt.datetime(2022,11,20,17,20,00))
cs = CubicSpline(Epoch_Lshell_low_T0, data_dict_Lshell_low['L-Shell'][0])
data_dict_EEPAA_current_low = {**data_dict_EEPAA_current_low, **{'L-Shell':[cs(Epoch_Jparallel_low_T0), {}]}}


# LOW FLYER: interpolate L-Shell into DERPA data
Epoch_DERPA_low_T0 = stl.EpochTo_T0_Rocket(data_dict_DERPA_low['Epoch1_temp'][0], T0=dt.datetime(2022,11,20,17,20,00))
data_dict_DERPA_low = {**data_dict_DERPA_low, **{'L-Shell':[cs(Epoch_DERPA_low_T0), {}]}}


# MPI DATA



############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

stl.prgMsg('Plotting Data')

fig, ax = plt.subplots(6, height_ratios=[1, 1, 0.5, 1, 1, 1], sharex=False)
fig.set_figwidth(Figure_width)
fig.set_figheight(Figure_height)

### HIGH FLYER DATA ###

# --- HF EEPAA/IEPAA PARALLEL CURRENT---
axNo = 0
ax[axNo].plot(data_dict_Lshell_high['L-Shell'][0], data_dict_EEPAA_current_high['J_parallel'][0]/Jpara_scale, label='EEPAA/IEPAA $J_{\parallel}$')
ax[axNo].set_ylabel('Parallel Current [nA]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
ax[axNo].set_ylim(-1E-6, 1E-7)

# --- HF Electron Temperatue (DERPA/LP) ---
axNo +=1
ax[axNo].plot(data_dict_DERPA_high['L-Shell'][0],data_dict_DERPA_high['ERPA1_temp'][0], color='tab:red', linewidth=Plot_LineWidth+2, label='ERPA1 $T_{e}$')
ax[axNo].plot(data_dict_DERPA_high['L-Shell'][0],data_dict_DERPA_high['ERPA2_temp'][0], color='tab:blue', linewidth=Plot_LineWidth+2, label='ERPA2 $T_{e}$')
ax[axNo].set_ylabel('Temperature [eV]', fontsize=Label_FontSize,labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='both', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
leg=ax[axNo].legend(fontsize=Legend_FontSize,loc='upper right')

# --- BREAK AXIS ---
axNo +=1
ax[axNo].axis('off')


### LOW FLYER DATA ###

# --- LF EEPAA/IEPAA PARALLEL CURRENT---
axNo +=1
ax[axNo].plot(data_dict_Lshell_low['L-Shell'][0], data_dict_EEPAA_current_low['J_parallel'][0]/Jpara_scale, label='EEPAA/IEPAA $J_{\parallel}$')
ax[axNo].set_ylabel('Parallel Current [nA]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
ax[axNo].set_ylim(-500,100)

# --- LF Electron Temperatue (DERPA/LP) ---
axNo +=1
ax[axNo].plot(data_dict_DERPA_low['L-Shell'][0],data_dict_DERPA_low['ERPA1_temp'][0], color='tab:red', linewidth=Plot_LineWidth+2, label='ERPA1 $T_{e}$')
ax[axNo].plot(data_dict_DERPA_low['L-Shell'][0],data_dict_DERPA_low['ERPA2_temp'][0], color='tab:blue', linewidth=Plot_LineWidth+2, label='ERPA2 $T_{e}$')
ax[axNo].set_ylabel('Temperature [eV]', fontsize=Label_FontSize,labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='both', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
leg=ax[axNo].legend(fontsize=Legend_FontSize,loc='upper right')

# ---Low Flyer MPI ---
axNo +=1
ax[axNo].plot(data_dict_DERPA_low['L-Shell'][0],data_dict_DERPA_low['ERPA2_temp'][0], color='tab:blue', linewidth=Plot_LineWidth+2, label='ERPA2 $T_{e}$')
ax[axNo].set_ylabel('Temp', fontsize=Label_FontSize,labelpad=Label_Padding)



# # --- get L-Shell labels and Alt Labels together ---
xTickLabels = ax[1].axes.get_xticklabels()
print(xTickLabels)
xTick_Locations = [float(tickVal.get_text()) for tickVal in xTickLabels]
xTick_newLabels_high = [f'{LshellVal}\n{round(data_dict_Lshell_high["Alt"][0][np.abs(data_dict_Lshell_high["L-Shell"][0] - LshellVal).argmin()]/stl.m_to_km)}' for LshellVal in xTick_Locations]
xTick_newLabels_low = [f'{LshellVal}\n{round(data_dict_Lshell_low["Alt"][0][np.abs(data_dict_Lshell_low["L-Shell"][0] - LshellVal).argmin()]/stl.m_to_km)}' for LshellVal in xTick_Locations]
ax[1].set_xticks(xTick_Locations,labels=xTick_newLabels_high)
ax[5].set_xticks(xTick_Locations,labels=xTick_newLabels_low)

# -- Do some minor adjustments to labels/margins/limits ---
for i in range(6):
    ax[i].margins(0)
    ax[i].set_xlim(simLShell_min, simLShell_max)
fig.align_ylabels(ax[:])

for i in [0, 1, 3, 4, 5]:
    ax[i].grid(alpha=0.9, which='both')

fig.subplots_adjust(left=0.13, bottom=0.06, right=0.88, top=0.98,hspace=0)  # remove the space between plots
plt.savefig(rf'{path_to_plots}\Plot2\Plot2_rkt_derived_data.png', dpi=dpi)
stl.Done(start_time)






