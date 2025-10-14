# --- Plot1_rkt_data_stackplot.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Stack Plot detailing the inputs into the simulation


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
print(stl.color.UNDERLINE + f'Plot1_data_stackplot' + stl.color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
plt.rcParams["font.family"] = "Arial"
dpi = 800
Escale = 1000 # what to scale the deltaE field by
# LP_scale = 1E5
LP_scale = 1
E_limits = [-125, 265]
B_limits = [-1500, 3000]

# --- Cbar ---
cbarMin, cbarMax = 1E8, 1E11
cbar_TickLabelSize = 14
# LP_limit = [0,2.75]
LP_limit = [1E4, 7E5]
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

# DC B-Field
data_dict_BField_high = stl.loadDictFromFile(r'C:\Data\ACESII\L2\high\ACESII_36359_RingCore_auroral_median_filter.cdf')
data_dict_BField_low = stl.loadDictFromFile(r'C:\Data\ACESII\L2\low\ACESII_36364_RingCore_auroral_median_filter.cdf')

# E-Field Data
data_dict_Efield_low = stl.loadDictFromFile('C:\Data\ACESII\L2\low\ACESII_36364_l2_EFI_auroral.cdf')

# EEPAA Particle Data
data_dict_flux_low = stl.loadDictFromFile(glob('C:\Data\ACESII\L3\Energy_Flux\low\ACESII_36364_l3_eepaa_flux.cdf')[0])
data_dict_flux_high = stl.loadDictFromFile(glob('C:\Data\ACESII\L3\Energy_Flux\high\ACESII_36359_l3_eepaa_flux.cdf')[0])

# L-Shell Data
data_dict_Lshell_low = stl.loadDictFromFile(glob('C:\Data\ACESII\coordinates\Lshell\low\ACESII_36364_Lshell.cdf')[0])
data_dict_Lshell_high = stl.loadDictFromFile(glob('C:\Data\ACESII\coordinates\Lshell\high\ACESII_36359_Lshell.cdf')[0])

# Langmuir Probe Density Data
data_dict_LP_low = stl.loadDictFromFile(glob('C:\Data\ACESII\L3\Langmuir\low\ACESII_36364_l3_langmuir_fixed_fullCal.cdf')[0])
data_dict_LP_high = stl.loadDictFromFile(glob('C:\Data\ACESII\L3\Langmuir\high\ACESII_36359_l3_langmuir_fixed_fullCal.cdf')[0])

# EISCAT calibration data
data_dict_EISCAT_cal_high = stl.loadDictFromFile('C:\Data\ACESII\calibration\LP_postFlight_calibration\high\ACESII_36359_postFlight_cal.cdf')
data_dict_EISCAT_cal_low = stl.loadDictFromFile('C:\Data\ACESII\calibration\LP_postFlight_calibration\low\ACESII_36364_postFlight_cal.cdf')
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

# HIGH FLYER: reduce the EEPAA data
low_idx = np.abs(data_dict_Lshell_high['L-Shell'][0] - simLShell_min).argmin()
high_idx = np.abs(data_dict_Lshell_high['L-Shell'][0] - simLShell_max).argmin()
for key in data_dict_flux_high.keys():
    if key in ['Epoch','varPhi_E_Parallel']:
        data_dict_flux_high[key][0] = data_dict_flux_high[key][0][low_idx:high_idx+1]

# HIGH FLYER: reduce the Langmuir Probe Data
time_target_low = data_dict_Lshell_high['Epoch'][0][low_idx]
time_target_high = data_dict_Lshell_high['Epoch'][0][high_idx]
low_idx = np.abs(data_dict_LP_high['Epoch'][0] - time_target_low).argmin()
high_idx = np.abs(data_dict_LP_high['Epoch'][0] - time_target_high).argmin()
for key in data_dict_LP_high.keys():
    if key in ['Epoch', 'ni']:
        data_dict_LP_high[key][0] = data_dict_LP_high[key][0][low_idx:high_idx+1]

# HIGH FLYER: remove some bad indices from the LP High Data
temp = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_LP_high['Epoch'][0]])
bad_idxs = np.where(temp < 0)
for key in data_dict_LP_high.keys():
    if key in ['Epoch','ni']:
        data_dict_LP_high[key][0] = np.delete(data_dict_LP_high[key][0],bad_idxs)

# HIGH FLYER: reduce the L-Shell data
low_idx = np.abs(data_dict_Lshell_high['L-Shell'][0] - simLShell_min).argmin()
high_idx = np.abs(data_dict_Lshell_high['L-Shell'][0] - simLShell_max).argmin()
for key in data_dict_Lshell_high.keys():
    if key in ['Epoch', 'L-Shell','Alt']:
        data_dict_Lshell_high[key][0] = data_dict_Lshell_high[key][0][low_idx:high_idx+1]

# HIGH FLYER: interpolate L-Shell into LP data
Epoch_LP_high_T0 = stl.EpochTo_T0_Rocket(data_dict_LP_high['Epoch'][0],T0=dt.datetime(2022,11,20,17,20,00))
Epoch_LShell_high_T0 = stl.EpochTo_T0_Rocket(data_dict_Lshell_high['Epoch'][0],T0=dt.datetime(2022,11,20,17,20,00))
cs = CubicSpline(Epoch_LShell_high_T0,data_dict_Lshell_high['L-Shell'][0])
data_dict_LP_high = {**data_dict_LP_high,
                     **{'L-Shell':[cs(Epoch_LP_high_T0),{}]}}

# HIGH FLYER: interpolate L-Shell into B-Field median data
Epoch_B_high_T0 = stl.EpochTo_T0_Rocket(data_dict_BField_high['Epoch'][0],T0=dt.datetime(2022,11,20,17,20,00))
cs = CubicSpline(Epoch_LShell_high_T0,data_dict_Lshell_high['L-Shell'][0])
data_dict_BField_high = {**data_dict_BField_high, **{'L-Shell':[cs(Epoch_B_high_T0),{}]}}

# HIGH FLYER: reduce B-Field data to only the relevant range
low_idx = np.abs(data_dict_BField_high['L-Shell'][0] - simLShell_min).argmin()
high_idx = np.abs(data_dict_BField_high['L-Shell'][0] - simLShell_max).argmin()
for key in data_dict_BField_high.keys():
    data_dict_BField_high[key][0] = data_dict_BField_high[key][0][low_idx:high_idx+1]

# LOW FLYER: reduce the EEPAA data
low_idx = np.abs(data_dict_Lshell_low['L-Shell'][0] - simLShell_min).argmin()
high_idx = np.abs(data_dict_Lshell_low['L-Shell'][0] - simLShell_max).argmin()
for key in data_dict_flux_low.keys():
    if key in ['Epoch','varPhi_E_Parallel']:
        data_dict_flux_low[key][0] = data_dict_flux_low[key][0][low_idx:high_idx+1]

# LOW FLYER: reduce the Langmuir Probe Data
time_target_low = data_dict_Lshell_low['Epoch'][0][low_idx]
time_target_high = data_dict_Lshell_low['Epoch'][0][high_idx]
low_idx = np.abs(data_dict_LP_low['Epoch'][0] - time_target_low).argmin()
high_idx = np.abs(data_dict_LP_low['Epoch'][0] - time_target_high).argmin()
for key in data_dict_LP_low.keys():
    if key in ['Epoch', 'ni']:
        data_dict_LP_low[key][0] = data_dict_LP_low[key][0][low_idx:high_idx+1]

# LOW FLYER: reduce the EFI Data
low_idx = np.abs(data_dict_Efield_low['Epoch'][0] - time_target_low).argmin()
high_idx = np.abs(data_dict_Efield_low['Epoch'][0] - time_target_high).argmin()
for key in data_dict_Efield_low.keys():
    if key in ['Epoch', 'E_N','E_T','E_p']:
        data_dict_Efield_low[key][0] = data_dict_Efield_low[key][0][low_idx:high_idx+1]

# LOW FLYER: reduce the L-Shell data
low_idx = np.abs(data_dict_Lshell_low['L-Shell'][0] - simLShell_min).argmin()
high_idx = np.abs(data_dict_Lshell_low['L-Shell'][0] - simLShell_max).argmin()
for key in data_dict_Lshell_low.keys():
    if key in ['Epoch', 'L-Shell','Alt']:
        data_dict_Lshell_low[key][0] = data_dict_Lshell_low[key][0][low_idx:high_idx+1]
stl.Done(start_time)

# LOW FLYER: interpolate L-Shell into LP data
Epoch_LP_low_T0 = stl.EpochTo_T0_Rocket(data_dict_LP_low['Epoch'][0],T0=dt.datetime(2022,11,20,17,20,00))
Epoch_LShell_low_T0 = stl.EpochTo_T0_Rocket(data_dict_Lshell_low['Epoch'][0],T0=dt.datetime(2022,11,20,17,20,00))
cs = CubicSpline(Epoch_LShell_low_T0,data_dict_Lshell_low['L-Shell'][0])
data_dict_LP_low = {**data_dict_LP_low,
                     **{'L-Shell':[cs(Epoch_LP_low_T0),{}]}}

# LOW FLYER: interpolate L-Shell into LP data
Epoch_EFI_low_T0 = stl.EpochTo_T0_Rocket(data_dict_Efield_low['Epoch'][0],T0=dt.datetime(2022,11,20,17,20,00))
Epoch_LShell_low_T0 = stl.EpochTo_T0_Rocket(data_dict_Lshell_low['Epoch'][0],T0=dt.datetime(2022,11,20,17,20,00))
cs = CubicSpline(Epoch_LShell_low_T0,data_dict_Lshell_low['L-Shell'][0])
data_dict_Efield_low = {**data_dict_Efield_low,
                     **{'L-Shell':[cs(Epoch_EFI_low_T0),{}]}}

# HIGH FLYER: interpolate L-Shell into B-Field median data
Epoch_B_low_T0 = stl.EpochTo_T0_Rocket(data_dict_BField_low['Epoch'][0],T0=dt.datetime(2022,11,20,17,20,00))
cs = CubicSpline(Epoch_LShell_low_T0,data_dict_Lshell_low['L-Shell'][0])
data_dict_BField_low = {**data_dict_BField_low, **{'L-Shell':[cs(Epoch_B_low_T0),{}]}}

# HIGH FLYER: reduce B-Field data to only the relevant range
low_idx = np.abs(data_dict_BField_low['L-Shell'][0] - simLShell_min).argmin()
high_idx = np.abs(data_dict_BField_low['L-Shell'][0] - simLShell_max).argmin()
for key in data_dict_BField_low.keys():
    data_dict_BField_low[key][0] = data_dict_BField_low[key][0][low_idx:high_idx+1]

############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

stl.prgMsg('Plotting Data')

fig, ax = plt.subplots(8, height_ratios=[1.5, 1, 0.5, 0.5, 1.5, 1, 1, 0.5],sharex=False)
fig.set_figwidth(Figure_width)
fig.set_figheight(Figure_height)


### HIGH FLYER DATA ###

# --- HF EEPAA 0 to 180 DEG---
axNo = 0
cmap = ax[axNo].pcolormesh(data_dict_Lshell_high['L-Shell'][0],data_dict_flux_high['Energy'][0],data_dict_flux_high['varPhi_E_Parallel'][0].T, cmap=my_cmap,vmin=cbarMin,vmax=cbarMax, norm='log')
ax[axNo].set_ylabel('Parallel Energy Flux\nEnergy [eV]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
ax[axNo].set_yscale('log')
ax[axNo].set_ylim(28,1E4)

# --- HF DC B-Field ---
axNo +=1
ax[axNo].plot(data_dict_BField_high['L-Shell'][0],data_dict_BField_high['B_p'][0], color='tab:blue', linewidth=Plot_LineWidth+2, label='$\Delta B_{p}$')
ax[axNo].plot(data_dict_BField_high['L-Shell'][0],data_dict_BField_high['B_T'][0], color='tab:red', linewidth=Plot_LineWidth+2, label='$\Delta B_{T}$')
ax[axNo].plot(data_dict_BField_high['L-Shell'][0],data_dict_BField_high['B_N'][0], color='tab:green', linewidth=Plot_LineWidth+2, label='$\Delta B_{N}$')
ax[axNo].set_ylabel('DC B-Field\n[nT]', fontsize=Label_FontSize,labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='both', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].set_ylim(B_limits)
leg=ax[axNo].legend(fontsize=Legend_FontSize,loc='upper right')
for line in leg.get_lines():
    line.set_linewidth(4)


# --- LP High---
axNo +=1
colorChoice = 'black'
filtered = stl.butterFilter().butter_filter(data= data_dict_LP_high['ni'][0]/LP_scale,
                                                              lowcutoff=0.03,
                                                              highcutoff=0.03,
                                                              fs=100,
                                                              filtertype='lowpass',
                                                              order=4
                                                              )

ax[axNo].plot(data_dict_EISCAT_cal_high['L-Shell'][0], data_dict_EISCAT_cal_high['ne'][0]/(np.power(stl.cm_to_m,3)),color='tab:purple', linewidth=Plot_LineWidth+1, label='EISCAT Background n$_{e}$')
ax[axNo].plot(data_dict_LP_high['L-Shell'][0], filtered/(np.power(stl.cm_to_m,3)),color=colorChoice, linewidth=Plot_LineWidth+1)
ax[axNo].set_ylabel('LP\n[cm$^{-3}$]', fontsize=Label_FontSize-2, color=colorChoice, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize-3, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize - 6, length=Tick_Length-2, width=Tick_Width)
ax[axNo].tick_params(axis='x', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length + 4, width=Tick_Width, pad=Tick_Padding)
ax[axNo].tick_params(axis='x', which='minor', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width, pad=Tick_Padding)
ax[axNo].set_ylim(LP_limit[0], LP_limit[1])
ax[axNo].minorticks_on()
ax[axNo].set_yscale('log')
ax[axNo].xaxis.set_tick_params(labelbottom=True)
ax[axNo].set_xlabel('L-Shell \n Alt [km]', fontsize=Tick_FontSize-2, weight='bold')
ax[axNo].xaxis.set_label_coords(-0.085, -0.26)
ax[axNo].legend()

# --- BREAK AXIS ---
axNo +=1
ax[axNo].axis('off')

# # --- LF EEPAA---
axNo +=1
cmap = ax[axNo].pcolormesh(data_dict_Lshell_low['L-Shell'][0], data_dict_flux_low['Energy'][0], data_dict_flux_low['varPhi_E_Parallel'][0].T, cmap=my_cmap, vmin=cbarMin, vmax=cbarMax, norm='log')
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
ax[axNo].set_ylabel('Parallel Energy Flux\nEnergy [eV]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].set_yscale('log')
ax[axNo].set_ylim(28,1E4)

# # --- DC Delta B LF---
axNo +=1
ax[axNo].plot(data_dict_BField_low['L-Shell'][0], data_dict_BField_low['B_p'][0], color='tab:blue', linewidth=Plot_LineWidth+2, label='$\Delta B_{p}$')
ax[axNo].plot(data_dict_BField_low['L-Shell'][0], data_dict_BField_low['B_T'][0], color='tab:red', linewidth=Plot_LineWidth+2, label='$\Delta B_{T}$')
ax[axNo].plot(data_dict_BField_low['L-Shell'][0], data_dict_BField_low['B_N'][0], color='tab:green', linewidth=Plot_LineWidth+2, label='$\Delta B_{N}$')
ax[axNo].set_ylabel('DC B-Field\n[nT]', fontsize=Label_FontSize,labelpad=Label_Padding)
ax[axNo].tick_params(axis='y',which='both', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
# ax[axNo].set_ylim(-700, 1000)
ax[axNo].set_ylim(-750, 1050)
leg=ax[axNo].legend(fontsize=Legend_FontSize,loc='upper right')
for line in leg.get_lines():
    line.set_linewidth(4)

# ---Low Flyer E-Field ---
axNo +=1
ax[axNo].plot(data_dict_Efield_low['L-Shell'][0], data_dict_Efield_low['E_p'][0]*Escale, linewidth=Plot_LineWidth, color='tab:blue', label=r'$E_{p}$')
ax[axNo].plot(data_dict_Efield_low['L-Shell'][0], data_dict_Efield_low['E_T'][0]*Escale, linewidth=Plot_LineWidth, color='tab:red', label=r'$E_{T}$')
ax[axNo].plot(data_dict_Efield_low['L-Shell'][0], data_dict_Efield_low['E_N'][0]*Escale, linewidth=Plot_LineWidth, color='tab:green', label=r'$E_{N}$')
ax[axNo].set_ylabel('E-Field\n[mV/m]', fontsize=Label_FontSize, color='black', labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='both', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=0, length=0, width=0)
ax[axNo].set_ylim(E_limits[0], E_limits[1])
ax[axNo].yaxis.set_ticks(np.arange(-75,300,75))
leg = ax[axNo].legend(loc='upper right',fontsize=Legend_FontSize)
for line in leg.get_lines():
    line.set_linewidth(4)


# --- LP Low ---
axNo +=1
filtered= stl.butterFilter().butter_filter(data= data_dict_LP_low['ni'][0]/LP_scale,
                                                              lowcutoff=0.03,
                                                              highcutoff=0.03,
                                                              fs=100,
                                                              filtertype='lowpass',
                                                              order=4
                                                              )
ax[axNo].plot(data_dict_EISCAT_cal_low['L-Shell'][0], data_dict_EISCAT_cal_low['ne'][0]/(np.power(stl.cm_to_m,3)),color='tab:purple', linewidth=Plot_LineWidth+1, label='EISCAT Background n$_{e}$')
ax[axNo].plot(data_dict_LP_low['L-Shell'][0],filtered/(np.power(stl.cm_to_m,3)), color=colorChoice, linewidth=Plot_LineWidth+1)
# ax[axNo].set_ylabel('[10$^{5}$ cm$^{-3}$]', fontsize=Label_FontSize-2, color='black', labelpad=Label_Padding )
ax[axNo].set_ylabel('LP\n[cm$^{-3}$]', fontsize=Label_FontSize-2, color='black', labelpad=Label_Padding )
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize-3, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize - 6, length=Tick_Length-2, width=Tick_Width)
ax[axNo].tick_params(axis='x', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length + 4, width=Tick_Width, pad=Tick_Padding)
ax[axNo].tick_params(axis='x', which='minor', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width, pad=Tick_Padding)
ax[axNo].set_ylim(LP_limit[0], LP_limit[1])
ax[axNo].set_xlabel('L-Shell \n Alt [km]', fontsize=Tick_FontSize-2, weight='bold')
ax[axNo].xaxis.set_label_coords(-0.085, -0.26)
ax[axNo].set_yscale('log')
ax[axNo].minorticks_on()
ax[axNo].legend()

# # --- get L-Shell labels and Alt Labels together ---
xTickLabels = ax[axNo].axes.get_xticklabels()
xTick_Locations = [float(tickVal.get_text()) for tickVal in xTickLabels]
xTick_newLabels_high = [f'{LshellVal}\n{round(data_dict_Lshell_high["Alt"][0][np.abs(data_dict_Lshell_high["L-Shell"][0] - LshellVal).argmin()]/stl.m_to_km)}' for LshellVal in xTick_Locations]
xTick_newLabels_low = [f'{LshellVal}\n{round(data_dict_Lshell_low["Alt"][0][np.abs(data_dict_Lshell_low["L-Shell"][0] - LshellVal).argmin()]/stl.m_to_km)}' for LshellVal in xTick_Locations]
ax[2].set_xticks(xTick_Locations,labels=xTick_newLabels_high)
ax[7].set_xticks(xTick_Locations,labels=xTick_newLabels_low)

# -- Do some minor adjustments to labels/margins/limits ---
for i in range(8):
    ax[i].margins(0)
    ax[i].set_xlim(simLShell_min, simLShell_max)
fig.align_ylabels(ax[:])

for i in [1,2,5,6,7]:
    ax[i].grid(alpha=0.9, which='both')

# for i in [0,1,4,5,6]:
#     ax[i].get_xaxis().set_visible(False)

# --- cbar 1---
cax = fig.add_axes([0.885, 0.796, 0.03, 0.184])
cbar = plt.colorbar(cmap, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbar_TickLabelSize + 5)
cbar.set_label(r'[eV-cm$^{-2}$-s$^{-1}$]', fontsize=cbar_FontSize)

# --- cbar 2---
cax = fig.add_axes([0.885, 0.367, 0.03, 0.184])
cbar = plt.colorbar(cmap, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbar_TickLabelSize + 5)
cbar.set_label(r'[eV-cm$^{-2}$-s$^{-1}$]', fontsize=cbar_FontSize)

fig.subplots_adjust(left=0.13, bottom=0.06, right=0.88, top=0.98,hspace=0)  # remove the space between plots
plt.savefig(rf'{path_to_plots}\Plot1\Plot1_rkt_data_stackplot_base.png', dpi=dpi)
stl.Done(start_time)






