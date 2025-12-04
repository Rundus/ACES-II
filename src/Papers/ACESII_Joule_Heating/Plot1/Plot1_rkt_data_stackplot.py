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
LP_scale = 1E6
E_limits = [-110, 110]
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
data_dicts_BField = [stl.loadDictFromFile(glob(rf'C:/Data/ACESII/L3/B_median/{ACESII.fliers[i]}/*.cdf*')[0]) for i in range(2)]

# E-Field Data
data_dict_Efield_low = stl.loadDictFromFile('C:/Data/ACESII/L2/low/ACESII_36364_l2_EFI_auroral_fullCal.cdf')

# EEPAA Particle Data
data_dicts_flux = [stl.loadDictFromFile(glob(f'C:/Data/ACESII/L3/Energy_Flux/{ACESII.fliers[i]}/*.cdf*')[0]) for i in range(2)]

# Langmuir Probe Density Data
data_dicts_LP = [stl.loadDictFromFile(glob(f'C:/Data/ACESII/L3/Langmuir/{ACESII.fliers[i]}/*fixed_fullCal.cdf*')[0]) for i in range(2)]

data_dict_Lshells = [stl.loadDictFromFile(glob(f'C:/Data/ACESII/coordinates/Lshell/{ACESII.fliers[i]}/*.cdf*')[0]) for i in range(2)]
stl.Done(start_time)


##########################
# --- --- --- --- --- ---
# --- PREPARE THE DATA ---
# --- --- --- --- --- ---
##########################
stl.prgMsg('Down-sampling Data')
# [1] --- Downsample data to only the specified L-Shell range ---

# # Get the simulation range from the spatial environment data
data_dict_sim_spatial = stl.loadDictFromFile('C:/Data/physicsModels/ionosphere/spatial_environment/spatial_environment.cdf')
simLShell_min = deepcopy(data_dict_sim_spatial['simLShell'][0][0])
simLShell_max = deepcopy(data_dict_sim_spatial['simLShell'][0][-1])
stl.Done(start_time)

############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

stl.prgMsg('Plotting Data')

fig, ax = plt.subplots(8, height_ratios=[1.5, 1, 0.5, 0.5, 1.5, 1, 1, 0.5],sharex=False)
fig.set_figwidth(Figure_width)
fig.set_figheight(Figure_height)

axNo = 0
for rkt_idx in range(2):

    # --- EEPAA 0 to 180 DEG---
    if rkt_idx ==1:
        axNo += 1
    cmap = ax[axNo].pcolormesh(data_dicts_flux[rkt_idx]['L-Shell'][0],data_dicts_flux[rkt_idx]['Energy'][0],data_dicts_flux[rkt_idx]['varPhi_E_Parallel'][0].T, cmap=my_cmap,vmin=cbarMin,vmax=cbarMax, norm='log')
    ax[axNo].set_ylabel('Parallel Energy Flux\nEnergy [eV]', fontsize=Label_FontSize, labelpad=Label_Padding)
    ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
    ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
    ax[axNo].set_yscale('log')
    ax[axNo].set_ylim(28, 1E4)

    # --- DC B-Field ---
    axNo +=1
    ax[axNo].plot(data_dicts_BField[rkt_idx]['L-Shell'][0],data_dicts_BField[rkt_idx]['B_p_median'][0], color='tab:blue', linewidth=Plot_LineWidth+2, label='$\Delta B_{p}$')
    ax[axNo].plot(data_dicts_BField[rkt_idx]['L-Shell'][0],data_dicts_BField[rkt_idx]['B_T_median'][0], color='tab:red', linewidth=Plot_LineWidth+2, label='$\Delta B_{T}$')
    ax[axNo].plot(data_dicts_BField[rkt_idx]['L-Shell'][0],data_dicts_BField[rkt_idx]['B_N_median'][0], color='tab:green', linewidth=Plot_LineWidth+2, label='$\Delta B_{N}$')
    ax[axNo].set_ylabel('DC B-Field\n[nT]', fontsize=Label_FontSize,labelpad=Label_Padding)
    ax[axNo].tick_params(axis='y', which='both', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
    ax[axNo].set_ylim(B_limits)
    ax[axNo].yaxis.set_ticks(np.arange(B_limits[0],B_limits[1], 1000))
    leg=ax[axNo].legend(fontsize=Legend_FontSize,loc='upper right')
    for line in leg.get_lines():
        line.set_linewidth(4)

    if rkt_idx==1:
        # ---Low Flyer E-Field ---
        axNo += 1
        ax[axNo].plot(data_dict_Efield_low['L-Shell'][0], data_dict_Efield_low['E_p'][0] * Escale,
                      linewidth=Plot_LineWidth, color='tab:blue', label=r'$E_{p}$')
        ax[axNo].plot(data_dict_Efield_low['L-Shell'][0], data_dict_Efield_low['E_T'][0] * Escale,
                      linewidth=Plot_LineWidth, color='tab:red', label=r'$E_{T}$')
        ax[axNo].plot(data_dict_Efield_low['L-Shell'][0], data_dict_Efield_low['E_N'][0] * Escale,
                      linewidth=Plot_LineWidth, color='tab:green', label=r'$E_{N}$')
        ax[axNo].set_ylabel('E-Field\n[mV/m]', fontsize=Label_FontSize, color='black', labelpad=Label_Padding)
        ax[axNo].tick_params(axis='y', which='both', colors='black', labelsize=Tick_FontSize, length=Tick_Length,
                             width=Tick_Width)
        ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=0, length=0, width=0)
        ax[axNo].set_ylim(E_limits[0], E_limits[1])
        ax[axNo].yaxis.set_ticks(np.arange(-100,100, 50))
        leg = ax[axNo].legend(loc='upper right', fontsize=Legend_FontSize)
        for line in leg.get_lines():
            line.set_linewidth(4)


    # --- LP ---
    axNo +=1
    colorChoice = 'black'
    filtered = stl.butterFilter().butter_filter(data= data_dicts_LP[rkt_idx]['ni'][0],
                                                                  lowcutoff=0.03,
                                                                  highcutoff=0.03,
                                                                  fs=100,
                                                                  filtertype='lowpass',
                                                                  order=4
                                                                  )

    ax[axNo].plot(data_dicts_LP[rkt_idx]['L-Shell'][0], filtered/LP_scale,color=colorChoice, linewidth=Plot_LineWidth+1)
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
    if rkt_idx == 0:
        axNo +=1
        ax[axNo].axis('off')

# # --- get L-Shell labels and Alt Labels together ---
import re
xTickLabels = ax[axNo].axes.get_xticklabels()
# xTick_Locations = [float(re.sub(u"/u2212", "-", i.get_text())) for i in xTickLabels] # negative signs in xTickLabel locations cause issues. Fix them
# xTick_Locations = [float(tickVal.get_text()) for tickVal in xTickLabels]
# xTick_newLabels_high = [f'{LshellVal}\n{round(data_dict_Lshells[0]["Alt"][0][np.abs(data_dict_Lshells[0]["L-Shell"][0] - LshellVal).argmin()]/stl.m_to_km)}' for LshellVal in xTick_Locations]
# xTick_newLabels_low = [f'{LshellVal}\n{round(data_dict_Lshells[1]["Alt"][0][np.abs(data_dict_Lshells[1]["L-Shell"][0] - LshellVal).argmin()]/stl.m_to_km)}' for LshellVal in xTick_Locations]
# ax[2].set_xticks(xTick_Locations,labels=xTick_newLabels_high)
# ax[7].set_xticks(xTick_Locations,labels=xTick_newLabels_low)

# -- Do some minor adjustments to labels/margins/limits ---
for i in range(8):
    ax[i].margins(0)
    ax[i].set_xlim(simLShell_min, simLShell_max)
fig.align_ylabels(ax[:])

for i in [1,2,5,6,7]:
    ax[i].grid(alpha=0.9, which='both')

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
plt.savefig(rf'{path_to_plots}/Plot1/Plot1_rkt_data_stackplot_base.png', dpi=dpi)
stl.Done(start_time)






