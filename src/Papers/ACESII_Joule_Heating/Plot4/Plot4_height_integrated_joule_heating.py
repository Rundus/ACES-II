# --- Plot3_simulation_conductivities.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Stack Plot detailing the outputs of the simulation


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
print(stl.color.UNDERLINE + f'Plot3_simulation_conductivity' + stl.color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
plt.rcParams["font.family"] = "Arial"
dpi = 800

# --- Cbar ---
cbar_TickLabelSize = 14
my_cmap = stl.blue_green_white_yellow_red_cmap()


# --- Plot toggles ---
Figure_width = 12 # in inches
Figure_height =12# in inches
Text_FontSize = 20
Label_FontSize = 20
Tick_FontSize = 20
Tick_Length = 5
Tick_Width = 2
Plot_LineWidth = 1.5
Label_Padding = 15
Tick_Padding = 10
Legend_FontSize = 20
cbar_FontSize = 25

# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---
stl.prgMsg('Loading Data')

# Simulation - E-Field Data
data_dict_Efield = stl.loadDictFromFile('C:\Data\physicsModels\ionosphere\electricField\electric_field.cdf')

# Simulation - Potential
data_dict_potential = stl.loadDictFromFile('C:\Data\physicsModels\ionosphere\electrostaticPotential\electrostaticPotential.cdf')

# Simulation - Conductivity
data_dict_conductivity = stl.loadDictFromFile('C:\Data\physicsModels\ionosphere\conductivity\conductivity.cdf')

# Simulation - Spatial Grid
data_dict_spatial = stl.loadDictFromFile('C:\Data\physicsModels\ionosphere\spatial_environment\spatial_environment.cdf')

stl.Done(start_time)

############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

stl.prgMsg('Plotting Data')

fig, ax = plt.subplots(2)
fig.set_figwidth(Figure_width)
fig.set_figheight(Figure_height)

# Define some variables
simAlt_grid = deepcopy(data_dict_spatial['grid_alt'][0])/stl.m_to_km
simLShell_grid = deepcopy(data_dict_spatial['grid_LShell'][0])
simAlt = deepcopy(data_dict_spatial['simAlt'][0])
simLShell = deepcopy(data_dict_spatial['simLShell'][0])

# --- HEIGHT INTEGRATED PEDERSEN CONDUCTIVITY ---
axNo = 0
ax[axNo].set_ylabel('$\Sigma_{P}$\n[S]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
ax[axNo].plot(simLShell, data_dict_conductivity['Sigma_P'][0],linewidth=Plot_LineWidth,color='black', label=r'Simulation')
ax[axNo].plot(simLShell, data_dict_conductivity['Sigma_P_Robinson'][0],linewidth=Plot_LineWidth,color='tab:blue', label=r'Robinson 1987')
ax[axNo].plot(simLShell, data_dict_conductivity['Sigma_P_K10'][0],linewidth=Plot_LineWidth,color='tab:red', label=r'Kaeppler 2023')
ax[axNo].set_ylim(-0.5, 8)
ax[axNo].legend(fontsize=Legend_FontSize)

# --- HEIGHT INTEGRATED HALL CONDUCTIVITY ---
axNo += 1
ax[axNo].set_ylabel('$\Sigma_{H}$\n[S]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
ax[axNo].plot(simLShell, data_dict_conductivity['Sigma_H'][0],linewidth=Plot_LineWidth,color='black', label=r'Simulation')
ax[axNo].plot(simLShell, data_dict_conductivity['Sigma_H_Robinson'][0],linewidth=Plot_LineWidth,color='tab:blue', label=r'Robinson 1987')
ax[axNo].plot(simLShell, data_dict_conductivity['Sigma_H_K10'][0],linewidth=Plot_LineWidth,color='tab:red', label=r'Kaeppler 2023')
ax[axNo].legend(fontsize=Legend_FontSize)
ax[axNo].set_ylim(-0.5, 13)

# -- Do some minor adjustments to labels/margins/limits ---
for i in range(axNo+1):
    ax[i].margins(0)
    ax[i].set_xlim(data_dict_spatial['simLShell'][0][0], data_dict_spatial['simLShell'][0][-1])
    ax[i].tick_params(axis='x', which='major', colors='black', labelsize=Tick_FontSize+10, length=Tick_Length,width=Tick_Width)

    if i <= 0:
        ax[i].set_xticklabels([])

fig.align_ylabels(ax[:])

# fig.subplots_adjust(left=0.12, bottom=0.06, right=0.88, top=0.98,hspace=0.08)  # remove the space between plots
plt.tight_layout()
plt.savefig(r'C:\Users\cfelt\Desktop\Research\ACESII\Feltman2025_ACESII_JouleHeating\PLOTS\Plot4\Plot4_simulation_HI_conductivity_base.png', dpi=dpi)
stl.Done(start_time)






