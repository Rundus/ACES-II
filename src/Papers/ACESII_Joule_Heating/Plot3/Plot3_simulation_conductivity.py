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
my_cmap = stl.apl_rainbow_black0_cmap()
my_cmap.set_bad(color=(0,0,0))

# --- Plot toggles ---
Figure_width = 13 # in inches
Figure_height =18# in inches
Text_FontSize = 20
Label_FontSize = 25
Tick_FontSize = 20
Tick_Length = 7
Tick_Width = 3
Plot_LineWidth = 1.5
Label_Padding = 15
Tick_Padding = 10
Legend_FontSize = 20
cbar_FontSize = 25

# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---
stl.prgMsg('Loading Data')

# Simulation - Neutral Environemnt
data_dict_neutral = stl.loadDictFromFile(r'C:\Data\physicsModels\ionosphere\neutral_environment\neutral_environment.cdf')

# Simulation - Conductivity
data_dict_conductivity = stl.loadDictFromFile('C:\Data\physicsModels\ionosphere\conductivity\conductivity.cdf')

# Simulation - plasma environment
data_dict_plasma = stl.loadDictFromFile('C:\Data\physicsModels\ionosphere\plasma_environment\plasma_environment.cdf')

# Simulation - Spatial Grid
data_dict_spatial = stl.loadDictFromFile('C:\Data\physicsModels\ionosphere\spatial_environment\spatial_environment.cdf')

stl.Done(start_time)

############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

stl.prgMsg('Plotting Data')

fig, ax = plt.subplots(6)
fig.set_figwidth(Figure_width)
fig.set_figheight(Figure_height)

# Define some variables
simAlt_grid = deepcopy(data_dict_spatial['grid_alt'][0])/stl.m_to_km
simLShell_grid = deepcopy(data_dict_spatial['grid_LShell'][0])

# --- PEDERSEN CONDUCTIVITY ---
cbarMin, cbarMax = 4E-8, 4E-4
axNo = 0
cmap1 = ax[axNo].pcolormesh(simLShell_grid,simAlt_grid,data_dict_conductivity['sigma_P'][0], cmap=my_cmap,vmin=cbarMin,vmax=cbarMax, norm='log')
ax[axNo].set_ylabel('Pedersen\nAltitude [km]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)


# --- HALL CONDUCTIVITY ---
cbarMin, cbarMax = 4E-8, 4E-4
axNo += 1
cmap2 = ax[axNo].pcolormesh(simLShell_grid,simAlt_grid,data_dict_conductivity['sigma_H'][0], cmap=my_cmap,vmin=cbarMin,vmax=cbarMax, norm='log')
ax[axNo].set_ylabel('Hall\nAltitude [km]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)

# --- PARALLEL CONDUCTIVITY ---
cbarMin, cbarMax = 1E-2, 1E2
axNo += 1
cmap3 = ax[axNo].pcolormesh(simLShell_grid,simAlt_grid,data_dict_conductivity['sigma_D'][0], cmap=my_cmap,vmin=cbarMin,vmax=cbarMax, norm='log')
ax[axNo].set_ylabel('Parallel\nAltitude [km]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)


# --- ELECTRON DENSITY ---
cbarMin, cbarMax = 3E3, 3E6
axNo +=1
cmap4 = ax[axNo].pcolormesh(simLShell_grid,simAlt_grid,data_dict_conductivity['ne_total'][0]/(1E6), cmap=my_cmap,vmin=cbarMin,vmax=cbarMax,norm='log')
ax[axNo].set_ylabel('Electron Density \nAltitude [km]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)


# --- ELECTRON TEMPERATURE ---
axNo +=1
# cbarMin, cbarMax = 0, 0.2
cbarMin, cbarMax = 0, 1800
cmap5 = ax[axNo].pcolormesh(simLShell_grid,simAlt_grid,data_dict_plasma['Te'][0], cmap=my_cmap,vmin=cbarMin,vmax=cbarMax)
ax[axNo].set_ylabel('Electron Temp.\nAltitude [km]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)


# --- ION TEMPERATURE ---
axNo +=1
# cbarMin, cbarMax = 0, 0.2
cbarMin, cbarMax = 0, 1800
# cmap6 = ax[axNo].pcolormesh(simLShell_grid,simAlt_grid,data_dict_plasma['Ti'][0], cmap=my_cmap,vmin=cbarMin,vmax=cbarMax)
cmap6 = ax[axNo].pcolormesh(simLShell_grid,simAlt_grid,data_dict_neutral['Tn'][0], cmap=my_cmap,vmin=cbarMin,vmax=cbarMax)
ax[axNo].set_ylabel('Neutral Temp.\nAltitude [km]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
ax[axNo].set_xlabel('L-Shell', fontsize=Label_FontSize+4, weight='bold')
ax[axNo].xaxis.set_label_coords(-0.085, -0.06)

# -- Do some minor adjustments to labels/margins/limits ---
for i in range(axNo+1):
    ax[i].margins(0)
    ax[i].set_xlim(data_dict_spatial['simLShell'][0][0], data_dict_spatial['simLShell'][0][-1])
    ax[i].tick_params(axis='x', which='major', colors='black', labelsize=Tick_FontSize+10, length=Tick_Length,width=Tick_Width)

    if i <= 4:
        ax[i].set_xticklabels([])

fig.align_ylabels(ax[:])

###################
# --- COLORBARS ---
###################
# --- cbar 1: Pedersen ---
cax = fig.add_axes([0.885, 0.836, 0.03, 0.144])
cbar = plt.colorbar(cmap1, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbar_TickLabelSize + 5)
cbar.set_label(r'[S/m]', fontsize=cbar_FontSize)

# --- cbar 2: Hall---
cax = fig.add_axes([0.885, 0.681, 0.03, 0.144 ])
cbar = plt.colorbar(cmap2, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbar_TickLabelSize + 5)
cbar.set_label(r'[S/m]', fontsize=cbar_FontSize)

# --- cbar 3: Parallel---
cax = fig.add_axes([0.885, 0.525, 0.03, 0.144])
cbar = plt.colorbar(cmap3, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbar_TickLabelSize + 5)
cbar.set_label(r'[S/m]', fontsize=cbar_FontSize)

# --- cbar 4: Electron density ---
cax = fig.add_axes([0.885, 0.37, 0.03, 0.144])
cbar = plt.colorbar(cmap4, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbar_TickLabelSize + 5)
cbar.set_label(r'[cm$^{-3}$]', fontsize=cbar_FontSize)

# --- cbar 5: Electron Temperature ---
cax = fig.add_axes([0.885, 0.215, 0.03, 0.144])
cbar = plt.colorbar(cmap5, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbar_TickLabelSize + 5)
cbar.set_label(r'[K]', fontsize=cbar_FontSize)

# --- cbar 6: Ion Temperature ---
cax = fig.add_axes([0.885, 0.06, 0.03, 0.144])
cbar = plt.colorbar(cmap6, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbar_TickLabelSize + 5)
cbar.set_label(r'[K]', fontsize=cbar_FontSize)

fig.subplots_adjust(left=0.12, bottom=0.06, right=0.88, top=0.98,hspace=0.08)  # remove the space between plots
plt.savefig(r'C:\Users\cfelt\Desktop\Research\ACESII\Feltman2025_ACESII_JouleHeating\PLOTS\Plot3\Plot3_simulation_conductivity_base.png', dpi=dpi)
stl.Done(start_time)






