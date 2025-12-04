# --- Plot3_simulation_conductivities.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Stack Plot detailing the outputs of the simulation


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

from src.ACESII.my_imports import *

start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import matplotlib.pyplot as plt

print(stl.color.UNDERLINE + f'Plot3_simulation_conductivity' + stl.color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
plt.rcParams["font.family"] = "Arial"
dpi = 800
altScale = stl.m_to_km

# --- Cbar ---
cbar_TickLabelSize = 14
my_cmap = stl.blue_green_white_yellow_red_cmap()


# --- Plot toggles ---
Figure_width = 15 # in inches
Figure_height =17# in inches
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

# Simulation - Currents
data_dict_currents = stl.loadDictFromFile('C:\Data\physicsModels\ionosphere\currents\currents_01.cdf')

# Simulation - Joule Heating
data_dict_joule = stl.loadDictFromFile('C:\Data\physicsModels\ionosphere\joule_heating\joule_heating.cdf')

stl.Done(start_time)


############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

stl.prgMsg('Plotting Data')

fig, ax = plt.subplots(6, height_ratios=[1, 1, 1, 1, 1, 1], sharex=False)
fig.set_figwidth(Figure_width)
fig.set_figheight(Figure_height)

# Define some variables
simAlt_grid = deepcopy(data_dict_spatial['grid_alt'][0])/stl.m_to_km
simLShell_grid = deepcopy(data_dict_spatial['grid_LShell'][0])
simAlt = deepcopy(data_dict_spatial['simAlt'][0])
simLShell = deepcopy(data_dict_spatial['simLShell'][0])

# --- HEIGHT INTEGRATED PEDERSEN CONDUCTIVITY ---
axNo = 0
ax[axNo].set_ylabel('[S]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
ax[axNo].plot(simLShell, data_dict_conductivity['Sigma_P'][0],linewidth=Plot_LineWidth,color='black', label=r'$\Sigma_{P}$ Simulation')
ax[axNo].plot(simLShell, data_dict_conductivity['Sigma_P_Robinson'][0],linewidth=Plot_LineWidth,color='tab:blue', label=r'$\Sigma_{P}$ Robinson 1987')
ax[axNo].plot(simLShell, data_dict_conductivity['Sigma_P_K10'][0],linewidth=Plot_LineWidth,color='tab:red', label=r'$\Sigma_{P}$ Kaeppler 2023')
ax[axNo].set_ylim(-0.5, 8)
ax[axNo].legend(fontsize=Legend_FontSize)

# --- HEIGHT INTEGRATED HALL CONDUCTIVITY ---
axNo += 1
ax[axNo].set_ylabel('$\Sigma_{H}$\n[S]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
ax[axNo].plot(simLShell, data_dict_conductivity['Sigma_H'][0],linewidth=Plot_LineWidth,color='black', label=r'$\Sigma_{P}$ Simulation')
ax[axNo].plot(simLShell, data_dict_conductivity['Sigma_H_Robinson'][0],linewidth=Plot_LineWidth,color='tab:blue', label=r'$\Sigma_{P}$ Robinson 1987')
ax[axNo].plot(simLShell, data_dict_conductivity['Sigma_H_K10'][0],linewidth=Plot_LineWidth,color='tab:red', label=r'$\Sigma_{P}$ Kaeppler 2023')
ax[axNo].legend(fontsize=Legend_FontSize)
ax[axNo].set_ylim(-0.5, 13)

# --- DC PEDERSEN CURRENT (HEIGHT-INTEGRATED) ---
axNo += 1
ax[axNo].set_ylabel('[A/m$^{2}$]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
ax[axNo].plot(simLShell, data_dict_currents['J_P_HI_DC'][0],linewidth=Plot_LineWidth,color='tab:blue', label=r'$J_{P}$ (Height Integrated) DC')
ax[axNo].plot(simLShell, data_dict_currents['J_P_HI_AC'][0],linewidth=Plot_LineWidth,color='tab:green', label=r'$J_{P}$ (Height Integrated) AC')
ax[axNo].legend(fontsize=Legend_FontSize)
ax[axNo].set_ylim(-0.15,0.15)

# --- DC HALL CURRENT (HEIGHT-INTEGRATED) ---
axNo += 1
ax[axNo].set_ylabel('[A/m$^{2}$]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
ax[axNo].plot(simLShell, data_dict_currents['J_H_HI_DC'][0],linewidth=Plot_LineWidth,color='tab:blue', label=r'$J_{H}$ (Height Integrated) DC')
ax[axNo].plot(simLShell, data_dict_currents['J_H_HI_AC'][0],linewidth=Plot_LineWidth,color='tab:green', label=r'$J_{H}$ (Height Integrated) AC')
ax[axNo].legend(fontsize=Legend_FontSize)
ax[axNo].set_ylim(-0.15,0.15)

# --- DC ELECTRIC FIELD  ---
axNo += 1
ax[axNo].set_ylabel('$E_{T}$\nAltitude [km]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
cmap1 = ax[axNo].pcolormesh(simLShell, simAlt/altScale, data_dict_Efield['E_N'][0].T, cmap='bwr', vmin=-0.4, vmax=0.4)

# --- DC JOULE HEATING RATE  ---
axNo += 1
ax[axNo].set_ylabel('$q_{j}$\nAltitude [km]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
ax[axNo].tick_params(axis='x', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width, pad=Tick_Padding)
cmap2 = ax[axNo].pcolormesh(simLShell, simAlt/altScale, data_dict_joule['q_j_DC'][0].T, cmap=stl.apl_rainbow_black0_cmap())
ax[axNo].set_xlabel('L-Shell', fontsize=Label_FontSize, weight='bold')
ax[axNo].xaxis.set_label_coords(-0.085, -0.03)

# --- cbar 1---
cax = fig.add_axes([0.885, 0.2, 0.03, 0.15])
cbar = plt.colorbar(cmap1, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbar_TickLabelSize + 5)
cbar.set_label(r'[V/m]', fontsize=cbar_FontSize)

# --- cbar 2---
cax = fig.add_axes([0.885, 0.05, 0.03, 0.15])
cbar = plt.colorbar(cmap2, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbar_TickLabelSize + 5)
cbar.set_label(r'[W/$m^{3}$]', fontsize=cbar_FontSize)

# -- Do some minor adjustments to labels/margins/limits ---
simLShell_min = deepcopy(data_dict_spatial['simLShell'][0][0])
simLShell_max = deepcopy(data_dict_spatial['simLShell'][0][-1])
for i in range(6):
    ax[i].margins(0)
    ax[i].set_xlim(simLShell_min, simLShell_max)
fig.align_ylabels(ax[:])

for i in [0, 1, 2, 3]:
    ax[i].grid(alpha=0.9, which='both')

fig.subplots_adjust(left=0.13, bottom=0.06,right=0.85, top=0.98,hspace=0)  # remove the space between plots
plt.savefig(r'C:\Users\cfelt\OneDrive - University of Iowa\Research\ACESII\Feltman2025_ACESII_JouleHeating\PLOTS\Plot6\Plot6_simulation_results_base.png', dpi=dpi)
stl.Done(start_time)






