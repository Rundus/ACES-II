# --- Plots6_pitchAnglePlots.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Recreation of DiffEFlux pitch angle Plots, focusing on
# a few particle signatures

# TODO: sort the white line on the spectrogram plots


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

from myImports import *
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import spaceToolsLib as stl
from Science.AlfvenSingatureAnalysis.Particles.dispersionAttributes import dispersionAttributes

print(color.UNDERLINE + f'Plot7p5_DetermineSourceDistributionCurves' + color.END)

#################
# --- TOGGLES ---
#################
wPitchsToFit = [2, 3, 4, 5]

######################
# --- PLOT TOGGLES ---
######################
figure_height = 6*(2.5)
figure_width = (16)

dpi = 200

Label_FontSize = 24
Label_Padding = 0.5
velSpace_Label_FontSize = 25

Text_FontSize = 18
Title_FontSize = 18
Tick_FontSize = 25
Tick_Width = 3
Tick_Length = 7
Plot_Linewidth = 6.5
vertical_lineStyles = ['dashdotdotted','dashdot','dashdotdotted','dashdot']
cbar_Fontsize = 24

my_cmap = stl.apl_rainbow_black0_cmap()
my_cmap.set_bad(color=(0,0,0))
my_cmap.set_under('black')
cbarMin,cbarMax = 1E6, 1E9


# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---

prgMsg('Loading Data')

# Dispersive region fit data
inputFile = rf'C:\Data\ACESII\science\invertedV\TempDensityPotential_Fitting\DispersiveRegion\invertedVFitdata_DispersiveRegion.cdf'
data_dict_dispersiveFitData = loadDictFromFile(inputFilePath=inputFile)

# primary inverted-V data
inputFile = rf'C:\Data\ACESII\science\invertedV\TempDensityPotential_Fitting\PrimaryInvertedV\invertedVFitdata_PrimaryV.cdf'
data_dict_invertedVFitData = loadDictFromFile(inputFilePath=inputFile)

data_dicts = [data_dict_dispersiveFitData, data_dict_invertedVFitData]

# DiffNFlux data
inputFile = 'C:\Data\ACESII\L2\high\ACESII_36359_l2_eepaa_fullCal.cdf'
data_dict_diffNFlux = loadDictFromFile(inputFilePath=inputFile)

# DistributionFunc Data
inputFile = 'C:\Data\ACESII\L3\DistFunc\high\ACESII_36359_distFunc_eepaa_fullCal.cdf'
data_dict_dist = loadDictFromFile(inputFilePath=inputFile)


# Attitude Data
inputFile = glob(r'C:\Data\ACESII\attitude\low\*.cdf')[0]
data_dict_attitude = loadDictFromFile(inputFilePath=inputFile)
Done(start_time)


############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

prgMsg('Beginning Plot')
fig = plt.figure()
fig.set_size_inches(figure_width, figure_height)

gs0 = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[5/6,1/6])
gs00 = gs0[0].subgridspec(5, 2)
gs01 = gs0[1].subgridspec(1, 1)




for tmeRangeIdx in range(2):

    # get the data dict set
    data_dict = data_dicts[tmeRangeIdx]

    # pull out the individual data
    Epoch = data_dict['Epoch'][0]
    Density = data_dict['Density'][0]
    Temperature = data_dict['Temperature'][0]
    Potential = data_dict['potential'][0]
    ChiSquare = data_dict['chiSquare'][0]
    peakE = data_dict['peakE'][0]
    ptchVal = data_dict['pitchAngle'][0]

    # --- Do some initial processing on the data ---

    ### find the average fitted peakE for each time##

    ptchChangeOverIndicies = [0]+[idx+1 for idx in range(len(ptchVal)-1) if 1E2 >np.abs(ptchVal[idx+1] - ptchVal[idx]) > 1] + [len(ptchVal)]
    peakEarrays = np.array([np.array( peakE[ptchChangeOverIndicies[i]:ptchChangeOverIndicies[i+1]] ) for i in range(len(ptchChangeOverIndicies)-1)]).T
    avgPeakE = np.sum(peakEarrays,axis=1)/len(wPitchsToFit)
    avgPeakE_Epoch = Epoch[0:ptchChangeOverIndicies[1]]

    # --- [1] make the spectrogram from 0 to 40 deg ---
    axSpec = fig.add_subplot(gs00[0, tmeRangeIdx])
    diffFlux_avg = np.array(data_dict_diffNFlux['Differential_Energy_Flux'][0])
    diffFlux_avg[diffFlux_avg < 0] = 0
    diffFlux_avg = diffFlux_avg[:, wPitchsToFit, :]
    diffFlux_avg = np.sum(diffFlux_avg, axis=1) / (len(wPitchsToFit))

    # reduce the spectrogram to only the fit region
    lowIdx, highIdx = np.abs(data_dict_diffNFlux['Epoch'][0] - Epoch[0]).argmin(), np.abs(data_dict_diffNFlux['Epoch'][0] - Epoch[-1]).argmin()

    cmap = axSpec.pcolormesh(data_dict_diffNFlux['Epoch'][0][lowIdx:highIdx+1], data_dict_diffNFlux['Energy'][0], diffFlux_avg[lowIdx:highIdx+1].T, cmap=my_cmap, vmin=cbarMin, vmax=cbarMax, norm='log') # plot spectrogram
    axSpec.scatter(avgPeakE_Epoch, avgPeakE, marker='o',color='white')
    axSpec.set_ylabel(rf'Energy [eV]', fontsize=Label_FontSize,labelpad=Label_Padding) if tmeRangeIdx == 0 else None
    axSpec.set_yscale('log')
    axSpec.set_ylim(28.22, 2500)

    # --- [2] Make the Chi Square Plot ---
    axChi = fig.add_subplot(gs00[1, tmeRangeIdx])
    axChi.scatter(Epoch,ChiSquare, marker='o',color='black')
    axChi.set_ylim(0.1, 100)
    axChi.set_ylabel(r'$\chi ^{2}_{\nu}$', fontsize=Label_FontSize,labelpad=Label_Padding) if tmeRangeIdx == 0 else None
    axChi.set_yscale('log')

    # --- [3] Make the Density Plot ---
    axDensity = fig.add_subplot(gs00[2, tmeRangeIdx])
    axDensity.scatter(Epoch, Density, marker='o',color='black')
    axDensity.set_ylim(0, 7)
    axDensity.set_ylabel('$n_{0}$  [cm$^{-3}$]', fontsize=Label_FontSize,labelpad=Label_Padding) if tmeRangeIdx == 0 else None

    # --- [4] Make the Parallel Potential Plot ---
    axV0 = fig.add_subplot(gs00[3, tmeRangeIdx])
    axV0.scatter(Epoch, Potential, marker='o',color='black')
    axV0.set_ylim(10, 500)
    axV0.set_ylabel('V [eV]', fontsize=Label_FontSize,labelpad=Label_Padding) if tmeRangeIdx == 0 else None
    axV0.set_yscale('log')

    # --- [4] Make the Temperature Plot ---
    axTemp = fig.add_subplot(gs00[4, tmeRangeIdx])
    axTemp.scatter(Epoch, Temperature, marker='o',color='black')
    axTemp.set_ylim(10, 1000)
    axTemp.set_ylabel('$T_{e}$ [eV]', fontsize=Label_FontSize,labelpad=Label_Padding) if tmeRangeIdx == 0 else None
    axTemp.set_yscale('log')

    # Make general adjustments
    axes = [axSpec,axChi, axDensity, axV0, axTemp]
    fig.align_ylabels(axes[:])
    for i,ax in enumerate(axes):
        ax.minorticks_on()
        ax.tick_params(axis='y', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
        ax.tick_params(axis='x', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
        if i > 0:
            ax.grid(True,which='Major', alpha=0.5)

        if i < 4:
            ax.tick_params(axis='x', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width,labelbottom=False,bottom=False,top=False)

cax = fig.add_axes([0.91, 0.06, 0.02, 0.103])
cbar = plt.colorbar(cmap, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbar_Fontsize)

fig.subplots_adjust(left=0.05, bottom=0.06, right=0.89, top=0.985, wspace=None,hspace=0.05)  # remove the space between plots
plt.savefig(rf'C:\Users\cfelt\Desktop\Research\Feltman2024_ACESII_Alfven_Observations\PLOTS\Plot7p5\Plot7p5_base.png', dpi=dpi)
Done(start_time)