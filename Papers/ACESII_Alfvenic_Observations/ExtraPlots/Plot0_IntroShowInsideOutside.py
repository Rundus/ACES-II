# --- Plot0_IntroShowInsideOutside.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Isolate Alfvenic Dispersions and perform cross-correlation analysis
# to determine values that are fitted linearly to get the height of the source region

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---
modifier = ''
inputPath_modifier = 'l2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder

justPrintFileNames = False
wFile = 0


# --- --- --- --- --- -
# --- PLOT TOGGLES ---
# --- --- --- --- ----
targetTimes = [[dt.datetime(2022,11,20,17,24,55,000000),dt.datetime(2022,11,20,17,25,00,000000)], # isolated
                [dt.datetime(2022,11,20,17,25,22,500000),dt.datetime(2022,11,20,17,25,27,500000)], # edge type
               [dt.datetime(2022,11,20,17,24,16,000000),dt.datetime(2022,11,20,17,24,21,000000)], # underneath/auroral
               [dt.datetime(2022,11,20,17,25,3,000000),dt.datetime(2022,11,20,17,25,8,500000)] # unclear mix
               ]
wPitch = 2
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
mycmap = apl_rainbow_black0_cmap()
cbar_low, cbar_high = 5E6, 5E9

# plotparams
figure_height = 8
figure_width = 15
labelPadding = -0.25
textFontSize = 10
titleFontSize = 18
labelsFontSize = 22
SidelabelAdjust = 2.6
tickFontSize = 16
cbartickFontSize = 20
tickWidth = 3
tickLength = 8
cbarFont = 15
dpi = 200


# ---------------------------


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.Science.AlfvenSingatureAnalysis.Particles.dispersionAttributes import dispersionAttributes




def Plot0_IntroShowInsideOutside(wRocket, rocketFolderPath):

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}{modifier}\*eepaa_fullcal*')
    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}{modifier}\\', '') for ifile in inputFiles]
    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return

    # --- get the data from the ESA file ---
    data_dict = loadDictFromFile(inputFiles[0])

    Energy = data_dict['Energy'][0]
    Pitch = data_dict['Pitch_Angle'][0]
    Epoch = data_dict['Epoch'][0]
    diffEFlux = data_dict['Differential_Energy_Flux'][0][:,wPitch,:]

    #################################################
    # --- PLOT THE DISPERSION IN ORDER TO ISOLATE ---
    #################################################

    fig, ax = plt.subplots(nrows=2,ncols=2)
    fig.set_figwidth(figure_width)
    fig.set_figheight(figure_height)


    for pltIdx,ttime in enumerate(targetTimes):

        # --- GET TH DATA ---
        lowIdx, highIdx = np.abs(Epoch - targetTimes[pltIdx][0]).argmin(),np.abs(Epoch - targetTimes[pltIdx][1]).argmin()
        Epoch_dis = Epoch[lowIdx:highIdx]
        fluxData = diffEFlux[lowIdx:highIdx].T


        # --- PLOT ---
        colIdx = 0 if pltIdx in [0,2] else 1
        rowIdx = 0 if pltIdx in [0,1] else 1

        cmap = ax[rowIdx][colIdx].pcolormesh(Epoch_dis, Energy, fluxData, cmap=mycmap, vmin=cbar_low, vmax=cbar_high,norm='log')
        ax[rowIdx][colIdx].set_yscale('log')
        ax[rowIdx][colIdx].set_ylim(28, 1E4)

        # LABELS AND TICKS

        # ticks
        ax[rowIdx][colIdx].tick_params(axis='y', which='major', labelsize=tickFontSize, width=tickWidth-1,length=tickLength)
        ax[rowIdx][colIdx].tick_params(axis='y', which='minor', labelsize=tickFontSize, width=tickWidth-1,length=tickLength / 2)
        ax[rowIdx][colIdx].tick_params(axis='x', which='major', labelsize=tickFontSize, width=tickWidth,length=tickLength)
        ax[rowIdx][colIdx].tick_params(axis='x', which='minor', labelsize=tickFontSize, width=tickWidth,length=tickLength / 2)

        # labels
        if colIdx in [0]:
            ax[rowIdx][colIdx].set_ylabel('Energy [eV]',fontsize=labelsFontSize, weight='bold')
        if rowIdx in [1]:
            ax[rowIdx][colIdx].set_xlabel('Time [UTC]', fontsize=labelsFontSize,weight='bold')
        if colIdx in [1]:
            ax[rowIdx][colIdx].set_yticklabels([])



    # --- ADD THE COLORBAR ---
    cax = fig.add_axes([0.94, 0.09, 0.02, 0.88])
    cbar = plt.colorbar(cmap, cax=cax)
    cbar.ax.minorticks_on()
    cbar.ax.tick_params(labelsize=cbartickFontSize)
    fig.subplots_adjust(left=0.07, bottom=0.09, right=0.93, top=0.97, wspace=0.05, hspace=0.14)  # remove the space between plots
    plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot0\Plot_0.png',dpi=dpi)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
Plot0_IntroShowInsideOutside(wRocket=4,rocketFolderPath = ACES_data_folder)