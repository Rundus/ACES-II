# --- PlasmaDensity_vs_Alt.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Using the derived plasma density from the Fixed Langmuir Probe, we can
# see how this density looks when plotted vs altitude



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
from ACESII_code.myImports import *

start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = []

modifier = ''
inputPath_modifier_LP = 'L3\Langmuir' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_EISCAT = r'science\EISCAT\tromso\UHF'
wLP_File = 0
wEISCAT_file = 0
outputPath_modifier = '' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
outputData = False

EISCAT_targetTimes = [dt.datetime(2022,11,20,17,22,00,000000),dt.datetime(2022,11,20,17,27,00,000000)]

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.interpolate import CubicSpline


def PlasmaDensity_vs_Alt(wRocket, wFile, rocketFolderPath, justPrintFileNames):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wRocket-4]
    globalAttrsMod = rocketAttrs.globalAttributes[wRocket-4]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'

    inputFiles_LP = glob(f'{rocketFolderPath}{inputPath_modifier_LP}\{fliers[wRocket-4]}{modifier}\*langmuir_fixed_LowPass*')
    inputFiles_EISCAT = r'C:\Data\ACESII\science\EISCAT\tromso\UHF\MAD6400_2022-11-20_beata_ant@uhfa.cdf'

    input_names_LP = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier_LP}\{fliers[wRocket-4]}{modifier}\\', '') for ifile in inputFiles_LP]
    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier_LP.lower() +'_', '').replace('_v00', '') for ifile in input_names_LP]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles_LP):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return

    print('\n')
    print(color.UNDERLINE + f'Plotting Density vs Altutde information for ACESII_{rocketID}' + color.END)

    # --- get LP file ---
    data_dict_LP = loadDictFromFile(inputFiles_LP[wLP_File])
    Epoch = np.array(data_dict_LP['Epoch'][0])
    ni = np.array(data_dict_LP['ni'][0])
    Alt = np.array(data_dict_LP['Alt'][0])
    Done(start_time)

    # --- get the EISCAT file ---
    data_dict_EISCAT = loadDictFromFile(inputFiles_EISCAT)
    Done(start_time)

    # --- Split the data into upleg and downleg -LP DATA ---
    splitIndex = np.abs(Alt - Alt.max()).argmin()

    newAlt_upleg = Alt[:splitIndex]
    newNi_upleg = ni[:splitIndex]

    newAlt_downleg = Alt[splitIndex:]
    newNi_downleg = ni[splitIndex:]

    # --- Average the EISCAT data over a series of times and altitude ---
    lowIdx, highIdx = np.abs(data_dict_EISCAT['Epoch'][0]-EISCAT_targetTimes[0]).argmin(),np.abs(data_dict_EISCAT['Epoch'][0]-EISCAT_targetTimes[1]).argmin()

    alt_EISCAT = data_dict_EISCAT['range'][0]
    ne_EISCAT_raw = data_dict_EISCAT['ne'][0][lowIdx:highIdx].T
    ne_EISCAT_raw[np.where(np.isnan(ne_EISCAT_raw))] = 0
    ne_EISCAT_avg = np.flip((1E-6)*np.array([np.average(arr[np.nonzero(arr)]) for arr in ne_EISCAT_raw]))


    # --- Density Model ---
    xDataModel = np.linspace(70,10000, 5000)
    R_REF = 6378
    h = 0.06 * (R_REF )  # in km from E's surface
    n0 = 6E4
    n1 = 1.34E7
    z0 = 0.05 * (R_REF )  # in km from E's surface
    n_density = np.array([(n0 * np.exp(-1 * (alt  - z0) / h) + n1 * ((alt ) ** (-1.55))) for alt in xDataModel])  # calculated density (in m^-3)



    #################################
    # --- Plot the Density vs Alt ---
    #################################

    xData = [newNi_upleg, newNi_downleg]
    yData = [newAlt_upleg, newAlt_downleg]

    labels= ['Upleg', 'Downleg']
    fig, ax = plt.subplots(1)
    fig.suptitle(f'ACESII {rocketID} {labels[0]}')
    ax.scatter(xData[0], yData[0], color='red')
    ax.scatter(ne_EISCAT_avg, alt_EISCAT, color='blue')
    ax.plot(n_density, xDataModel,color='green')
    ax.set_xlabel('Density [cm$^{-3}$]')
    ax.set_ylabel('Altitude [km]')
    ax.set_ylim(50, 3000)
    ax.set_yscale('log')
    ax.set_xlim(1E1, 1E7)
    ax.set_xscale('log')
    ax.grid(True)
    plt.show()


    #    # for i in range(2):
    #     fig.suptitle(f'ACESII {rocketID} {labels[i]}')
    #     ax[i].scatter(xData[i], yData[i],color='red')
    #     ax[i].scatter(ne_EISCAT_avg, alt_EISCAT,color='blue')
    #     ax[i].set_xlabel('Density [cm$^{-3}$]')
    #     ax[i].set_ylabel('Altitude [km]')
    #     ax[i].set_ylim(70, 800) if wRocket ==4 else ax[i].set_ylim(70, 200)
    #     ax[i].set_xlim(1E1, 1E7)
    #     ax[i].set_xscale('log')
    #     ax[i].grid(True)
    # plt.show()






    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---








# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder


if len(glob(f'{rocketFolderPath}{inputPath_modifier_LP}\{fliers[wRocket-4]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    PlasmaDensity_vs_Alt(wRocket, 0, rocketFolderPath, justPrintFileNames)