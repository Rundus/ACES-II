# --- SignatureIdentificationPlot.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: In order to talk about the dispersion signatures, we need to create a
# helpful look plot



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

# --- --- --- --- ---

import time

import matplotlib.pyplot as plt

from ACESII_code.class_var_func import Done, setupPYCDF

start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False
wPitch = 2
EnergyBounds = [16, 42]

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

modifier = ''
inputPath_modifier = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_trajectory = 'trajectories'
outputPath_modifier = 'science\AlfvenSignatureAnalysis' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import datetime as dt
from ACESII_code.missionAttributes import ACES_mission_dicts
from ACESII_code.data_paths import ACES_data_folder, fliers
from ACESII_code.class_var_func import color, prgMsg, L2_TRICE_Quick,loadDictFromFile
from glob import glob
from os.path import getsize

setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)

def IdentificationPlot(wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf')
    inputFiles_trajectory = glob(f'{rocketFolderPath}{inputPath_modifier_trajectory}\{fliers[wflyer]}{modifier}\*_ILat_ILong*')
    outputFiles = glob(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    input_names_trajectory = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles_trajectory]
    output_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    output_names_searchable = [ofile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(outputPath_modifier.lower() +'_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '')
    dataFile_name_trajectory = inputFiles_trajectory[0].replace(f'{rocketFolderPath}{inputPath_modifier_trajectory}\{fliers[wflyer]}{modifier}\\', '')

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in dataFile_name:
            descriptiorNam = ['EEPAA', 'LEESA', 'IEPAA', 'Langmuir_Probe']
            wInstr = [index, instr, descriptiorNam[index]]

    fileoutName = f'ACESII_{rocketID}_{wInstr[1]}_InvertedV_removed.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the ESA file ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict = loadDictFromFile(inputFiles[wFile],{})
        data_dict['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_esa'][0][i]) for i in (range(len(data_dict['Epoch_esa'][0])))])
        Done(start_time)

        # --- get data from Trajectory Files ---
        prgMsg(f'Loading data from {inputPath_modifier_trajectory} Files')
        data_dict_traj = loadDictFromFile(inputFiles_trajectory[0], {})
        data_dict_traj['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_traj['Epoch_esa'][0][i]) for i in (range(len(data_dict_traj['Epoch_esa'][0])))])
        Done(start_time)

        ########################################
        # --- CREATE THE IDENTIFICATION PLOT ---
        ########################################
        prgMsg('creating identification plot')

        # --- prepare the data ---
        from ACESII_code.Science.AlfvenSingatureAnalysis.Particles.dispersionAttributes import dispersionAttributes
        Energy = data_dict['Energy'][0][EnergyBounds[0]:EnergyBounds[1]]
        Pitch = data_dict['Pitch_Angle'][0]
        targetTimes = [pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 24, 55, 000000)),
                       pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 25, 9, 000000))]
        limitIndexes = [np.abs(data_dict['Epoch_esa'][0] - targetTimes[0]).argmin(), np.abs(data_dict['Epoch_esa'][0] - targetTimes[1]).argmin()]
        Epoch_tt2000 = np.array([data_dict['Epoch_esa'][0][i] for i in range(limitIndexes[0], limitIndexes[1])])
        Epoch = np.array([pycdf.lib.tt2000_to_datetime(Epoch_tt2000[i]) for i in range(len(Epoch_tt2000))])
        esaData_allpitches = np.array(data_dict[wInstr[1]][0][limitIndexes[0]:limitIndexes[1]])
        esaData = np.array([esaData_allpitches[tme][wPitch][EnergyBounds[0]:EnergyBounds[1]] for tme in range(len(esaData_allpitches))], dtype='int')

        # --- get the dispersion information ---

        # delta Times
        deltaT = [(pycdf.lib.datetime_to_tt2000(time[1]) - pycdf.lib.datetime_to_tt2000(time[0]))/1E9 for time in dispersionAttributes.keyDispersionTimes]

        # delta Energy
        deltaE = [ data_dict['Energy'][0][engylimit] - data_dict['Energy'][0][-1]  for engylimit in dispersionAttributes.keyDispersionEnergyLimits]

        # delta km
        deltakm = []
        for time in dispersionAttributes.keyDispersionTimes: # need to find the closest index to the dispersion times
            indexLower = np.abs(data_dict_traj['Epoch_esa'][0] - pycdf.lib.datetime_to_tt2000(time[0])).argmin()
            indexUpper = np.abs(data_dict_traj['Epoch_esa'][0] - pycdf.lib.datetime_to_tt2000(time[1])).argmin()
            deltakm.append(data_dict_traj['distanceFromLaunchPoint_km'][0][indexUpper] - data_dict_traj['distanceFromLaunchPoint_km'][0][indexLower])

        # dispersion names
        dispNames = [f's{i+1}' for i in range(len(deltaE))]
        dispersionInfo = [[dispNames[i], deltaT[i], deltaE[i], deltakm[i]] for i in range(len(dispNames))]

        # find the xticks to make it look nice
        ticks = [
            pycdf.lib.tt2000_to_datetime(targetTimes[0]),
            dt.datetime(2022, 11, 20, 17, 24, 56, 00),
            dt.datetime(2022, 11, 20, 17, 24, 57, 00),
            dt.datetime(2022, 11, 20, 17, 24, 58, 00),
            dt.datetime(2022, 11, 20, 17, 24, 59, 00),
            dt.datetime(2022, 11, 20, 17, 25, 00, 00),
            dt.datetime(2022, 11, 20, 17, 25, 1, 00),
            dt.datetime(2022, 11, 20, 17, 25, 2, 00),
            dt.datetime(2022, 11, 20, 17, 25, 3, 00),
            dt.datetime(2022, 11, 20, 17, 25, 4, 00),
            dt.datetime(2022, 11, 20, 17, 25, 5, 00),
            dt.datetime(2022, 11, 20, 17, 25, 6, 00),
            dt.datetime(2022, 11, 20, 17, 25, 7, 00),
            dt.datetime(2022, 11, 20, 17, 25, 8, 00),
            pycdf.lib.tt2000_to_datetime(targetTimes[1])
        ]

        # --- make the plot ---
        fig, ax = plt.subplots(dpi=100)
        fig.set_figwidth(18)
        fig.set_figheight(10)
        cmap = ax.pcolormesh(Epoch, Energy, esaData.T, vmin=0, vmax=40, cmap='turbo')
        cbar = plt.colorbar(cmap)
        cbar.set_label(label='counts', weight='bold', size=12)
        ax.set_xticks(ticks)
        ax.set_xticklabels([tme.strftime("%M:%S") for tme in ticks], size=15)
        ax.set_xlabel('Epoch', size=15)
        ax.set_ylabel('Energy [eV]', size=15)
        fig.suptitle(f'ACESII 36359 \n {wInstr[1].upper()}\n' + rf' Pitch {Pitch[wPitch]}$^\circ$')

        xPosAdjust = 0.5

        # plot the dispersion information
        EnergyHeights = np.array([
        24,  # s1
        27,  # s2
        26,  # s3
        20,  # s4
        19,  # s5
        18,  # s6
        25,  # s7
        24,  # s8
        23,  # s9
        22,  # s10
        24,  # s11
        22,  # s12
        27,  # s13
        31  # s14
        ])

        for i, info in enumerate(dispersionInfo):
            textlabel= f'{info[0]} \n' \
                       f'$\Delta t$: {info[1]}s\n' \
                       f'$\Delta E$: {round(info[2])}eV\n' \
                       f'$\Delta x$: {round(info[3],2)}km\n'
            xlocation = pycdf.lib.tt2000_to_datetime(int(pycdf.lib.datetime_to_tt2000(dispersionAttributes.keyDispersionTimes[i][0]) + xPosAdjust*(info[1]*1E9)/2))

            # make minor adjustments to height of labels
            if i == 3- 1:
                yPosAdjust = 1.3
            elif i == 7-1:
                yPosAdjust = 0.6
            elif i == 9-1:
                yPosAdjust = 0.6
            elif i == 11-1:
                yPosAdjust = 0.6
            else:
                yPosAdjust = 1

            ylocation = (yPosAdjust)*data_dict['Energy'][0][EnergyHeights[i]]
            ax.text(x=xlocation, y=ylocation, s=textlabel, color='white', ha='center', fontsize=12)

        plt.savefig(r'C:\Users\cfelt\OneDrive\Desktop\AlfvenSignature.png')
        Done(start_time)
        # plt.show()





# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
wflyer = 0

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        IdentificationPlot(0, rocketFolderPath, justPrintFileNames, wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            IdentificationPlot(fileNo, rocketFolderPath, justPrintFileNames, wflyer)
    else:
        for filesNo in wFiles:
            IdentificationPlot(filesNo, rocketFolderPath, justPrintFileNames,wflyer)