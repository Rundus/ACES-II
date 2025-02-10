# --- L_shell.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: calculate the dipolar L-Shell traversed by the rocket throughout its journey



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
wRocket = 5

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

modifier = ''
inputPath_modifier = 'attitude' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = '\science\L_shell' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from spacepy import coordinates as coord
coord.DEFAULTS.set_values(use_irbem=False, itol=5)  # maximum separation, in seconds, for which the coordinate transformations will not be recalculated. To force all transformations to use an exact transform for the time, set ``itol`` to zero.
from spacepy.time import Ticktock #used to determine the time I'm choosing the reference geomagentic field

def ENU_to_Field_Aligned(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]

    fileoutName = f'ACESII_{rocketID}_Lshell.cdf'


    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:

        # --- get the data from the attitude file ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict_attitude = loadDictFromFile(inputFiles[wFile], {})
        data_dict_attitude['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_attitude['Epoch'][0]])
        data_dict_attitude['Alt'][0] = data_dict_attitude['Alt'][0]/1000
        Epoch = data_dict_attitude['Epoch'][0]
        Done(start_time)

        ############################################
        # --- CONVERT TO GEOMAGNETIC COORDINATES ---
        ############################################
        prgMsg('Converting to geomagnetic coordinates')

        # Trajectory Data

        geodeticPos = np.array([ [data_dict_attitude['Alt'][0][i], data_dict_attitude['Latgd'][0][i], data_dict_attitude['Long'][0][i]] for i in range(len(Epoch))])
        ISOtime = [pycdf.lib.tt2000_to_datetime(tme).isoformat() for tme in Epoch]
        cvals_GDZ = coord.Coords(geodeticPos, 'GDZ', 'sph')
        cvals_GDZ.ticks = Ticktock(ISOtime, 'ISO')
        cvals_GDZ_MAG = cvals_GDZ.convert('MAG', 'sph')

        geomagAlt = np.array([cvals_GDZ_MAG[i].radi for i in range(len(Epoch))])
        geomagLat = np.array([cvals_GDZ_MAG[i].lati for i in range(len(Epoch))])
        geomagLong = np.array([cvals_GDZ_MAG[i].long for i in range(len(Epoch))])

        Done(start_time)



        # ###################################################
        # --- Calculate L-Shell Parameter over the flight ---
        # ###################################################
        prgMsg('Calculating L-Shell')
        L_shell = np.array([geomagAlt[i]/((np.cos(np.radians(geomagLat[i])))**2) for i in range(len(data_dict_attitude['Epoch'][0]))])
        Done(start_time)

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('creating output file')

            varAttrs = {'LABLAXIS': 'L-shell','DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None,
                                                       'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                       'UNITS': None,
                                                       'VALIDMIN': L_shell.min(), 'VALIDMAX': L_shell.max(),
                                                       'VAR_TYPE': 'data', 'SCALETYP': 'linear'}
            Epoch_output = deepcopy(data_dict_attitude['Epoch'])
            Epoch_output[1]['VAR_TYPE'] = 'support_data'


            data_dict_output = {'Epoch':Epoch_output,'L-Shell':[L_shell,varAttrs]}

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            outputCDFdata(outputPath, data_dict_output, outputModelData, globalAttrsMod, 'Attitude')

            Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        ENU_to_Field_Aligned(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            ENU_to_Field_Aligned(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            ENU_to_Field_Aligned(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)