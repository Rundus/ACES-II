# --- ACESII_EISCAT_DensitySlice.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Uses the ACESII attitde data to determine an expected density profile for the
# Rockets

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt

from ACESII_code.myImports import *
from scipy.interpolate import CubicSpline
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
wRocket = 4
modifier = ''
useTromso = True # if false USE SVALBARD
inputPath_modifier = r'science\EISCAT\tromso\UHF' if useTromso else r'science\EISCAT\svalbard\UKFI_radar' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\EISCAT_ACESII_Slice' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
# ---------------------------
outputData = True
# ---------------------------

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none


def ACESII_EISCAT_DensitySlice(wRocket, rocketFolderPath):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wRocket-4]
    globalAttrsMod = rocketAttrs.globalAttributes[wRocket-4]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    inputFiles_EISCAT = glob(rf'{rocketFolderPath}{inputPath_modifier}\*.cdf')[0]
    inputFiles_attitude = glob(f'{rocketFolderPath}attitude\{fliers[wRocket-4]}{modifier}\*.cdf')[0]
    fileoutName = f'ACESII_{rocketID}_EISCAT_Tromso_rktSlice.cdf' if useTromso else f'ACESII_{rocketID}_EISCAT_Svalbard_rktSlice.cdf'


    # --- --- --- --- --- -
    # --- LOAD THE DATA ---
    # --- --- --- --- --- -

    # --- get the data from the file ---
    prgMsg(f'Loading data EISCAT Files')
    data_dict_attitude = loadDictFromFile(inputFiles_attitude)
    data_dict_EISCAT = loadDictFromFile(inputFiles_EISCAT, targetVar=[[data_dict_attitude['Epoch'][0][0],data_dict_attitude['Epoch'][0][-1]],'Epoch'])
    Done(start_time)

    exampleAttrs = {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
     'UNITS': None, 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data', 'SCALETYP': 'linear',
     'LABLAXIS': None}

    # --- --- --- --- --- --- --- --
    # --- INTERPOLATE EISCAT DATA---
    # --- --- --- --- --- --- --- --
    # description: For each slice in time, take all the real values of ne, Ti, T_ratio, Nu_in, Ion_Comp, Op_Comp
    # and spline interpolate them. Evaluate them on a new altitude base (HF/LF altitude). Then for each Rkt epoch val
    # find the nearest relevant value to keep
    prgMsg('Interpolating EISCAT')

    # --- determine the range of altitude to interpolate over ---
    minVal = 70 if wRocket == 4 else 90
    maxVal = 500 if wRocket == 4 else 200

    altMin = np.abs(data_dict_attitude['Alt'][0]/1000 - minVal).argmin() # in kilometers
    altMax = np.abs(data_dict_attitude['Alt'][0]/1000 - maxVal).argmin() # in kilometers
    altForInterp = deepcopy(data_dict_attitude['Alt'][0][altMin:altMax+1]/1000)

    ######################################
    # --- EISCAT Spatial Interpolation ---
    ######################################

    data_dict_interp = {
        'ne': [[], deepcopy(exampleAttrs)],
        'ti': [[], deepcopy(exampleAttrs)],
        'tr': [[], deepcopy(exampleAttrs)],
        'co': [[], deepcopy(exampleAttrs)],
        'pm': [[], deepcopy(exampleAttrs)], # The number of ions with molecular weight 28 to 32 AS A FRACTION OF TOTAL ELECTRON DENSITY (# of ions/Ne)
        'po+': [[], deepcopy(exampleAttrs)], # composition [O+]/Ne
    }

    for key in data_dict_interp.keys():
        for tme, EpochVal in tqdm(enumerate(data_dict_EISCAT['Epoch'][0])):
            # get each variable and ONLY the real data
            indicies = np.where(np.isnan(data_dict_EISCAT[key][0][tme])==False)[0]
            val = data_dict_EISCAT[key][0][tme][indicies]
            Alt = data_dict_EISCAT['gdalt'][0][tme][indicies]

            splCub = CubicSpline(Alt, val)
            data_dict_interp[key][0].append(np.array([splCub(altRktVal) for altRktVal in altForInterp]))

    # add in the Epoch and Range variable, then convert everything to numpy array
    data_dict_interp = {**data_dict_interp, **{'Epoch':[deepcopy(data_dict_EISCAT['Epoch'][0]),deepcopy(exampleAttrs)], 'range':[altForInterp,deepcopy(exampleAttrs)]}}
    for key in data_dict_interp.keys():
        data_dict_interp[key][0] = np.array(data_dict_interp[key][0])


    # --- Sample the interpolated Values on ACESII Timebase ---
    data_dict = {
        'ne': [ np.zeros(shape=(len(data_dict_attitude['Epoch'][0]))),deepcopy(exampleAttrs)],
        'Ti': [np.zeros(shape=(len(data_dict_attitude['Epoch'][0]))), deepcopy(exampleAttrs)],
        'T_ratio': [np.zeros(shape=(len(data_dict_attitude['Epoch'][0]))), deepcopy(exampleAttrs)],
        'Nu_in': [np.zeros(shape=(len(data_dict_attitude['Epoch'][0]))), deepcopy(exampleAttrs)],
        'Ion_Comp': [np.zeros(shape=(len(data_dict_attitude['Epoch'][0]))), deepcopy(exampleAttrs)], # The number of ions with molecular weight 28 to 32 AS A FRACTION OF TOTAL ELECTRON DENSITY (# of ions/Ne)
        'Op_Comp': [np.zeros(shape=(len(data_dict_attitude['Epoch'][0]))), deepcopy(exampleAttrs)], # composition [O+]/Ne
    }

    Done(start_time)

    # loop over epoch of rocket
    prgMsg('Sampling Tromso EISCAT Data')
    for tme in range(len(data_dict_attitude['Epoch'][0])):

        # get the rockets altitude and epoch
        altRkt = data_dict_attitude['Alt'][0][tme]/1000
        epochValRkt = data_dict_attitude['Epoch'][0][tme]

        # get the relevant point in the EISCAT data
        closestAlt_Index = np.abs(data_dict_interp['range'][0] - altRkt).argmin()
        closestTime_Index = np.abs(data_dict_interp['Epoch'][0] - epochValRkt).argmin()
        closestAlt = data_dict_interp['range'][0][closestAlt_Index]

        # check if we're 1km away from the cloests EISCAT point. If so, don't include this point
        if np.abs(altRkt - closestAlt) > 10:

            # set everything to zero if we're too far away
            for key in data_dict.keys():
                data_dict[key][0][tme] = 0

        else:
            data_dict['ne'][0][tme] = data_dict_interp['ne'][0][closestTime_Index][closestAlt_Index]
            data_dict['Ti'][0][tme] = (kB/q0)*data_dict_interp['ti'][0][closestTime_Index][closestAlt_Index]
            data_dict['T_ratio'][0][tme] = data_dict_interp['tr'][0][closestTime_Index][closestAlt_Index]
            data_dict['Nu_in'][0][tme] = data_dict_interp['co'][0][closestTime_Index][closestAlt_Index]
            data_dict['Ion_Comp'][0][tme] = data_dict_interp['pm'][0][closestTime_Index][closestAlt_Index]
            data_dict['Op_Comp'][0][tme] = data_dict_interp['po+'][0][closestTime_Index][closestAlt_Index]

    Done(start_time)
    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        prgMsg('Creating output file')


        # update the data_dict with other known variables
        data_dict = {**data_dict, **{
            'Epoch': deepcopy(data_dict_attitude['Epoch']),
            'Alt' : deepcopy(data_dict_attitude['Alt']),
            'Long': deepcopy(data_dict_attitude['Long']),
            'Lat': deepcopy(data_dict_attitude['Lat']),
            'ILat': deepcopy(data_dict_attitude['ILat']),
            'ILong': deepcopy(data_dict_attitude['ILong'])
        }}

        # update the attributs
        data_dict['ne'][0] = data_dict['ne'][0]/(cm_to_m**3)
        data_dict['ne'][1]['LABLAXIS'] ='Density'
        data_dict['ne'][1]['UNITS'] = 'cm!A-3!N'

        data_dict['Ti'][1]['LABLAXIS'] = 'Ion Temperature'
        data_dict['Ti'][1]['UNITS'] = 'eV'

        data_dict['T_ratio'][1]['LABLAXIS'] = 'Temperature Ratio'
        data_dict['T_ratio'][1]['UNITS'] = 'Te/Ti'

        data_dict['Nu_in'][1]['LABLAXIS'] = 'Ion-Neutral Collision Freq'
        data_dict['Nu_in'][1]['UNITS'] = 'Hz'

        data_dict['Ion_Comp'][1]['LABLAXIS'] = 'Ion Weight Ratio'
        data_dict['Ion_Comp'][1]['UNITS'] = '(Ions mol 28 to 32)/Ne'

        data_dict['Op_Comp'][1]['LABLAXIS'] = 'Oxygen Ratio'
        data_dict['Op_Comp'][1]['UNITS'] = '[O+]/Ne'

        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wRocket-4]}\\{fileoutName}'

        outputCDFdata(outputPath, data_dict, instrNam='EISCAT')

        Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
if len(glob(rf'{rocketFolderPath}\science\EISCAT\tromso\UHF\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    ACESII_EISCAT_DensitySlice(wRocket, rocketFolderPath)
