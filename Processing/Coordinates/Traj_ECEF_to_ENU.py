# --- Traj_ECEF_to_ENU.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:


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
justPrintFileNames = False
wRocket = 5
modifier = ''
wFileTraj = 0
wFileAttitude = 0
inputPath_modifier = 'trajectories' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\cross_track_velocities' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# ---------------------------
outputData = True
# ---------------------------

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.class_var_func import ENUtoECEF


def Traj_ECEF_to_ENU(wRocket, rocketFolderPath, justPrintFileNames, wflyer, wFileAtt, wFile_Traj):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    ModelData = L2_TRICE_Quick(wflyer)

    inputFiles, input_names, input_names_searchable = getInputFiles(rocketFolderPath=rocketFolderPath,wRocket=wRocket,inputPath_modifier=inputPath_modifier)
    inputFiles_attitude, input_names_attitude, input_names_searchable_attitude = getInputFiles(rocketFolderPath=rocketFolderPath, wRocket=wRocket,inputPath_modifier='attitude')

    fileoutName = f'ACESII_{rocketID}_crossTrackVel.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))

        print('\n')
        for i, file in enumerate(inputFiles_attitude):
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable_attitude[i], round(getsize(file) / (10 ** 6), 1)))
        return

    print('\n')

    # --- --- --- --- --- -
    # --- LOAD THE DATA ---
    # --- --- --- --- --- -

    # --- get the data from the file ---
    prgMsg(f'Loading data from {inputPath_modifier} Files')
    data_dict_traj, globalAttrsMod  = loadDictFromFile(inputFiles[wFile_Traj], getGlobalAttrs=True)
    data_dict_attitude = loadDictFromFile(inputFiles_attitude[wFileAtt])
    Done(start_time)

    exampleAttrs = {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal,
                    'FORMAT': 'E12.2', 'UNITS': None, 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data', 'SCALETYP': 'linear',
                    'LABLAXIS': None}

    # --- --- --- --- ---
    # --- INTERPOLATE ---
    # --- --- --- --- ---
    # interpolate trajectory data onto attitude timebase
    data_dict_traj_interp = InterpolateDataDict(InputDataDict=data_dict_traj,
                                                InputEpochArray=data_dict_traj['Epoch'][0],
                                                targetEpochArray=data_dict_attitude['Epoch'][0],
                                                wKeys=['ECEFXPOS',
                                                       'ECEFYPOS',
                                                       'ECEFZPOS',
                                                       'ECEFXVEL',
                                                       'ECEFYVEL',
                                                       'ECEFZVEL',
                                                       'Lat',
                                                       'Long',
                                                       'Alt',
                                                       'Epoch'])

    # --- --- --- --- --- --- --- --- --
    # --- CREATE THE ROTATION MATRIX ---
    # --- --- --- --- --- --- --- --- --

    # Form the ECEF vector
    ECEF_POS = np.array([data_dict_traj_interp['ECEFXPOS'][0], data_dict_traj_interp['ECEFYPOS'][0], data_dict_traj_interp['ECEFZPOS'][0]]).T
    ECEF_VEL = np.array([data_dict_traj_interp['ECEFXVEL'][0], data_dict_traj_interp['ECEFYVEL'][0], data_dict_traj_interp['ECEFZVEL'][0]]).T

    transformMatrix = np.array([ENUtoECEF(Lat=data_dict_traj_interp['Lat'][0][i], Long=data_dict_traj_interp['Long'][0][i]).T for i in range(len(data_dict_traj_interp['Epoch'][0]))])
    ENU_VEL = np.array([np.matmul(mat, vec) for mat, vec in zip(transformMatrix, ECEF_VEL)])/1000
    ENU_POS = np.array([np.matmul(mat, vec) for mat, vec in zip(transformMatrix, ECEF_POS)])
    ENU_Vperp = np.array([(ENU_VEL[i][0]**2 + ENU_VEL[i][1]**2)**(0.5) for i in range(len(data_dict_traj_interp['Epoch'][0]))])

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        prgMsg('Creating output file')

        data_dict_output = {'ENU_E_VEL' : [ENU_VEL[:, 0], deepcopy(exampleAttrs)],
                            'ENU_N_VEL' : [ENU_VEL[:, 1], deepcopy(exampleAttrs)],
                            'ENU_U_VEL' : [ENU_VEL[:, 2], deepcopy(exampleAttrs)],
                            'ENU_E_POS' : [ENU_POS[:, 0], deepcopy(exampleAttrs)],
                            'ENU_N_POS' : [ENU_POS[:, 1], deepcopy(exampleAttrs)],
                            'ENU_U_POS' : [ENU_POS[:, 2], deepcopy(exampleAttrs)],
                            'ENU_Vperp' : [ENU_Vperp, deepcopy(exampleAttrs)],
                            'Epoch' : deepcopy(data_dict_traj_interp['Epoch']),
                            'Alt' : deepcopy(data_dict_traj_interp['Alt']),
                            'Lat' : deepcopy(data_dict_traj_interp['Lat']),
                            'Long' : deepcopy(data_dict_traj_interp['Long'])}

        data_dict_output['ENU_Vperp'][1]['UNITS'] = 'km/s'
        data_dict_output['ENU_E_VEL'][1]['UNITS'] = 'km/s'
        data_dict_output['ENU_N_VEL'][1]['UNITS'] = 'km/s'
        data_dict_output['ENU_U_VEL'][1]['UNITS'] = 'km/s'

        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

        outputCDFdata(outputPath, data_dict_output, globalAttrsMod=globalAttrsMod,instrNam='Traj')

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
    Traj_ECEF_to_ENU(wRocket, rocketFolderPath, justPrintFileNames,wflyer,wFileAttitude, wFileTraj)
