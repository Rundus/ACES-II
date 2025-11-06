# --- ECEF_to_ENU_transformation_matrix.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Input the ACES-II trajectory GPS data and convert it to ENU and Field-aligned coordinates

# OUTPUT:
# [X_ECEF, Y_ECEF, Z_ECEF] *transform_matrix = [E, N, U]

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np
from src.my_imports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False
wRocket = 4
inputPath_modifier = 'trajectories' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'coordinates\\transforms' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# ---------------------------
outputData = True
# ---------------------------

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.interpolate import CubicSpline
def ENU_coordinates_transformation_matrix(wflyer, file_idx, justPrintFileNames):

    # --- FILE I/O ---
    rocket_folder_path = DataPaths.ACES_data_folder
    rocketID = ACESII.payload_IDs[wflyer]

    # get the input data
    inputFiles, input_names, input_names_searchable = stl.getInputFiles(rocketFolderPath=rocketFolderPath,
                                                                        wRocket=wRocket,
                                                                        inputPath_modifier=inputPath_modifier)

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return

    # --- --- --- --- --- -
    # --- LOAD THE DATA ---
    # --- --- --- --- --- -

    # --- get the data from the file ---
    stl.prgMsg(f'Loading data from {inputPath_modifier} Files')
    data_dict_traj, globalAttrsMod = stl.loadDictFromFile(inputFiles[file_idx], getGlobalAttrs=True)
    data_dict_attitude = stl.loadDictFromFile(rf'C:\Data\ACESII\attitude\{ACESII.fliers[wflyer]}\ACESII_{rocketID}_Attitude_Solution_fullCal.cdf')
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
                        'a11': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'a12': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'a13': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'a21': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'a22': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'a23': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'a31': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'a32': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'a33': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'ECEF_to_ENU':[np.zeros(shape=(len(data_dict_traj['Epoch'][0]),3,3)), {}],
                        'ENU_to_ECEF': [np.zeros(shape=(len(data_dict_traj['Epoch'][0]), 3, 3)), {}],
                        'Epoch': deepcopy(data_dict_traj['Epoch']),
                        'Alt': deepcopy(data_dict_traj['Alt']),
                        'Lat': deepcopy(data_dict_traj['Lat']),
                        'Long': deepcopy(data_dict_traj['Long'])
                        }

    #################################################
    # --- INTERPOLATE DATA ONTO ATTITUDE SOLUTION ---
    #################################################
    stl.prgMsg('Interpolating Data on attitude')
    # cubic spline interpolate the trajectory data then evaluate on the Attitude solution Epoch
    Epoch_tt2000_traj = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_traj['Epoch'][0]])
    Epoch_tt2000_attitude = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_attitude['Epoch'][0]])
    for key, val in data_dict_traj.items():
        if key in ['Alt','Lat','Long']:
            data = deepcopy(data_dict_traj[key][0])
            cs = CubicSpline(x=Epoch_tt2000_traj,y=data)
            data_dict_output[f'{key}'][0] = cs(Epoch_tt2000_attitude)

    data_dict_output['Epoch'][0] = deepcopy(data_dict_attitude['Epoch'][0])
    stl.Done(start_time)

    #####################
    # --- ECEF TO ENU ---
    #####################

    # --- form the rotation matrix ---
    ECEFtoENU_transform = np.array([stl.ENUtoECEF(Lat=data_dict_output['Lat'][0][i], Long=data_dict_output['Long'][0][i]).T for i in range(len(data_dict_output['Epoch'][0]))])

    # store the outputs
    data_dict_output['a11'][0] = ECEFtoENU_transform[:, 0, 0]
    data_dict_output['a12'][0] = ECEFtoENU_transform[:, 0, 1]
    data_dict_output['a13'][0] = ECEFtoENU_transform[:, 0, 2]
    data_dict_output['a21'][0] = ECEFtoENU_transform[:, 1, 0]
    data_dict_output['a22'][0] = ECEFtoENU_transform[:, 1, 1]
    data_dict_output['a23'][0] = ECEFtoENU_transform[:, 1, 2]
    data_dict_output['a31'][0] = ECEFtoENU_transform[:, 2, 0]
    data_dict_output['a32'][0] = ECEFtoENU_transform[:, 2, 1]
    data_dict_output['a33'][0] = ECEFtoENU_transform[:, 2, 2]
    data_dict_output['ECEF_to_ENU'][0] = ECEFtoENU_transform
    data_dict_output['ENU_to_ECEF'][0] = np.array([ECEFtoENU_transform[i].T for i in range(len(data_dict_output['Epoch'][0]))])

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:
        stl.prgMsg('Creating output file')

        exampleAttrs = {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': ACESII.epoch_fillVal,
                        'FORMAT': 'E12.2', 'UNITS': None, 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data',
                        'SCALETYP': 'linear',
                        'LABLAXIS': None}


        # update the data dict attrs
        for key, val in data_dict_output.items():
            newAttrs = deepcopy(exampleAttrs)

            for subKey, subVal in data_dict_output[key][1].items():
                newAttrs[subKey] = subVal

            data_dict_output[key][1] = newAttrs


        fileoutName = f'ACESII_{rocketID}_ECEF_to_ENU.cdf'
        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{ACESII.fliers[wflyer]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output, globalAttrsMod=globalAttrsMod)
        stl.Done(start_time)






# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

rocketFolderPath = DataPaths.ACES_data_folder

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}\*.cdf')) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    ENU_coordinates_transformation_matrix(wRocket-4, 0, justPrintFileNames)
