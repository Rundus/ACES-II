# --- traj_to_auroral_coordinates.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert the payload ENU coordinates to auroral coordinates by first
# converting them to field-aligned and then rotating them by the amount indicated in
# src.science.auroral_Coordinates.L2_to_auroral_Coordiantes


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
wFile = 0
wRocket = 5
inputPath_modifier = 'trajectories' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = r'science\auroral_coordinates' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# ---------------------------
outputData = True
# ---------------------------

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.interpolate import CubicSpline
import spaceToolsLib as stl



def traj_to_auroral_coordinates(wflyer, wFile, justPrintFileNames):

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
    data_dict_traj, globalAttrsMod  = stl.loadDictFromFile(inputFiles[wFile], getGlobalAttrs=True)
    data_dict_attitude = stl.loadDictFromFile(rf'C:\Data\ACESII\attitude\{ACESII.fliers[wflyer]}\ACESII_{rocketID}_Attitude_Solution.cdf')
    data_dict_auroral_transform = stl.loadDictFromFile(r'C:\Data\ACESII\science\auroral_coordinates\low\ACESII_36364_E_Field_Auroral_Coordinates.cdf')
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
                        'Epoch': deepcopy(data_dict_traj['Epoch']),
                        'Alt': deepcopy(data_dict_traj['Alt']),
                        'Lat': deepcopy(data_dict_traj['Lat']),
                        'Long': deepcopy(data_dict_traj['Long']),
                        'Vel_T': deepcopy(data_dict_traj['ECEFXVEL']),
                        'Vel_N': deepcopy(data_dict_traj['ECEFYVEL']),
                        'Vel_p': deepcopy(data_dict_traj['ECEFZVEL']),
                        'Vel_auroral_perp': deepcopy(data_dict_traj['ECEFZVEL']),
                        'dx_T': deepcopy(data_dict_traj['ECEFXPOS']),
                        'dx_N': deepcopy(data_dict_traj['ECEFXPOS']),
                        'dx_p': deepcopy(data_dict_traj['ECEFXPOS']),
                        }

    #################################################
    # --- INTERPOLATE DATA ONTO ATTITUDE SOLUTION ---
    #################################################
    stl.prgMsg('Interpolating Data on attitude')
    # cubic spline interpolate the trajectory data then evaluate on the Attitude solution Epoch
    Epoch_tt2000_traj = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_traj['Epoch'][0]])
    Epoch_tt2000_attitude = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_attitude['Epoch'][0]])
    data_dict_traj_interp = {}
    for key, val in data_dict_traj.items():
        if key in ['Alt', 'Lat', 'Long']:
            data = deepcopy(data_dict_traj[key][0])
            cs = CubicSpline(x=Epoch_tt2000_traj, y=data)
            data_dict_output[key][0] = cs(Epoch_tt2000_attitude)
        elif key in ['ECEFXPOS', 'ECEFYPOS', 'ECEFZPOS', 'ECEFXVEL', 'ECEFYVEL', 'ECEFZVEL']:
            data = deepcopy(data_dict_traj[key][0])
            cs = CubicSpline(x=Epoch_tt2000_traj,y=data)
            data_dict_traj_interp = {**data_dict_traj_interp,
                                        **{key: [cs(Epoch_tt2000_attitude), data_dict_traj[key][1]]}
                                        }

    data_dict_output['Epoch'][0] = deepcopy(data_dict_attitude['Epoch'][0])
    stl.Done(start_time)

    #####################
    # --- ECEF TO ENU ---
    #####################
    # Form the ECEF vector
    ECEF_POS = np.array([data_dict_traj_interp['ECEFXPOS'][0], data_dict_traj_interp['ECEFYPOS'][0], data_dict_traj_interp['ECEFZPOS'][0]]).T
    ECEF_VEL = np.array([data_dict_traj_interp['ECEFXVEL'][0], data_dict_traj_interp['ECEFYVEL'][0], data_dict_traj_interp['ECEFZVEL'][0]]).T

    ##############################################################
    # --- Convert ECEF coordinates to Field Aligned coordinates ---
    ##############################################################
    stl.prgMsg(f'Loading CHAOS model')

    # Get the Data
    B_model = stl.CHAOS(lat=data_dict_output['Lat'][0],
                    long=data_dict_output['Long'][0],
                    alt=data_dict_output['Alt'][0],
                    times=data_dict_output['Epoch'][0])  # CHAOS in ENU coordinates

    # --- Convert B-Data to GEO (ECEF) XYZ coordinates ---
    ENUtoECEF_transform = np.array([stl.ENUtoECEF(Lat=data_dict_output['Lat'][0][i], Long=data_dict_output['Long'][0][i]) for i in range(len(data_dict_output['Epoch'][0]))])
    B_CHAOS_ECEF = np.array([np.matmul(ENUtoECEF_transform[i], B_model[i]) for i in range(len(data_dict_output['Epoch'][0]))])

    # --- determine the Payload's Position Vector in GEO (ECEF) coordinate XYZ ---
    R_REF = 6371.2  # earth Radius in km
    Radius = deepcopy(data_dict_output['Alt'][0])/1000 + R_REF
    coLatRad = [np.radians(90 - lat) for lat in data_dict_output['Lat'][0]]
    LongRad = [np.radians(long) for long in data_dict_output['Long'][0]]
    Rsc = np.array([
        [Radius[i] * np.sin(coLatRad[i]) * np.cos(LongRad[i]),
         Radius[i] * np.sin(coLatRad[i]) * np.sin(LongRad[i]),
         Radius[i] * np.cos(coLatRad[i])] for i in range(len(data_dict_output['Epoch'][0]))])
    stl.Done(start_time)

    # --- calculate Field Aligned unit vectors over the duration of the flight ---
    stl.prgMsg('Converting to Field Aligned Coordinates')

    # pHat comes from the CHAOS model direction of B in GEO
    pHat = np.array([B_CHAOS_ECEF[i] / np.linalg.norm(B_CHAOS_ECEF[i]) for i in range(len(data_dict_output['Epoch'][0]))])

    # e-hat comes from the cross of pHat and the Rocket's radius vector (in geomagnetic coordinates)
    eHat = np.array([np.cross(pHat[i], Rsc[i]) / np.linalg.norm(np.cross(pHat[i], Rsc[i])) for i in range(len(data_dict_output['Epoch'][0]))])

    # rHat comes from the cross of eHat and pHat
    rHat = np.array([np.cross(eHat[i], pHat[i]) for i in range(len(data_dict_output['Epoch'][0]))])

    # form the transformation matrix FROM ECEF TO FIELD ALIGNED
    FAC_transform = np.array([[eHat[i], pHat[i], rHat[i]] for i in range(len(data_dict_output['Epoch'][0]))])

    # --- Transform Rocket Coordinates into EPR ---
    FAC_VEL = np.array([np.matmul(mat, vec) for mat, vec in zip(FAC_transform, ECEF_VEL)])
    FAC_POS = np.array([np.matmul(mat, vec) for mat, vec in zip(FAC_transform, ECEF_POS)])

    stl.Done(start_time)

    #######################################################
    # --- ROTATE EPR coordinates to Auroral Coordinates ---
    #######################################################
    stl.prgMsg('Rotating into Auroral coordinates')
    rotation_matrix_auroral = np.array([ stl.Ry(data_dict_auroral_transform['rotation_Angle'][0]) for i in range(len(FAC_VEL))])

    VEL_auroral_coordinates = np.array([np.matmul(mat, vec) for mat, vec in zip(rotation_matrix_auroral, FAC_VEL)])
    POS_auroral_coordinates= np.array([np.matmul(mat, vec) for mat, vec in zip(rotation_matrix_auroral, FAC_POS)])

    # store the outputs
    data_dict_output['dx_T'][0] = POS_auroral_coordinates[:, 0] - POS_auroral_coordinates[:, 0][0]
    data_dict_output['dx_N'][0] = POS_auroral_coordinates[:, 1] - POS_auroral_coordinates[:, 0][1]
    data_dict_output['dx_p'][0] = POS_auroral_coordinates[:, 2] - POS_auroral_coordinates[:, 0][2]

    data_dict_output['Vel_T'][0] = VEL_auroral_coordinates[:, 0]
    data_dict_output['Vel_N'][0] = VEL_auroral_coordinates[:, 1]
    data_dict_output['Vel_p'][0] = VEL_auroral_coordinates[:, 2]

    data_dict_output['Vel_auroral_perp'][0] = np.linalg.norm(np.array([VEL_auroral_coordinates[:, 0], VEL_auroral_coordinates[:, 1]]).T, axis=1)
    stl.Done(start_time)

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

        for key in data_dict_output.keys():
            if 'Vel' in key or 'Vperp' in key:
                data_dict_output[key][1]['UNITS'] = 'm/s'
                data_dict_output[key][1]['LABLAXIS'] = key
            elif 'dx_' in key:
                data_dict_output[key][1]['LABLAXIS'] = key

        fileoutName = f'ACESII_{rocketID}_auroral_coordinates_rkt_velocity.cdf'
        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{ACESII.fliers[wflyer]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output, globalAttrsMod=globalAttrsMod,instrNam='Traj')
        stl.Done(start_time)






# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

rocketFolderPath = DataPaths.ACES_data_folder

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}\*.cdf')) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    traj_to_auroral_coordinates(wRocket-4, wFile, justPrintFileNames)

