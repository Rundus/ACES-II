# --- traj_ECEF_to_ENU_to_Field_Aligned.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Input the ACES-II trajectory GPS data and convert it to ENU and Field-aligned coordinates


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
wRocket = 5
inputPath_modifier = 'trajectories' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\cross_track_velocities' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# ---------------------------
outputData = True
# ---------------------------

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.interpolate import CubicSpline



def Traj_ECEF_to_ENU_to_Field_Aligned(wflyer, file_idx, justPrintFileNames):

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
    data_dict_traj, globalAttrsMod  = stl.loadDictFromFile(inputFiles[file_idx], getGlobalAttrs=True)
    data_dict_attitude = stl.loadDictFromFile(rf'C:\Data\ACESII\attitude\{ACESII.fliers[wflyer]}\ACESII_{rocketID}_Attitude_Solution.cdf')
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {'ECEFXVEL': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'ECEFYVEL': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'ECEFZVEL': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'ECEFXPOS': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'ECEFYPOS': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'ECEFZPOS': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],

                        'ENU_E_VEL': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'ENU_N_VEL': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'ENU_U_VEL': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'ENU_E_POS': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'ENU_N_POS': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'ENU_U_POS': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'ENU_Vperp': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],

                        'FAC_E_VEL': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'FAC_P_VEL': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'FAC_R_VEL': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'FAC_E_POS': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'FAC_P_POS': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'FAC_R_POS': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],
                        'FAC_Vperp': [np.zeros(len(data_dict_traj['Epoch'][0])), {}],

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
        if key in ['Alt','Lat','Long','ECEFXPOS','ECEFYPOS','ECEFZPOS','ECEFXVEL','ECEFYVEL','ECEFZVEL']:
            data = deepcopy(data_dict_traj[key][0])
            cs = CubicSpline(x=Epoch_tt2000_traj,y=data)
            data_dict_output[f'{key}'][0] = cs(Epoch_tt2000_attitude)

    data_dict_output['Epoch'][0] = deepcopy(data_dict_attitude['Epoch'][0])
    stl.Done(start_time)

    #####################
    # --- ECEF TO ENU ---
    #####################
    # Form the ECEF vector
    ECEF_POS = np.array([data_dict_output['ECEFXPOS'][0], data_dict_output['ECEFYPOS'][0], data_dict_output['ECEFZPOS'][0]]).T
    ECEF_VEL = np.array([data_dict_output['ECEFXVEL'][0], data_dict_output['ECEFYVEL'][0], data_dict_output['ECEFZVEL'][0]]).T

    # --- form the rotation matrix ---
    ECEFtoENU_transform = np.array([stl.ENUtoECEF(Lat=data_dict_output['Lat'][0][i], Long=data_dict_output['Long'][0][i]).T for i in range(len(data_dict_output['Epoch'][0]))])
    ENU_VEL = np.array([np.matmul(mat, vec) for mat, vec in zip(ECEFtoENU_transform, ECEF_VEL)])
    ENU_POS = np.array([np.matmul(mat, vec) for mat, vec in zip(ECEFtoENU_transform, ECEF_POS)])

    # store the outputs
    data_dict_output['ENU_E_POS'][0] = ENU_POS[:, 0]
    data_dict_output['ENU_N_POS'][0] = ENU_POS[:, 1]
    data_dict_output['ENU_U_POS'][0] = ENU_POS[:, 2]

    data_dict_output['ENU_E_VEL'][0] = ENU_VEL[:, 0]
    data_dict_output['ENU_N_VEL'][0] = ENU_VEL[:, 1]
    data_dict_output['ENU_U_VEL'][0] = ENU_VEL[:, 2]

    data_dict_output['ENU_Vperp'][0] = np.linalg.norm(np.array([ENU_VEL[:, 0], ENU_VEL[:, 1]]).T, axis=1)


    ##############################################################
    # --- Convert ENU coordinates to Field Aligned coordinates ---
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

    # store the outputs
    data_dict_output['FAC_E_POS'][0] = FAC_POS[:, 0]
    data_dict_output['FAC_P_POS'][0] = FAC_POS[:, 1]
    data_dict_output['FAC_R_POS'][0] = FAC_POS[:, 2]

    data_dict_output['FAC_E_VEL'][0] = FAC_VEL[:, 0]
    data_dict_output['FAC_P_VEL'][0] = FAC_VEL[:, 1]
    data_dict_output['FAC_R_VEL'][0] = FAC_VEL[:, 2]

    data_dict_output['FAC_Vperp'][0] = np.linalg.norm(np.array([FAC_VEL[:, 0], FAC_VEL[:, 2]]).T, axis=1)
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


        data_dict_output['ENU_Vperp'][1]['UNITS'] = 'm/s'
        data_dict_output['ENU_E_VEL'][1]['UNITS'] = 'm/s'
        data_dict_output['ENU_N_VEL'][1]['UNITS'] = 'm/s'
        data_dict_output['ENU_U_VEL'][1]['UNITS'] = 'm/s'
        data_dict_output['FAC_Vperp'][1]['UNITS'] = 'm/s'
        data_dict_output['FAC_E_VEL'][1]['UNITS'] = 'm/s'
        data_dict_output['FAC_P_VEL'][1]['UNITS'] = 'm/s'
        data_dict_output['FAC_R_VEL'][1]['UNITS'] = 'm/s'

        fileoutName = f'ACESII_{rocketID}_cross_track_vel.cdf'
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
    file_idx = 0
    Traj_ECEF_to_ENU_to_Field_Aligned(wRocket-4, file_idx, justPrintFileNames)
