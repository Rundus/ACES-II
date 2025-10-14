# --- cal1_L1_to_spinUp_vxB_to_rktFrm.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: take the CHAOS geomagnetic field model and
# spin it into the rocket frame. Then calculate |E| and |vxB|. Perform a
# timing analysis on the peaks between vxB to find the temporal offset

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"

import matplotlib.pyplot as plt
import numpy as np
import spaceToolsLib

from src.my_imports import *
import time
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

# --- OutputData ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.my_imports import *
import spaceToolsLib as stl
from scipy.interpolate import CubicSpline

#######################
# --- MAIN FUNCTION ---
#######################
def cal1_L1_to_spinUp_vxB_to_rktFrm(wRocket, justPrintFileNames):

    # --- FILE I/O ---
    data_folder_path = rf'{DataPaths.ACES_data_folder}\\L1\\{ACESII.fliers[wRocket-4]}\\'

    # get the EFI files
    input_files = glob(data_folder_path + '*EFI_rktFrm*')
    input_names = [ifile.replace(data_folder_path, '') for ifile in input_files]

    # get the attitude files
    data_folder_path = rf'{DataPaths.ACES_data_folder}\\attitude\\{ACESII.fliers[wRocket-4]}\\'
    input_files_attitude = glob(data_folder_path +'*.cdf*')[0]

    # get the ACES-II ECEF trajectory files
    data_folder_path = rf'{DataPaths.ACES_data_folder}\trajectories\\{ACESII.fliers[wRocket - 4]}\\'
    input_files_traj = glob(data_folder_path + '*ECEF.cdf*')[0]

    # get the ACES-II RingCore Data
    data_folder_path = rf'{DataPaths.ACES_data_folder}\L1\\{ACESII.fliers[wRocket - 4]}\\'
    input_files_mag = glob(data_folder_path + '*RingCore_rktFrm.cdf*')[0]

    if justPrintFileNames:
        for i, file in enumerate(input_files):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names[i], round(os.path.getsize(file) / (10 ** 6), 1)))
        return

    # --- get the data from the L2 file ---
    stl.prgMsg(f'Loading data from L1 Files')
    data_dict_EFI = stl.loadDictFromFile(input_files[0])
    data_dict_attitude = stl.loadDictFromFile(input_files_attitude)
    data_dict_traj = stl.loadDictFromFile(input_files_traj)
    data_dict_mag = stl.loadDictFromFile(input_files_mag)
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}

    ############################################################
    # --- Interpolate the Attitude DCM onto the EFI timebase ---
    ############################################################

    # Attempt 1: Interpolate the DCM onto EFI timebase.
    stl.prgMsg('Interpolating DCM')
    T0 = dt.datetime(2022, 11, 20, 17, 20)
    T0_EFI = stl.EpochTo_T0_Rocket(data_dict_EFI['Epoch'][0],T0=T0)
    T0_attitude = stl.EpochTo_T0_Rocket(data_dict_attitude['Epoch'][0],T0=T0)
    DCM_rkt_to_ENU = np.zeros(shape=(len(data_dict_EFI['Epoch'][0]), 3, 3))
    for i in range(1, 4):
        for j in range(1, 4):
            cs = CubicSpline(T0_attitude, data_dict_attitude[f'a{i}{j}'][0])
            DCM_rkt_to_ENU[:, i - 1, j - 1] = cs(T0_EFI)
    DCM_ENU_to_rkt = np.array([mat.T for mat in DCM_rkt_to_ENU])
    stl.Done(start_time)

    # get B in ENU coordinates
    stl.prgMsg('Getting B-Field in ENU')
    B_ENU = stl.CHAOS(lat=data_dict_traj['Lat'][0],
                      long=data_dict_traj['Long'][0],
                      alt=data_dict_traj['Alt'][0],
                      times=data_dict_traj['Epoch'][0])
    stl.Done(start_time)


    # Interpolate the trajectory Velocity (ECEF)
    stl.prgMsg('Interpolating Trajectory')
    T0_traject = stl.EpochTo_T0_Rocket(data_dict_traj['Epoch'][0], T0=T0)
    for key in ['ECEFXVEL', 'ECEFYVEL', 'ECEFZVEL', 'Lat', 'Long', 'Alt']:
        cs = CubicSpline(T0_traject, data_dict_traj[key][0])
        data_dict_traj[key][0] = cs(T0_EFI)
    stl.Done(start_time)

    # Interpolate the B-Field onto EFI timebase
    stl.prgMsg('Interpolating B')
    # CHAOS ENU coordinates
    B_ENU_new = [[],[],[]]
    for idx in range(3):
        cs = CubicSpline(T0_traject, B_ENU[:, idx])
        B_ENU_new[idx] = cs(T0_EFI)
    B_ENU = (1E-9)*np.array(B_ENU_new).T

    # RingCore rkt XYZ coordinates
    T0_mag = stl.EpochTo_T0_Rocket(data_dict_mag['Epoch'][0], T0=T0)
    for key in ['Bx', 'By', 'Bz']:
        cs = CubicSpline(T0_mag, data_dict_mag[key][0])
        data_dict_mag[key][0] = cs(T0_EFI)
    stl.Done(start_time)

    # Rotate the Trajectory velocity in ENU
    stl.prgMsg('Rotating Trajectory into RktFrm')
    V_rkt_ECEF = np.array([data_dict_traj['ECEFXVEL'][0], data_dict_traj['ECEFYVEL'][0], data_dict_traj['ECEFZVEL'][0]]).T
    ENU_to_ECEF = np.array([stl.ENUtoECEF(Lat=data_dict_traj['Lat'][0][i], Long=data_dict_traj['Long'][0][i]) for i in range(len(V_rkt_ECEF))])
    V_rkt_ENU = np.array([ENU_to_ECEF[i]@V_rkt_ECEF[i] for i in range(len(V_rkt_ECEF))])
    V_rkt_rktFrm = np.array([DCM_ENU_to_rkt[i]@V_rkt_ENU[i] for i in range(len(V_rkt_ECEF))])
    stl.Done(start_time)

    # Calculate vxB in rocket Coordinates
    stl.prgMsg('Calculating vxB')
    B = 1E-9*(np.array([data_dict_mag['Bx'][0], data_dict_mag['By'][0], data_dict_mag['Bz'][0]]).T)
    vxB_ENU = -1*np.array([np.cross(V_rkt_ENU[i], B_ENU[i]) for i in range(len(V_rkt_ENU))])
    vxB = -1*np.array([np.cross(V_rkt_rktFrm[i], B[i]) for i in range(len(V_rkt_ECEF))])
    vxB_mag = np.array([np.linalg.norm(vxB[i]) for i in range(len(V_rkt_ECEF))])
    E = np.array([data_dict_EFI['E_x'][0], data_dict_EFI['E_x'][0], data_dict_EFI['E_z'][0]]).T
    E_mag = np.array([np.linalg.norm(E[i]) for i in range(len(V_rkt_ECEF))])
    stl.Done(start_time)

    data_dict_output = {**data_dict_output,
                        **{
                            'Epoch': deepcopy(data_dict_EFI['Epoch']),
                            'vxB_E': [np.array(vxB_ENU[:, 0]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS': 'vxB_E'}],
                            'vxB_N': [np.array(vxB_ENU[:, 1]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS': 'vxB_N'}],
                            'vxB_Up': [np.array(vxB_ENU[:, 2]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS': 'vxB_Up'}],
                            'vxB_X': [np.array(vxB[:, 0]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS':'vxB_X'}],
                            'vxB_Y': [np.array(vxB[:, 1]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS':'vxB_Y'}],
                            'vxB_Z': [np.array(vxB[:, 2]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS':'vxB_Z'}],
                            '|vxB|': [np.array(vxB_mag), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS':'|vxB|'}],
                            'E_X': deepcopy(data_dict_EFI['E_x']),
                            'E_Y': deepcopy(data_dict_EFI['E_y']),
                            'E_Z': deepcopy(data_dict_EFI['E_z']),
                            '|E|': [np.array(E_mag), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m','LABLAXIS':'|E|'}],
                            'DCM': [DCM_rkt_to_ENU, {'LABLAXIS':'DCM_rkt_to_ENU'}],
                        }
                        }


    # Include the DCM elements in the data_dict_output
    data_dict_output

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:
        stl.prgMsg('Creating output file')
        fileoutName = f'ACESII_{ACESII.payload_IDs[wRocket-4]}_EFI_vxB_rktFrm.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\\calibration\\EFI_cal1_spin_vxB_into_rktFrm_calibration\\{ACESII.fliers[wRocket-4]}\\'+ fileoutName
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
cal1_L1_to_spinUp_vxB_to_rktFrm(wRocket, justPrintFileNames)

