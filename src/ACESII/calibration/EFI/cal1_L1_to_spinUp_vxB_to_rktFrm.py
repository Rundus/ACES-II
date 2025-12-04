# --- cal1_L1_to_spinUp_vxB_to_rktFrm.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: take the CHAOS geomagnetic field model and
# spin it into the rocket frame. Then calculate |E| and |vxB|. Perform a
# timing analysis on the peaks between vxB to find the temporal offset

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"

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
from src.ACESII.my_imports import *
import spaceToolsLib as stl


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
    input_files_traj = glob(data_folder_path + '*ENU.cdf*')[0]

    # get the ACES-II RingCore Data - Rocket Frame
    data_folder_path = rf'{DataPaths.ACES_data_folder}\L1\\{ACESII.fliers[wRocket - 4]}\\'
    input_files_mag = glob(data_folder_path + '*RingCore_rktFrm.cdf*')[0] # units: nT

    # get the ACES-II RingCore Data - ENU Frame
    data_folder_path = rf'{DataPaths.ACES_data_folder}\L2\\{ACESII.fliers[wRocket - 4]}\\'
    input_files_mag_ENU = glob(data_folder_path + '*RingCore_ENU.cdf*')[0]  # units: nT

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
    data_dict_mag_ENU = stl.loadDictFromFile(input_files_mag_ENU)
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}

    ############################################################
    # --- Interpolate the Attitude DCM onto the EFI timebase ---
    ############################################################

    # Attempt 1: Interpolate the DCM onto EFI timebase.
    stl.prgMsg('Interpolating DCM')
    T0_EFI = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_EFI['Epoch'][0]])
    T0_attitude = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_attitude['Epoch'][0]])
    DCM_rkt_to_ENU = np.zeros(shape=(len(data_dict_EFI['Epoch'][0]), 3, 3))
    for i in range(1, 4):
        for j in range(1, 4):
            DCM_rkt_to_ENU[:, i - 1, j - 1] = np.interp(T0_EFI, T0_attitude, data_dict_attitude[f'a{i}{j}'][0])

    DCM_ENU_to_rkt = np.array([mat.T for mat in DCM_rkt_to_ENU])
    stl.Done(start_time)

    # get B in ENU coordinates
    stl.prgMsg('Getting B-Field in ENU')
    B_ENU = stl.CHAOS(lat=data_dict_traj['Lat'][0],
                      long=data_dict_traj['Long'][0],
                      alt=data_dict_traj['Alt'][0],
                      times=data_dict_traj['Epoch'][0])
    stl.Done(start_time)

    # Interpolate the trajectory Velocity (ENU) and B-Field onto EFI timebase
    stl.prgMsg('Interpolating Trajectory')
    T0_traject = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_traj['Epoch'][0]])
    for key in ['E_VEL', 'N_VEL', 'U_VEL', 'Lat', 'Long', 'Alt']:
        data_dict_traj[key][0] = np.interp(T0_EFI,T0_traject, data_dict_traj[key][0])
    stl.Done(start_time)

    stl.prgMsg('Interpolating B')
    temp_ENU_new = [[],[],[]] # CHAOS ENU coordinates
    for idx in range(3):
        temp_ENU_new[idx] = np.interp(T0_EFI,T0_traject, B_ENU[:, idx])
    B_ENU = (1E-9)*(np.array(temp_ENU_new).T)

    # # RingCore rkt XYZ coordinates
    # T0_mag = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_mag['Epoch'][0]])
    # for key in ['Bx', 'By', 'Bz']:
    #     data_dict_mag[key][0] = np.interp(T0_EFI,T0_mag, data_dict_mag[key][0])

    # RingCore rkt ENU coordinates
    T0_mag = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_mag_ENU['Epoch'][0]])
    for key in ['B_East', 'B_North', 'B_Up','B_model_East','B_model_North','B_model_Up']:
        data_dict_mag_ENU[key][0] = 1E-9*np.interp(T0_EFI, T0_mag, data_dict_mag_ENU[key][0])

    # Rotate the Trajectory velocity in ENU
    V_rkt_ENU = np.array([data_dict_traj['E_VEL'][0], data_dict_traj['N_VEL'][0], data_dict_traj['U_VEL'][0]]).T
    V_rkt_rktFrm = np.array([np.matmul(DCM_ENU_to_rkt[i], V_rkt_ENU[i]) for i in range(len(V_rkt_ENU))])
    stl.Done(start_time)

    # Calculate vxB in rocket Coordinates
    stl.prgMsg('Calculating vxB')
    B_rkt = np.array([np.matmul(DCM_ENU_to_rkt[i], B_ENU[i]) for i in range(len(T0_EFI))])
    vxB_rkt = np.array([np.cross(V_rkt_rktFrm[i], B_rkt[i]) for i in range(len(T0_EFI))])
    vxB_mag_rkt = np.array([np.linalg.norm(vxB_rkt[i]) for i in range(len(T0_EFI))])

    # Calculate vxB in ENU coordinates (model)
    vxB_ENU = np.array([np.cross(V_rkt_ENU[i], B_ENU[i]) for i in range(len(T0_EFI))])
    vxB_mag_ENU = np.array([np.linalg.norm(vxB_ENU[i]) for i in range(len(T0_EFI))])

    # Calculate vxB in ENU coordinates from RingCore (Not model)
    B_RingCore = np.array([
                            data_dict_mag_ENU['B_East'][0]+data_dict_mag_ENU['B_model_East'][0],
                            data_dict_mag_ENU['B_North'][0] + data_dict_mag_ENU['B_model_North'][0],
                            data_dict_mag_ENU['B_Up'][0] + data_dict_mag_ENU['B_model_Up'][0],
                           ]).T

    vxB_ENU_RingCore = np.array([np.cross(V_rkt_ENU[i], B_RingCore[i]) for i in range(len(T0_EFI))])


    # Calculate |E|, |B| (diagnostic)
    E = np.array([data_dict_EFI['E_x'][0], data_dict_EFI['E_y'][0], data_dict_EFI['E_z'][0]]).T
    E_mag = np.array([np.linalg.norm(E[i]) for i in range(len(T0_EFI))])
    B_mag = np.array([np.linalg.norm(B_rkt[i]) for i in range(len(T0_EFI))])
    stl.Done(start_time)

    data_dict_output = {**data_dict_output,
                        **{
                            'Epoch': deepcopy(data_dict_EFI['Epoch']),
                            'vxB_E': [np.array(vxB_ENU[:, 0]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS': 'vxB_E'}],
                            'vxB_N': [np.array(vxB_ENU[:, 1]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS': 'vxB_N'}],
                            'vxB_Up': [np.array(vxB_ENU[:, 2]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS': 'vxB_Up'}],
                            'vxB_E_RingCore': [np.array(vxB_ENU_RingCore[:, 0]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS': 'vxB_E_RingCore'}],
                            'vxB_N_RingCore': [np.array(vxB_ENU_RingCore[:, 1]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS': 'vxB_N_RingCore'}],
                            'vxB_Up_RingCore': [np.array(vxB_ENU_RingCore[:, 2]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS': 'vxB_Up_RingCore'}],
                            'vxB_X': [np.array(vxB_rkt[:, 0]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS':'vxB_X'}],
                            'vxB_Y': [np.array(vxB_rkt[:, 1]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS':'vxB_Y'}],
                            'vxB_Z': [np.array(vxB_rkt[:, 2]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS':'vxB_Z'}],
                            '|vxB|_rkt': [np.array(vxB_mag_rkt), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS':'|vxB|'}],
                            '|vxB|_ENU': [np.array(vxB_mag_ENU), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS': '|vxB|'}],
                            'E_X': deepcopy(data_dict_EFI['E_x']),
                            'E_Y': deepcopy(data_dict_EFI['E_y']),
                            'E_Z': deepcopy(data_dict_EFI['E_z']),
                            '|E|': [np.array(E_mag), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m','LABLAXIS':'|E|'}],
                            'DCM': [DCM_rkt_to_ENU, {'LABLAXIS':'DCM_rkt_to_ENU'}],
                            '|B|': [np.array(B_mag), {'DEPEND_0': 'Epoch', 'UNITS': 'T','LABLAXIS':'|B|'}],
                        }
                        }

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

