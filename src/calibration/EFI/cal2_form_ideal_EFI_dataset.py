# --- cal2_form_ideal_EFI_dataset ---
# --- Author: C. Feltman ---
# DESCRIPTION: After we have determined a good vxB from EFI cal1, we need to
# determining the timing offset for the DCM to correct the ENU misalignment.
# To do this, we need to determine an "ideal" dataset to calibrate to when
# applying timing offsets. Here we assume that virtually ALL the DC offset in our
# EFI data can be explained by vxB, so we form a "perfect" EFI dataset by
# fudging our calculations. We later determine the rotation required to
# make that perfect dataset.

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"

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
from scipy.interpolate import CubicSpline
import spaceToolsLib as stl

#######################
# --- MAIN FUNCTION ---
#######################
def cal2_form_ideal_EFI_dataset(wRocket):

    # --- FILE I/O ---

    # get the EFI files with vxB included
    data_folder_path = rf'C:\Data\ACESII\calibration\EFI\cal1_spin_vxB_into_rktFrm_calibration\\{ACESII.fliers[wRocket - 4]}\\'
    input_files = glob(data_folder_path + '*vxB_rktFrm*')
    data_dict_EFI = stl.loadDictFromFile(input_files[0])

    # get the attitude data
    data_dict_attitude = stl.loadDictFromFile(glob(rf"C:\Data\ACESII\\attitude\\{ACESII.fliers[wRocket-4]}\\*.cdf*")[0])

    # --- prepare the output ---
    data_dict_output = {}

    ######################################
    # --- FORM THE PERFECT EFI DATASET ---
    ######################################
    stl.prgMsg('Forming Dataset')

    # [1] Interpolate the DCM onto EFI dataset as-is
    DCM = np.zeros(shape=(len(data_dict_EFI['Epoch'][0]), 3, 3))
    T0 = dt.datetime(2022, 11, 20, 17, 20)
    T0_attitude = stl.EpochTo_T0_Rocket(data_dict_attitude['Epoch'][0], T0=T0)
    T0_EFI = stl.EpochTo_T0_Rocket(data_dict_EFI['Epoch'][0], T0=T0)
    for i in range(3):
        for j in range(3):
            cs = CubicSpline(T0_attitude,data_dict_attitude[f'a{i+1}{j+1}'][0])
            DCM[:, i, j] = cs(T0_EFI)

    # [2] Apply the DCM
    E_rkt = np.array([data_dict_EFI['E_X'][0],
                      data_dict_EFI['E_Y'][0],
                      data_dict_EFI['E_Z'][0]]).T

    E_scale = 1/np.sqrt(2)
    E_ENU = (E_scale)*np.array([np.matmul(DCM[i],E_rkt[i]) for i in range(len(DCM))])

    # [3] Fudge the rotated output to be what you want
    E_ENU[:, 0] = E_ENU[:, 0] - data_dict_EFI['vxB_E'][0]
    E_ENU[:, 1] = E_ENU[:, 1] + data_dict_EFI['vxB_N'][0]
    E_ENU[:, 2] = E_ENU[:, 2] + data_dict_EFI['vxB_Up'][0]
    Emag = np.array([np.linalg.norm(E_ENU[i]) for i in range(len(DCM))])

    # store everything
    data_dict_output = {**data_dict_output,
                        **{
                            'E_E': [E_ENU[:, 0], {'DEPEND_0':'Epoch','UNITS': 'V/m','LABLAXIS':'E_East'}],
                            'E_N': [E_ENU[:, 1], {'DEPEND_0':'Epoch','UNITS': 'V/m','LABLAXIS':'E_North'}],
                            'E_Up': [E_ENU[:, 2], {'DEPEND_0':'Epoch','UNITS': 'V/m','LABLAXIS':'E_Up'}],
                            '|E|': [Emag, {'DEPEND_0': 'Epoch', 'UNITS':'V/m',  'LABLAXIS':'|E|'}],
                            'Epoch': deepcopy(data_dict_EFI['Epoch']),
                        }
                        }
    stl.Done(start_time)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        stl.prgMsg('Creating output file')
        fileoutName = f'ACESII_{ACESII.payload_IDs[wRocket-4]}_cal2_form_ideal_EFI_dataset.cdf'
        outputPath = f'C:\Data\ACESII\calibration\EFI\cal2_form_ideal_EFI_dataset\low\\'+ fileoutName
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
cal2_form_ideal_EFI_dataset(wRocket)

