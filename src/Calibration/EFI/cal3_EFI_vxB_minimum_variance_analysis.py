# --- cal3_EFI_vxB_minimum_variace_analysis.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: minimize the variance between the
# EFI and vxB profiles through the ChiSquare value

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
E_field_scale = 1/(np.sqrt(2))

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

# --- OutputData ---
outputData = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.my_imports import *
import spaceToolsLib as stl
from tqdm import tqdm

#######################
# --- MAIN FUNCTION ---
#######################
def cal3_EFI_vxB_minimum_variace_analysis(wRocket):

    # --- FILE I/O ---
    data_folder_path = rf'{DataPaths.ACES_data_folder}calibration\\EFI_spin_vxB_into_rktFrm_calibration\\{ACESII.fliers[wRocket-4]}\\'

    # get the EFI files
    input_files = glob(data_folder_path + '*vxB_rktFrm.cdf*')

    # --- get the data from the L2 file ---
    stl.prgMsg(f'Loading data from L1 Files')
    data_dict_EFI = stl.loadDictFromFile(input_files[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}

    ##################################
    # --- Perform Minimum Variance ---
    ##################################

    # [1] Isolate only the data to performance variance on
    target_times = [
        dt.datetime(2022, 11, 20, 17, 23, 43),
        dt.datetime(2022, 11, 20, 17, 27, 51),
    ]
    low_idx, high_idx = np.abs(data_dict_EFI['Epoch'][0] - target_times[0]).argmin(),np.abs(data_dict_EFI['Epoch'][0] - target_times[1]).argmin()
    E = np.array([data_dict_EFI['E_X'][0][low_idx:high_idx+1], data_dict_EFI['E_Y'][0][low_idx:high_idx+1],data_dict_EFI['E_Z'][0][low_idx:high_idx+1]]).T
    vxB = np.array([data_dict_EFI['vxB_X'][0][low_idx:high_idx+1], data_dict_EFI['vxB_Y'][0][low_idx:high_idx+1],data_dict_EFI['vxB_Z'][0][low_idx:high_idx+1]]).T
    T0=dt.datetime(2022,11,20,17,20)
    Epoch_T0 = stl.EpochTo_T0_Rocket(data_dict_EFI['Epoch'][0][low_idx:high_idx+1],T0=T0)

    # [2] Define the minimum variance search space + outputs
    N = 10
    tolerance = 100/100
    slope_guess = -6.329794273708472E-4 # in vxB second/EFI second
    intercept_guess = 0.50112 # in seconds
    slope_space = np.linspace((1-tolerance)*slope_guess, (1+tolerance)*slope_guess, N)
    intercept_space = np.linspace((1-tolerance)*intercept_guess, (1+tolerance)*intercept_guess, N)
    ChiSquare = []
    slopes = []
    intercepts = []

    # [3] Loop over the search space to find the minimum ChiSquare Value
    def fitFunc(x, A, B):
        return A * x + B

    for idx1 in tqdm(range(len(slope_space))):
        for idx2 in range(len(intercept_space)):

            slp = slope_space[idx1]
            intp = intercept_space[idx2]
            params = [slp, intp]

            # Shift the timebase of vxB by the guess
            Epoch_vxB_new = Epoch_T0 - fitFunc(Epoch_T0, *params)
            vxB_Y_new = np.interp(Epoch_T0, Epoch_vxB_new, vxB[:,1])
            vxB_Z_new = np.interp(Epoch_T0, Epoch_vxB_new, vxB[:,2])

            vxB_new = np.array([vxB[:,0], vxB_Y_new, vxB_Z_new]).T

            # Calculate the ChiSquare value
            ChiSquare.append(np.sum(np.power(vxB_new - E, 2)))
            slopes.append(slp)
            intercepts.append(intp)

    best_idx = np.abs(np.array(ChiSquare) - np.min(ChiSquare)).argmin()
    print(ChiSquare[best_idx])
    print(slopes[best_idx])
    print(intercepts[best_idx])

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    # if outputData:
        # stl.prgMsg('Creating output file')
        # data_dict_output = {
        #                     'Epoch': deepcopy(data_dict_EFI['Epoch']),
        #                     'Epoch_EFI_T0': [Epoch_T0, deepcopy(data_dict_EFI['Epoch'][1])],
        #                     'Y_fitParams': [np.array([params_1]), {}],
        #                     'Z_fitParams': [np.array([params_1]), {}],
        #                     'deltaT_Y': [deltaT_epoch[1], {}],
        #                     'deltaT_Z': [deltaT_epoch[2], {}],
        #                     'Y_fit': [np.array(fitFunc(stl.EpochTo_T0_Rocket(deltaT_epoch[2], T0=T0), *params_1)), {'DEPEND_0':'deltaT_Y'}],
        #                     'Z_fit': [np.array(fitFunc(stl.EpochTo_T0_Rocket(deltaT_epoch[2], T0=T0), *params_2)), {'DEPEND_0':'deltaT_Z'}],
        #                     'vxB_X': deepcopy(data_dict_EFI['vxB_X']),
        #                     'vxB_Y': [vxB_Y_new, deepcopy(data_dict_EFI['vxB_Y'][1]), deepcopy(data_dict_EFI['vxB_Y'][1])],
        #                     'vxB_Z': [vxB_Z_new, deepcopy(data_dict_EFI['vxB_Z'][1]), deepcopy(data_dict_EFI['vxB_Z'][1])],
        #                     '|vxB|': [np.array([np.linalg.norm(vxB[i]) for i in range(len(vxB))]), deepcopy(data_dict_EFI['|vxB|'][1])]
        #                     }
        #
        #
        # fileoutName = f'ACESII_{ACESII.payload_IDs[wRocket-4]}_l1_EFI_vxB_rktFrm.cdf'
        # outputPath = f'{DataPaths.ACES_data_folder}\calibration\EFI_timing_offset_calibration\\{ACESII.fliers[wRocket-4]}\\'+ fileoutName
        # stl.outputCDFdata(outputPath, data_dict_output)
        # stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
cal3_EFI_vxB_minimum_variace_analysis(wRocket)

