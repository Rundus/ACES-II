# --- cal3_EFI_despin.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"

import numpy as np

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

#######################
# --- MAIN FUNCTION ---
#######################
def EFI_rkt_to_ENU_despin(wRocket):

    # --- FILE I/O ---

    # get the EFI files
    data_folder_path = rf'{DataPaths.ACES_data_folder}\\L1\\{ACESII.fliers[wRocket - 4]}\\'
    input_files = glob(data_folder_path + '*EFI_rktFrm*')

    # Get the fitted timing offset data
    data_folder_path = rf'{DataPaths.ACES_data_folder}\\calibration\\EFI_cal2_timing_offset_calibration\\{ACESII.fliers[wRocket-4]}\\'
    input_files_cal = glob(data_folder_path+'*.cdf*')[0]

    # --- get the data from the L2 file ---
    stl.prgMsg(f'Loading data from L1 EFI Files')
    data_dict_EFI = stl.loadDictFromFile(input_files[0])
    data_dict_cal = stl.loadDictFromFile(input_files_cal)
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}

    #############################
    # --- DESPIN THE EFI DATA ---
    #############################
    stl.prgMsg('De-spinning EFI data')

    # form the rocket-frame vectors
    EFI_rkt_vec = np.array([data_dict_EFI['E_x'][0], data_dict_EFI['E_y'][0], data_dict_EFI['E_z'][0]]).T
    vxB_rkt_vec = np.array([data_dict_cal['vxB_X'][0], data_dict_cal['vxB_Y'][0], data_dict_cal['vxB_Z'][0]]).T

    # Apply the DCM (rkt-->to_ENU)
    EFI_ENU_vec = np.array([np.matmul(data_dict_cal['DCM'][0][idx],vec) for idx,vec in enumerate(EFI_rkt_vec)])

    # TODO: Something is VERY wrong with the vxB below - Causing E_East to appear negative
    vxB_ENU_vec = np.array([np.matmul(data_dict_cal['DCM'][0][idx],vec) for idx,vec in enumerate(vxB_rkt_vec)])

    # Calculate some magnitudes for diagnostics
    E_mag_rkt = np.array([np.linalg.norm(vec) for vec in EFI_rkt_vec])
    E_mag_ENU = np.array([np.linalg.norm(vec) for vec in EFI_ENU_vec])

    # store everything
    data_dict_output = {**data_dict_output,
                        **{
                            'E_E': [EFI_ENU_vec[:, 0], {'DEPEND_0':'Epoch','UNITS': 'V/m','LABLAXIS':'E_East'}],
                            'E_N': [EFI_ENU_vec[:, 1], {'DEPEND_0':'Epoch','UNITS': 'V/m','LABLAXIS':'E_North'}],
                            'E_Up': [EFI_ENU_vec[:, 2], {'DEPEND_0':'Epoch','UNITS': 'V/m','LABLAXIS':'E_Up'}],
                            '|E|_ENU': [E_mag_ENU, {'DEPEND_0': 'Epoch', 'UNITS':'V/m',  'LABLAXIS':'|E|'}],
                            '|E|_rkt': [E_mag_rkt, {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS': '|E|'}],
                            'Epoch': deepcopy(data_dict_EFI['Epoch']),
                            'vxB_E': [np.array(vxB_ENU_vec[:, 0]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS': 'vxB_E'}],
                            'vxB_N': [np.array(vxB_ENU_vec[:, 1]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS': 'vxB_N'}],
                            'vxB_Up': [np.array(vxB_ENU_vec[:, 2]), {'DEPEND_0': 'Epoch', 'UNITS': 'V/m', 'LABLAXIS': 'vxB_Up'}],
                            '|vxB|': deepcopy(data_dict_cal['|vxB|'])
                        }
                        }
    stl.Done(start_time)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        stl.prgMsg('Creating output file')
        fileoutName = f'ACESII_{ACESII.payload_IDs[wRocket-4]}_EFI_ENU_withVxB.cdf'
        outputPath = f'C:\Data\ACESII\calibration\EFI_cal3_depsin_EFI\low\\'+ fileoutName
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
EFI_rkt_to_ENU_despin(wRocket)

