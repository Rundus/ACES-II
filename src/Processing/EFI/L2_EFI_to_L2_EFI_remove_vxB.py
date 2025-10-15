# --- L2_EFI_to_L2_EFI_remove_vxB.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Load in the De-spun EFI data, use the calibrated vxB values
# which were used to determine the time-offset then determine
# E = E' - vxB. Also re-scale the E-Field with whatever you want

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
E_scale = 1/np.sqrt(2)

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.my_imports import *

#######################
# --- MAIN FUNCTION ---
#######################
def L2_EFI_to_L2_EFI_remove_vxB(wRocket):

    # --- FILE I/O ---

    # get the EFI files
    stl.prgMsg('Loading EFI Data')
    data_folder_path = rf'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket - 4]}\\'
    input_files = glob(data_folder_path + '*_l2_EFI_ENU_withVxB.cdf*')
    data_dict_EFI = stl.loadDictFromFile(input_files[0])
    stl.Done(start_time)

    # Get the vxB data
    stl.prgMsg('Loading vxB Data')
    data_folder_path = rf'{DataPaths.ACES_data_folder}\\calibration\\\EFI_cal1_spin_vxB_into_rktFrm_calibration\\{ACESII.fliers[wRocket-4]}\\'
    input_files_cal = glob(data_folder_path+'*.cdf*')[0]
    data_dict_vxB = stl.loadDictFromFile(input_files_cal)
    stl.Done(start_time)

    #################################
    # --- SUBTRACT THE vxB effect ---
    #################################
    stl.prgMsg('Correcting E=E-vxB')
    E_Field = np.array([data_dict_EFI['E_E'][0],data_dict_EFI['E_N'][0],data_dict_EFI['E_Up'][0]]).T
    # vxB = np.array([data_dict_vxB['vxB_E'][0],data_dict_vxB['vxB_N'][0],data_dict_vxB['vxB_Up'][0]]).T # Note, this vector is already -vxB is is NOT |vxB|
    vxB = np.array([data_dict_EFI['vxB_E'][0], data_dict_EFI['vxB_N'][0], data_dict_EFI['vxB_Up'][0]]).T  # Note, this vector is already -vxB is is NOT |vxB|

    E_Field_corrected = E_scale*E_Field - vxB
    E_mag = np.array([np.linalg.norm(vec) for vec in E_Field_corrected])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}

    # store everything
    data_dict_output = {**data_dict_output,
                        **{
                            'E_E': [E_Field_corrected[:, 0], {'DEPEND_0':'Epoch','UNITS': 'V/m','LABLAXIS':'E_East'}],
                            'E_N': [E_Field_corrected[:, 1], {'DEPEND_0':'Epoch','UNITS': 'V/m','LABLAXIS':'E_North'}],
                            'E_Up': [E_Field_corrected[:, 2], {'DEPEND_0':'Epoch','UNITS': 'V/m','LABLAXIS':'E_up'}],
                            '|E|': [E_mag, {'DEPEND_0': 'Epoch', 'UNITS':'V/m',  'LABLAXIS':'|E|'}],
                            'Epoch': deepcopy(data_dict_EFI['Epoch']),
                        }
                        }

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        stl.prgMsg('Creating output file')
        fileoutName = f'ACESII_{ACESII.payload_IDs[wRocket-4]}_l2_EFI_ENU_fullCal.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket-4]}\\'+ fileoutName
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
L2_EFI_to_L2_EFI_remove_vxB(wRocket)

