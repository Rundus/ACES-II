# --- cal4_apply_EFI_despin.py ---
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
time_offset = -0.058 # in seconds. Determined from cal3/test2
time_slope = -0.000074 # in seconds. Determined from cal3/test2
E_scale = 1 / np.sqrt(2)

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.my_imports import *
from scipy.interpolate import CubicSpline

#######################
# --- MAIN FUNCTION ---
#######################
def cal4_apply_EFI_despin(wRocket):

    # --- FILE I/O ---
    stl.prgMsg('Loading data')

    # --- get the EFI ---
    data_folder_path = rf'{DataPaths.ACES_data_folder}\\calibration\\EFI\\\cal1_spin_vxB_into_rktFrm_calibration\\{ACESII.fliers[wRocket - 4]}\\'
    data_dict_EFI = stl.loadDictFromFile(glob(data_folder_path+'*.cdf*')[0])

    # get the attitude data
    data_dict_attitude = stl.loadDictFromFile(glob(rf"C:\Data\ACESII\\attitude\\{ACESII.fliers[wRocket - 4]}\\*.cdf*")[0])

    # get the L-Shell data
    data_dict_LShell = st

    # --- prepare the output ---
    data_dict_output = {}
    stl.Done(start_time)

    #############################
    # --- DESPIN THE EFI DATA ---
    #############################
    stl.prgMsg('De-spinning EFI data')
    # [1] Interpolate the DCM onto EFI dataset with the time corretions determiend in cal2/test2
    DCM = np.zeros(shape=(len(data_dict_EFI['Epoch'][0]), 3, 3))
    T0 = dt.datetime(2022, 11, 20, 17, 20)
    T0_attitude = stl.EpochTo_T0_Rocket(data_dict_attitude['Epoch'][0], T0=T0)
    T0_attitude_adjust = T0_attitude*(1+time_slope) + time_offset
    T0_EFI = stl.EpochTo_T0_Rocket(data_dict_EFI['Epoch'][0], T0=T0)

    for i in range(3):
        for j in range(3):
            cs = CubicSpline(T0_attitude_adjust, data_dict_attitude[f'a{i + 1}{j + 1}'][0])
            DCM[:, i, j] = cs(T0_EFI)

    # [2] Apply the DCM
    E_rkt = np.array([data_dict_EFI['E_X'][0],
                      data_dict_EFI['E_Y'][0],
                      data_dict_EFI['E_Z'][0]]).T

    vxB_ENU = np.array([data_dict_EFI['vxB_E'][0],
                        data_dict_EFI['vxB_N'][0],
                        data_dict_EFI['vxB_Up'][0]
                        ]).T

    E_ENU = (E_scale) * np.array([np.matmul(DCM[i], E_rkt[i]) for i in range(len(DCM))]) - vxB_ENU
    Emag = np.array([np.linalg.norm(E_ENU[i]) for i in range(len(E_ENU))])

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
        fileoutName = f'ACESII_{ACESII.payload_IDs[wRocket-4]}_L2_EFI_ENU_fullCal.cdf'
        outputPath = f'C:\Data\ACESII\L2\\{ACESII.fliers[wRocket-4]}\\'+ fileoutName
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
cal4_apply_EFI_despin(wRocket)

