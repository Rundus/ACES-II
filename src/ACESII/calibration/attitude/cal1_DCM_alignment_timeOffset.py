# --- cal1_DCM_alignment_timeOffset.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: The EFI data revealed a time-offset in the attitude data DCM
# matrix. This effectively means our DCM does NOT perfectly align our rocket
# coordinates to ENU. In this program we try to estimate the "true"
# ENU coordiantes using the trajectory (lat,long,alt) information
# and compare it to the DCM to see if angles exist.

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

from src.ACESII.my_imports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
wRocket = 4
time_offset = -0.058 # in seconds. Determined from test1
time_slope = -0.000074 # in seconds. Determined from test1

# ---------------------------
outputData = True
# ---------------------------

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
def cal1_ENU_alignment_timeOffset(wRocket):

    # --- --- --- --- --- -
    # --- LOAD THE DATA ---
    # --- --- --- --- --- -

    # --- get the data from the file ---
    stl.prgMsg(f'Loading data')
    data_dict_attitude = stl.loadDictFromFile(glob(rf'C:\Data\ACESII\attitude\{ACESII.fliers[wRocket-4]}\*.cdf*')[0])
    stl.Done(start_time)

    # Prepare some variables
    T0 = dt.datetime(2022, 11, 20, 17, 20)
    T0_attitude = stl.EpochTo_T0_Rocket(data_dict_attitude['Epoch'][0], T0=T0)
    T0_adjust = T0_attitude*(1+time_slope) + time_offset

    # --- Apply the calibration time offset for the DCM ---
    for key in data_dict_attitude.keys():
        if key.lower() not in ['epoch']:
            data_dict_attitude[key][0] = np.interp(T0_adjust, T0_attitude, deepcopy(data_dict_attitude[key][0]))


    if outputData:
        stl.prgMsg('Creating output file')
        fileoutName = f'ACESII_{ACESII.payload_IDs[wRocket-4]}_Attitude_Solution_fullCal.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\\attitude\\{ACESII.fliers[wRocket-4]}\\'+ fileoutName
        stl.outputCDFdata(outputPath, data_dict_attitude)
        stl.Done(start_time)






# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

cal1_ENU_alignment_timeOffset(wRocket)
