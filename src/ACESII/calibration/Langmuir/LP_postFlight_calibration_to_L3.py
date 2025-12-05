# --- LP_postFlight_calibration_to_L3.py ---
# Desciption: Determine the calibration function to apply to the ni_density
# and DERPA temperature profiles

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
from src.ACESII.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False  # Just print the names of files

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
from glob import glob
import spaceToolsLib as stl
from scipy.optimize import curve_fit
import scipy



#######################
# --- MAIN FUNCTION ---
#######################
def LP_postFlight_calibration(wRocket):

    # --- FILE I/O ---
    stl.prgMsg('Loading Data')

    # Load the calibration data
    data_path = glob(rf'{DataPaths.ACES_data_folder}/calibration/LP/postFlight_calibration/{ACESII.fliers[wRocket - 4]}/*_postFlight_cal*')[0]
    data_dict_calData = stl.loadDictFromFile(data_path)

    # load the rocket ion saturation current data
    data_path = glob(rf'{DataPaths.ACES_data_folder}/L3/Langmuir/{ACESII.fliers[wRocket-4]}/*l3_langmuir_fixed.cdf*')[0]
    data_dict_LP_current = stl.loadDictFromFile(data_path)

    # load the DERPA
    data_path = glob(rf'{DataPaths.ACES_data_folder}/L2/{ACESII.fliers[wRocket - 4]}/*_ERPA1.cdf*')[0]
    data_dict_DERPA1 = stl.loadDictFromFile(data_path)
    data_path = glob(rf'{DataPaths.ACES_data_folder}/L2/{ACESII.fliers[wRocket - 4]}/*_ERPA2.cdf*')[0]
    data_dict_DERPA2 = stl.loadDictFromFile(data_path)
    stl.Done(start_time)

    ################################################
    # --- INTERPOLATE L-SHELL INTO LP/DERPA DATA ---
    ################################################
    T0 = dt.datetime(2022,11,20,17,20)
    T0_LP = np.array(stl.EpochTo_T0_Rocket(data_dict_LP_current['Epoch'][0], T0=T0),dtype='float64')
    T0_DERPA1 = stl.EpochTo_T0_Rocket(data_dict_DERPA1['Epoch'][0], T0=T0)
    T0_DERPA2 = stl.EpochTo_T0_Rocket(data_dict_DERPA2['Epoch'][0], T0=T0)

    #############################################
    # --- DETERMINE THE CALIBRATION FUNCTIONS ---
    #############################################
    # Description: Isolate the regions of the data corresponding to "background"
    # and cross-calibrate between these

    def fitFunc(x,a):
        return a*x

    # --- ni density ---
    # reduce data to background
    target_times_bg = [
        [dt.datetime(2022,11,20,17,21,8),dt.datetime(2022,11,20,17,24,0)], # High Flyer
        [dt.datetime(2022,11,20,17,23,40),dt.datetime(2022,11,20,17,24,42)] # Low Flyer
    ]

    low_idx, high_idx = np.abs(data_dict_LP_current['Epoch'][0] - target_times_bg[wRocket-4][0]).argmin(), np.abs(data_dict_LP_current['Epoch'][0] - target_times_bg[wRocket-4][1]).argmin()
    ni_data = data_dict_LP_current['ni'][0][low_idx:high_idx+1]
    ni_EISCAT_data = data_dict_calData['ne'][0][low_idx:high_idx+1]
    good_idxs = np.where(np.isnan(ni_data)==False)
    ni_data = ni_data[good_idxs]
    ni_EISCAT_data = ni_EISCAT_data[good_idxs]
    params_LP, cov = curve_fit(fitFunc,ni_data,ni_EISCAT_data)

    # --- Te ---
    target_times_bg = [
        [dt.datetime(2022, 11, 20, 17, 21, 43), dt.datetime(2022, 11, 20, 17, 22, 16)],  # High Flyer
        [dt.datetime(2022, 11, 20, 17, 24, 4), dt.datetime(2022, 11, 20, 17, 24, 27)]  # Low Flyer
    ]

    # DERPA1
    low_idx, high_idx = np.abs(data_dict_calData['Epoch'][0] - target_times_bg[wRocket - 4][0]).argmin(), np.abs(data_dict_calData['Epoch'][0] - target_times_bg[wRocket - 4][1]).argmin()
    Te_data = data_dict_calData['Te_DERPA1'][0][low_idx:high_idx + 1]
    Te_EISCAT_data = data_dict_calData['Te'][0][low_idx:high_idx + 1]
    params_DERPA1, cov = curve_fit(fitFunc, Te_data, Te_EISCAT_data)

    # DERPA2
    low_idx, high_idx = np.abs(data_dict_calData['Epoch'][0] - target_times_bg[wRocket - 4][0]).argmin(), np.abs(data_dict_calData['Epoch'][0] - target_times_bg[wRocket - 4][1]).argmin()
    Te_data = data_dict_calData['Te_DERPA2'][0][low_idx:high_idx + 1]
    Te_EISCAT_data = data_dict_calData['Te'][0][low_idx:high_idx + 1]
    params_DERPA2, cov = curve_fit(fitFunc, Te_data, Te_EISCAT_data)

    print(params_LP)
    print(params_DERPA1)
    print(params_DERPA2)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:
        stl.prgMsg('Creating output file')

        # Write out the calibrated Langmuir Probe Data
        fileoutName = f'ACESII_{ACESII.payload_IDs[wRocket-4]}_l3_langmuir_fixed_fullCal.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/L3/Langmuir/{ACESII.fliers[wRocket-4]}/' + fileoutName
        data_dict_output = deepcopy(data_dict_LP_current)
        data_dict_output['ni'][0] = fitFunc(data_dict_output['ni'][0],*params_LP)

        # add the filtered_ni data via interpolation (to remove nans) then filter

        def fill_nan(A):
            '''
            interpolate to fill nan values
            '''
            inds = np.arange(A.shape[0])
            good = np.where(np.isfinite(A))
            f = scipy.interpolate.interp1d(inds[good], A[good], bounds_error=False)
            B = np.where(np.isfinite(A), A, f(inds))
            return B

        # remove nans (they will be filtered out)
        data_dict_output['ni'][0][np.where(np.isnan(data_dict_output['ni'][0]))] = 0

        ni_filtered = stl.butterFilter().butter_filter(data=deepcopy(data_dict_output['ni'][0]),
                                                       highcutoff=0.3,
                                                       lowcutoff=0,
                                                       fs=1000,
                                                       order=4,
                                                       filtertype='lowpass')
        data_dict_output ={**data_dict_output,
                           **{'ni_filtered':[ni_filtered,deepcopy(data_dict_output['ni'][1])]}
                           }

        stl.outputDataDict(outputPath, data_dict_output)

        # Write out the calibrated DERPA1 Data
        fileoutName = f'ACESII_{ACESII.payload_IDs[wRocket - 4]}_l3_ERPA1_fullCal.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/L3/DERPA/{ACESII.fliers[wRocket - 4]}/' + fileoutName
        data_dict_output = deepcopy(data_dict_DERPA1)
        data_dict_output['temperature'][0] = fitFunc(data_dict_output['temperature'][0], *params_DERPA1)
        stl.outputDataDict(outputPath, data_dict_output)

        # Write out the calibrated DERPA1 Data
        fileoutName = f'ACESII_{ACESII.payload_IDs[wRocket - 4]}_l3_ERPA2_fullCal.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/L3/DERPA/{ACESII.fliers[wRocket - 4]}/' + fileoutName
        data_dict_output = deepcopy(data_dict_DERPA2)
        data_dict_output['temperature'][0] = fitFunc(data_dict_output['temperature'][0], *params_DERPA2)
        stl.outputDataDict(outputPath, data_dict_output)


        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
LP_postFlight_calibration(wRocket)

