# --- L2_to_L2_RingCore_ECEF_to_auroral_coordinates.py ---
# Description: convert RingCore ENU coordinates to FAC and output as L2 data



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
from scipy.interpolate import CubicSpline

# --- Select the Rocket ---
# 4 -> ACES-II High Flier
# 5 -> ACES-II Low Flier
wRocket = 4
wFiles = [0]
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none

def L2_to_L2_RingCore_ECEF_to_auroal_coordinates(wRocket):

    # --- get the data from the B-Field file ---
    stl.prgMsg(f'Loading data')
    data_dict_mag = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}/L2/{ACESII.fliers[wRocket - 4]}/*RingCore_ECEF.cdf*')[0])
    data_dict_transform_FAC = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}/coordinates/transforms/{ACESII.fliers[wRocket - 4]}/*ECEF_to_auroral.cdf*')[0])
    data_dict_transform = deepcopy(data_dict_transform_FAC)
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
        'B_N': [np.zeros(shape=(len(data_dict_mag['Epoch'][0]))), data_dict_mag['B_X'][1]],
        'B_T': [np.zeros(shape=(len(data_dict_mag['Epoch'][0]))), data_dict_mag['B_Y'][1]],
        'B_p': [np.zeros(shape=(len(data_dict_mag['Epoch'][0]))), data_dict_mag['B_Z'][1]],
        'B_model_N': [np.zeros(shape=(len(data_dict_mag['Epoch'][0]))), data_dict_mag['B_model_X'][1]],
        'B_model_T': [np.zeros(shape=(len(data_dict_mag['Epoch'][0]))), data_dict_mag['B_model_Y'][1]],
        'B_model_p': [np.zeros(shape=(len(data_dict_mag['Epoch'][0]))), data_dict_mag['B_model_Z'][1]],
        'Bmag': data_dict_mag['Bmag'],
        'Epoch': data_dict_mag['Epoch']
    }

    ###########################################
    # --- Interpolate transformation matrix ---
    ###########################################
    Epoch_data_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_mag['Epoch'][0]])
    Epoch_transform_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_transform['Epoch'][0]])

    # interpolate transformation matrix onto data timebase
    interp_keys = ['a11','a12','a13','a21','a22','a23','a31','a32','a33']
    for key in interp_keys:
        cs = CubicSpline(Epoch_transform_tt2000, data_dict_transform[key][0])
        data_dict_transform[key][0] = deepcopy(cs(Epoch_data_tt2000))

    transform_matrix_interp = np.array([
        [[data_dict_transform['a11'][0][i], data_dict_transform['a12'][0][i], data_dict_transform['a13'][0][i]],
        [data_dict_transform['a21'][0][i], data_dict_transform['a22'][0][i], data_dict_transform['a23'][0][i]],
        [data_dict_transform['a31'][0][i], data_dict_transform['a32'][0][i], data_dict_transform['a33'][0][i]]]
        for i in range(len(data_dict_mag['Epoch'][0]))
    ])


    ###################################
    # --- Transform the Coordinates ---
    ###################################

    # form the mag vector
    mag_vec = np.array([data_dict_mag['B_X'][0], data_dict_mag['B_Y'][0], data_dict_mag['B_Z'][0]]).T
    mag_model_vec = np.array([data_dict_mag['B_model_X'][0], data_dict_mag['B_model_Y'][0], data_dict_mag['B_model_Z'][0]]).T

    mag_transformed = np.array([np.matmul(transform_matrix_interp[i], vec) for i, vec in enumerate(mag_vec)])
    mag_model_transformed = np.array([np.matmul(transform_matrix_interp[i], vec) for i, vec in enumerate(mag_model_vec)])

    data_dict_output['B_N'][0] = mag_transformed[:, 0]
    data_dict_output['B_N'][1]['LABLAXIS'] = 'Arc-Normal Component'
    data_dict_output['B_model_N'][0] = mag_model_transformed[:, 0]
    data_dict_output['B_model_N'][1]['LABLAXIS'] = 'Arc-Normal Model Component'

    data_dict_output['B_T'][0] = mag_transformed[:, 1]
    data_dict_output['B_T'][1]['LABLAXIS'] = 'Arc-Tangent Component'
    data_dict_output['B_model_T'][0] = mag_model_transformed[:, 1]
    data_dict_output['B_model_T'][1]['LABLAXIS'] = 'Arc-Tangent Model Component'

    data_dict_output['B_p'][0] = mag_transformed[:, 2]
    data_dict_output['B_p'][1]['LABLAXIS'] = 'Field-Aligned Component'
    data_dict_output['B_model_p'][0] = mag_model_transformed[:, 2]
    data_dict_output['B_model_p'][1]['LABLAXIS'] = 'Field-Aligned Model Component'

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')
        fileoutName = rf'ACESII_{ACESII.payload_IDs[wRocket-4]}_l2_RingCore_auroral.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/L2/{ACESII.fliers[wRocket-4]}/{fileoutName}'
        stl.outputDataDict(outputPath, data_dict_output, instrNam='RingCore')
        stl.Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
L2_to_L2_RingCore_ECEF_to_auroal_coordinates(wRocket)