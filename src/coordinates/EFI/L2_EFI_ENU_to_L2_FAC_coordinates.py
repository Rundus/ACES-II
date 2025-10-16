# --- L2_EFI_ENU_to_L2_FAC_coordinates.py ---
# Description: convert EFI ENU coordinates to FAC and output as L2 data



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
from src.my_imports import *
start_time = time.time()
# --- --- --- --- ---



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
from scipy.interpolate import CubicSpline

# --- Select the Rocket ---
# 4 -> ACES-II High Flier
# 5 -> ACES-II Low Flier
wRocket = 5
wFiles = [0]
outputData = True

Plot_correction_term = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none

def L2_to_L2_EFI_FAC_to_auroral_coordinates(wRocket):

    # --- get the data from the B-Field file ---
    stl.prgMsg(f'Loading data')
    data_dict_EFI = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket - 4]}\\*EFI_ENU_fullCal.cdf*')[0])
    data_dict_transform_ENU = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\transforms\\{ACESII.fliers[wRocket - 4]}\\*ECEF_to_ENU.cdf*')[0])
    data_dict_transform_FAC = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\transforms\\{ACESII.fliers[wRocket - 4]}\\*ECEF_to_FAC.cdf*')[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
        'E_r' : [np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))),data_dict_EFI['E_E'][1]],
        'E_e': [np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))),data_dict_EFI['E_N'][1]],
        'E_p': [np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))),data_dict_EFI['E_Up'][1]],
        'Epoch':data_dict_EFI['Epoch']
    }

    ###########################################
    # --- Interpolate trasnformation matrix ---
    ###########################################
    stl.prgMsg('Transforming Coordinates')

    Epoch_EFI_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_EFI['Epoch'][0]])
    Epoch_transform_ENU_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_transform_ENU['Epoch'][0]])
    Epoch_transform_FAC_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_transform_FAC['Epoch'][0]])

    # interpolate ENU_to_ECEF matrix onto EFI timebase
    ENU_to_ECEF_interp = np.zeros(shape=(len(data_dict_EFI['Epoch'][0]), 3, 3))
    for i in range(1, 4):
        for j in range(1, 4):
            cs = CubicSpline(Epoch_transform_ENU_tt2000, data_dict_transform_ENU['ENU_to_ECEF'][0][:, i - 1, j - 1])
            ENU_to_ECEF_interp[:, i - 1, j - 1] = cs(Epoch_EFI_tt2000)

    # interpolate ECEF_to_FAC matrix onto EFI timebase
    ECEF_to_FAC_interp = np.zeros(shape=(len(data_dict_EFI['Epoch'][0]), 3, 3))
    for i in range(1, 4):
        for j in range(1, 4):
            cs = CubicSpline(Epoch_transform_FAC_tt2000, data_dict_transform_FAC['ECEF_to_FAC'][0][:, i - 1, j - 1])
            ECEF_to_FAC_interp[:, i - 1, j - 1] = cs(Epoch_EFI_tt2000)

    ###################################
    # --- Transform the Coordinates ---
    ###################################

    # form the EFI ENU vector
    EFI_ENU_vec = np.array([data_dict_EFI['E_E'][0], data_dict_EFI['E_N'][0], data_dict_EFI['E_Up'][0]]).T
    EFI_ECEF_vec = np.array([np.matmul(ENU_to_ECEF_interp[i], vec) for i, vec in enumerate(EFI_ENU_vec)])
    EFI_FAC_vec = np.array([np.matmul(ECEF_to_FAC_interp[i], vec) for i, vec in enumerate(EFI_ECEF_vec)])

    data_dict_output['E_r'][0] = EFI_FAC_vec[:, 0]
    data_dict_output['E_r'][1]['LABLAXIS'] = 'North-like Component'

    data_dict_output['E_e'][0] = EFI_FAC_vec[:, 1]
    data_dict_output['E_e'][1]['LABLAXIS'] = 'East-like Component'

    data_dict_output['E_p'][0] = EFI_FAC_vec[:, 2]
    data_dict_output['E_p'][1]['LABLAXIS'] = 'Field-Aligned Component'


    vec = np.array([data_dict_output['E_r'][0],data_dict_output['E_e'][0],data_dict_output['E_p'][0]]).T
    data_dict_output = {**data_dict_output,
                        **{'|E|':[np.array([np.linalg.norm([v]) for v in vec]),deepcopy(data_dict_EFI['|E|'][1])]}}
    stl.Done(start_time)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')

        fileoutName = rf'ACESII_{ACESII.payload_IDs[wRocket-4]}_l2_EFI_FAC_fullCal.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output, instrNam='EFI')
        stl.Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
L2_to_L2_EFI_FAC_to_auroral_coordinates(wRocket)