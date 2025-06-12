# --- L2_to_L2_EFI_auroral_to_FAC_coordinates.py ---
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

def L2_to_L2_EFI_auroral_to_FAC_coordinates(wRocket):

    # --- get the data from the B-Field file ---
    stl.prgMsg(f'Loading data')
    data_dict_EFI = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket - 4]}\\*E_Field_auroral_fullCal.cdf*')[0])
    data_dict_transform_auroral = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\transforms\\{ACESII.fliers[wRocket - 4]}\\*ECEF_to_auroral.cdf*')[0])
    data_dict_transform_FAC = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\transforms\\{ACESII.fliers[wRocket - 4]}\\*ECEF_to_FAC.cdf*')[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
        'E_r' : [np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))),data_dict_EFI['E_T'][1]],
        'E_e': [np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))),data_dict_EFI['E_N'][1]],
        'E_p': [np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))),data_dict_EFI['E_p'][1]],
        'E_mag':data_dict_EFI['Emag'],
        'Epoch':data_dict_EFI['Epoch']
    }

    ###########################################
    # --- Interpolate trasnformation matrix ---
    ###########################################

    Epoch_EFI_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_EFI['Epoch'][0]])
    Epoch_ENU_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_transform_auroral['Epoch'][0]])
    Epoch_FAC_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_transform_FAC['Epoch'][0]])

    # interpolate ENU_to_ECEF matrix onto EFI timebase
    interp_keys = ['a11','a12','a13','a21','a22','a23','a31','a32','a33']
    for key in interp_keys:
        cs = CubicSpline(Epoch_ENU_tt2000, data_dict_transform_auroral[key][0])
        data_dict_transform_auroral[key][0] = deepcopy(cs(Epoch_EFI_tt2000))

    ECEF_to_ENU_matrix = np.array([
        [[data_dict_transform_auroral['a11'][0][i], data_dict_transform_auroral['a12'][0][i], data_dict_transform_auroral['a13'][0][i]],
        [data_dict_transform_auroral['a21'][0][i], data_dict_transform_auroral['a22'][0][i], data_dict_transform_auroral['a23'][0][i]],
        [data_dict_transform_auroral['a31'][0][i], data_dict_transform_auroral['a32'][0][i], data_dict_transform_auroral['a33'][0][i]]]
        for i in range(len(data_dict_EFI['Epoch'][0]))
    ])

    # interpolate auroral_to_ECEF matrix onto EFI timebase
    interp_keys = ['a11', 'a12', 'a13', 'a21', 'a22', 'a23', 'a31', 'a32', 'a33']
    for key in interp_keys:
        cs = CubicSpline(Epoch_FAC_tt2000, data_dict_transform_FAC[key][0])
        data_dict_transform_FAC[key][0] = deepcopy(cs(Epoch_EFI_tt2000))

    ECEF_to_FAC_matrix = np.array([
        [[data_dict_transform_FAC['a11'][0][i], data_dict_transform_FAC['a12'][0][i], data_dict_transform_FAC['a13'][0][i]],
         [data_dict_transform_FAC['a21'][0][i], data_dict_transform_FAC['a22'][0][i], data_dict_transform_FAC['a23'][0][i]],
         [data_dict_transform_FAC['a31'][0][i], data_dict_transform_FAC['a32'][0][i], data_dict_transform_FAC['a33'][0][i]]]
        for i in range(len(data_dict_EFI['Epoch'][0]))
    ])

    ###################################
    # --- Transform the Coordinates ---
    ###################################

    # form the EFI ENU vector
    EFI_auroral = np.array([data_dict_EFI['E_N'][0],data_dict_EFI['E_T'][0],data_dict_EFI['E_p'][0]]).T
    EFI_ECEF = np.array([np.matmul(ECEF_to_ENU_matrix[i].T, vec) for i, vec in enumerate(EFI_auroral)])
    EFI_FAC = np.array([np.matmul(ECEF_to_FAC_matrix[i], vec) for i, vec in enumerate(EFI_ECEF)])

    data_dict_output['E_r'][0] = EFI_FAC[:, 0]
    data_dict_output['E_r'][1]['LABLAXIS'] = 'North-like Component'

    data_dict_output['E_e'][0] = EFI_FAC[:, 1]
    data_dict_output['E_r'][1]['LABLAXIS'] = 'East-like Component'

    data_dict_output['E_p'][0] = EFI_FAC[:, 2]
    data_dict_output['E_r'][1]['LABLAXIS'] = 'Field-Aligned Component'


    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')

        fileoutName = rf'ACESII_{ACESII.payload_IDs[wRocket-4]}_l2_E_Field_FAC_fullCal.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output, instrNam='EFI')
        stl.Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
L2_to_L2_EFI_auroral_to_FAC_coordinates(wRocket)