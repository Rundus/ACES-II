# --- L2_to_L2_EFI_ENU_auroral_coordinates.py ---
# Description: convert EFI ENU coordinates to auroral and output as L2 data



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
remove_AC_flucuations = False
order= 4
freq_cutoff = 0.35
Plot_removed_AC = False

detrend_data = True
Plot_detrend = True
outputData = True


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none

def L2_to_L2_EFI_auroral_to_ENU_coordinates(wRocket):

    # --- get the data from the B-Field file ---
    stl.prgMsg(f'Loading data')
    data_dict_EFI = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket - 4]}\\*E_Field_auroral_fullCal.cdf*')[0])
    data_dict_transform_ENU = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\transforms\\{ACESII.fliers[wRocket - 4]}\\*ECEF_to_ENU.cdf*')[0])
    data_dict_transform_auroral = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\transforms\\{ACESII.fliers[wRocket - 4]}\\*ECEF_to_auroral.cdf*')[0])
    data_dict_LShell = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\Lshell\\{ACESII.fliers[wRocket - 4]}\\*Lshell.cdf*')[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
        'E_East' : [np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))),data_dict_EFI['E_N'][1]],
        'E_North': [np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))),data_dict_EFI['E_T'][1]],
        'E_Up': [np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))),data_dict_EFI['E_p'][1]],
        'E_mag':deepcopy(data_dict_EFI['Emag']),
        'Epoch':deepcopy(data_dict_EFI['Epoch'])
    }

    ###########################################
    # --- Interpolate trasnformation matrix ---
    ###########################################
    Epoch_EFI_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_EFI['Epoch'][0]])
    Epoch_ENU_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_transform_ENU['Epoch'][0]])
    Epoch_auroral_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_transform_auroral['Epoch'][0]])

    # interpolate auroral_to_ECEF matrix onto EFI timebase
    interp_keys = ['a11', 'a12', 'a13','a21','a22','a23','a31','a32','a33']
    for key in interp_keys:
        cs = CubicSpline(Epoch_ENU_tt2000, data_dict_transform_ENU[key][0])
        data_dict_transform_ENU[key][0] = deepcopy(cs(Epoch_EFI_tt2000))

    ECEF_to_ENU_matrix = np.array([
        [[data_dict_transform_ENU['a11'][0][i], data_dict_transform_ENU['a12'][0][i], data_dict_transform_ENU['a13'][0][i]],
        [data_dict_transform_ENU['a21'][0][i], data_dict_transform_ENU['a22'][0][i], data_dict_transform_ENU['a23'][0][i]],
        [data_dict_transform_ENU['a31'][0][i], data_dict_transform_ENU['a32'][0][i], data_dict_transform_ENU['a33'][0][i]]]
        for i in range(len(data_dict_EFI['Epoch'][0]))
    ])

    # interpolate ECEF_to_auroral matrix onto EFI timebase
    interp_keys = ['a11', 'a12', 'a13', 'a21', 'a22', 'a23', 'a31', 'a32', 'a33']
    for key in interp_keys:
        cs = CubicSpline(Epoch_auroral_tt2000, data_dict_transform_auroral[key][0])
        data_dict_transform_auroral[key][0] = deepcopy(cs(Epoch_EFI_tt2000))

    ECEF_to_auroral_matrix = np.array([
        [[data_dict_transform_auroral['a11'][0][i], data_dict_transform_auroral['a12'][0][i], data_dict_transform_auroral['a13'][0][i]],
         [data_dict_transform_auroral['a21'][0][i], data_dict_transform_auroral['a22'][0][i], data_dict_transform_auroral['a23'][0][i]],
         [data_dict_transform_auroral['a31'][0][i], data_dict_transform_auroral['a32'][0][i], data_dict_transform_auroral['a33'][0][i]]]
        for i in range(len(data_dict_EFI['Epoch'][0]))
    ])

    ###################################
    # --- Transform the Coordinates ---
    ###################################

    # form the EFI ENU vector
    EFI_auroral = np.array([data_dict_EFI['E_N'][0],data_dict_EFI['E_T'][0],data_dict_EFI['E_p'][0]]).T
    EFI_ECEF = np.array([np.matmul(ECEF_to_auroral_matrix[i].T, vec) for i, vec in enumerate(EFI_auroral)])
    EFI_ENU = np.array([np.matmul(ECEF_to_ENU_matrix[i], vec) for i, vec in enumerate(EFI_ECEF)])

    data_dict_output['E_East'][0] = EFI_ENU[:, 0]
    data_dict_output['E_East'][1]['LABLAXIS'] = 'E_East'

    data_dict_output['E_North'][0] = EFI_ENU[:, 1]
    data_dict_output['E_North'][1]['LABLAXIS'] = 'E_North'

    data_dict_output['E_Up'][0] = EFI_ENU[:, 2]
    data_dict_output['E_Up'][1]['LABLAXIS'] = 'E_Up'

    vec = np.array([data_dict_output['E_East'][0], data_dict_output['E_North'][0], data_dict_output['E_Up'][0]]).T
    data_dict_output = {**data_dict_output,
                        **{'|E|': [np.array([np.linalg.norm([v]) for v in vec]), deepcopy(data_dict_EFI['E_mag'][1])]}}

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')

        fileoutName = rf'ACESII_{ACESII.payload_IDs[wRocket-4]}_l2_E_Field_ENU_fullCal.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output, instrNam='EFI')
        stl.Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
L2_to_L2_EFI_auroral_to_ENU_coordinates(wRocket)