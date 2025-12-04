# --- ECEF_to_auroral_coordiantes_transformation_matrix.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: input the FAC transformation matrix and apply the
# east-north rotation angle indicated Plot_allsky_auroral_Coordinates_check

# OUTPUT:
# [X_ECEF, Y_ECEF, Z_ECEF] *transform_matrix = [N, T, p]

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

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = [4, 5]
outputData = True
outputPath_modifier = 'coordinates/transforms' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---

def FAC_to_auroral_coordiantes(wRocket):

    rocketID = ACESII.payload_IDs[wRocket-4]

    # --- ACES II Flight/Integration Data ---

    # --- get the data from the attitude file ---
    stl.prgMsg(f'Loading data')
    data_dict_transform = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}/coordinates/transforms/{ACESII.fliers[wRocket-4]}/*ECEF_to_FAC*')[0])
    data_dict_auroral_angle = stl.loadDictFromFile(glob(rf'{DataPaths.ACES_data_folder}/coordinates/auroral_coordinates/low/*.cdf*')[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
        'a11': [np.zeros(len(data_dict_transform['Epoch'][0])), {}],
        'a12': [np.zeros(len(data_dict_transform['Epoch'][0])), {}],
        'a13': [np.zeros(len(data_dict_transform['Epoch'][0])), {}],
        'a21': [np.zeros(len(data_dict_transform['Epoch'][0])), {}],
        'a22': [np.zeros(len(data_dict_transform['Epoch'][0])), {}],
        'a23': [np.zeros(len(data_dict_transform['Epoch'][0])), {}],
        'a31': [np.zeros(len(data_dict_transform['Epoch'][0])), {}],
        'a32': [np.zeros(len(data_dict_transform['Epoch'][0])), {}],
        'a33': [np.zeros(len(data_dict_transform['Epoch'][0])), {}],
        'ECEF_to_auroral': [np.zeros(shape=(len(data_dict_transform['Epoch'][0]), 3, 3)), {}],
        'auroral_to_ECEF': [np.zeros(shape=(len(data_dict_transform['Epoch'][0]), 3, 3)), {}],
        'Epoch': deepcopy(data_dict_transform['Epoch']),
        'Alt': deepcopy(data_dict_transform['Alt']),
        'Lat': deepcopy(data_dict_transform['Lat']),
        'Long': deepcopy(data_dict_transform['Long'])
    }

    ###############################
    # --- ROTATE THE FAC MATRIX ---
    ###############################

    FAC_to_auroral_transformation = np.array([np.matmul(stl.Rz(data_dict_auroral_angle['rotation_Angle'][0]),data_dict_transform['ECEF_to_FAC'][0][i]) for i in range(len(data_dict_transform['Epoch'][0]))])

    # store the outputs
    data_dict_output['a11'][0] = FAC_to_auroral_transformation[:, 0, 0]
    data_dict_output['a12'][0] = FAC_to_auroral_transformation[:, 0, 1]
    data_dict_output['a13'][0] = FAC_to_auroral_transformation[:, 0, 2]
    data_dict_output['a21'][0] = FAC_to_auroral_transformation[:, 1, 0]
    data_dict_output['a22'][0] = FAC_to_auroral_transformation[:, 1, 1]
    data_dict_output['a23'][0] = FAC_to_auroral_transformation[:, 1, 2]
    data_dict_output['a31'][0] = FAC_to_auroral_transformation[:, 2, 0]
    data_dict_output['a32'][0] = FAC_to_auroral_transformation[:, 2, 1]
    data_dict_output['a33'][0] = FAC_to_auroral_transformation[:, 2, 2]
    data_dict_output['ECEF_to_auroral'][0] = FAC_to_auroral_transformation
    data_dict_output['auroral_to_ECEF'][0] = np.array([FAC_to_auroral_transformation[i].T for i in range(len(data_dict_output['Epoch'][0]))])

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')

        exampleAttrs = {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': ACESII.epoch_fillVal,
                        'FORMAT': 'E12.2', 'UNITS': None, 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data',
                        'SCALETYP': 'linear',
                        'LABLAXIS': None}

        # update the data dict attrs
        for key, val in data_dict_output.items():
            newAttrs = deepcopy(exampleAttrs)

            for subKey, subVal in data_dict_output[key][1].items():
                newAttrs[subKey] = subVal

            data_dict_output[key][1] = newAttrs

        fileoutName = f'ACESII_{rocketID}_ECEF_to_auroral.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}{outputPath_modifier}/{ACESII.fliers[wRocket-4]}/{fileoutName}'
        stl.outputDataDict(outputPath, data_dict_output)
        stl.Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
for val in wRocket:
    FAC_to_auroral_coordiantes(val)