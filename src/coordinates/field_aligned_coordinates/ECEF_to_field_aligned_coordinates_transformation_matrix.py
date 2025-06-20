# --- ECEF_to_field_aligned_coordinates_transformation_matrix.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: input the ENU_to_Field_aligned transformation matrix file and
# output the FAC transformation matrix

# OUTPUT:
# [X_ECEF, Y_ECEF, Z_ECEF] *transform_matrix = [r, e, p]



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

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = [4, 5]
outputData = True
outputPath_modifier = 'coordinates' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---

def ECEF_to_Field_Aligned(wRocket):


    rocketID = ACESII.payload_IDs[wRocket-4]

    # --- ACES II Flight/Integration Data ---
    inputFiles = glob(f'{DataPaths.ACES_data_folder}\coordinates\{ACESII.fliers[wRocket-4]}\*ECEF_to_ENU*')[0]

    # --- get the data from the attitude file ---
    stl.prgMsg(f'Loading data')
    data_dict_transform = stl.loadDictFromFile(inputFiles)
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
        'ECEF_to_FAC': [np.zeros(shape=(len(data_dict_transform['Epoch'][0]), 3, 3)), {}],
        'FAC_to_ECEF': [np.zeros(shape=(len(data_dict_transform['Epoch'][0]), 3, 3)), {}],
        'Epoch': deepcopy(data_dict_transform['Epoch']),
        'Alt': deepcopy(data_dict_transform['Alt']),
        'Lat': deepcopy(data_dict_transform['Lat']),
        'Long': deepcopy(data_dict_transform['Long'])
    }

    # #############################################################
    # --- Convert ENU coordinates to Field Aligned coordinates ---
    # #############################################################
    stl.prgMsg(f'Loading CHAOS model')

    # Get the Data
    B_model = stl.CHAOS(lat=data_dict_transform['Lat'][0],
                    long=data_dict_transform['Long'][0],
                    alt=data_dict_transform['Alt'][0],
                    times=data_dict_transform['Epoch'][0])  # CHAOS in ENU coordinates

    B_CHAOS_ECEF = np.array([np.matmul(data_dict_transform['ENU_to_ECEF'][0][i], B_model[i]) for i in range(len(data_dict_transform['Epoch'][0]))])
    stl.Done(start_time)

    # --- determine the Payload's transformation matrix to ECEF ---
    R_REF = 6371.2  # earth Radius in km
    Radius = data_dict_transform['Alt'][0] + R_REF
    coLatRad = [np.radians(90 - lat) for lat in data_dict_transform['Lat'][0]]
    LongRad = [np.radians(long) for long in data_dict_transform['Long'][0]]
    Rsc = np.array([
        [Radius[i] * np.sin(coLatRad[i]) * np.cos(LongRad[i]),
         Radius[i] * np.sin(coLatRad[i]) * np.sin(LongRad[i]),
         Radius[i] * np.cos(coLatRad[i])] for i in range(len(data_dict_transform['Epoch'][0]))])

    stl.Done(start_time)


    # --- calculate Field Aligned unit vectors over the duration of the flight ---
    stl.prgMsg('Converting to Field Aligned Coordinates')

    # pHat comes from the CHAOS model direction of B in GEO
    pHat = np.array([B_CHAOS_ECEF[i] / np.linalg.norm(B_CHAOS_ECEF[i]) for i in range(len(data_dict_transform['Epoch'][0]))])

    # e-hat comes from the cross of pHat and the Rocket's radius vector (in geomagnetic coordinates)
    eHat = np.array([np.cross(pHat[i], Rsc[i]) / np.linalg.norm(np.cross(pHat[i], Rsc[i])) for i in range(len(data_dict_transform['Epoch'][0]))])

    # rHat comes from the cross of eHat and pHat
    rHat = np.array([np.cross(eHat[i], pHat[i]) for i in range(len(data_dict_transform['Epoch'][0]))])

    # form the transformation matrix FROM GEO TO FIELD ALIGNED
    ECEF_to_FAC_transformation = np.array([[rHat[i],eHat[i], pHat[i]] for i in range(len(data_dict_transform['Epoch'][0]))])
    stl.Done(start_time)

    # store the outputs
    data_dict_output['a11'][0] = ECEF_to_FAC_transformation[:, 0, 0]
    data_dict_output['a12'][0] = ECEF_to_FAC_transformation[:, 0, 1]
    data_dict_output['a13'][0] = ECEF_to_FAC_transformation[:, 0, 2]
    data_dict_output['a21'][0] = ECEF_to_FAC_transformation[:, 1, 0]
    data_dict_output['a22'][0] = ECEF_to_FAC_transformation[:, 1, 1]
    data_dict_output['a23'][0] = ECEF_to_FAC_transformation[:, 1, 2]
    data_dict_output['a31'][0] = ECEF_to_FAC_transformation[:, 2, 0]
    data_dict_output['a32'][0] = ECEF_to_FAC_transformation[:, 2, 1]
    data_dict_output['a33'][0] = ECEF_to_FAC_transformation[:, 2, 2]
    data_dict_output['ECEF_to_FAC'][0] = ECEF_to_FAC_transformation
    data_dict_output['FAC_to_ECEF'][0] = np.array([ECEF_to_FAC_transformation[i].T for i in range(len(data_dict_output['Epoch'][0]))])

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

        fileoutName = f'ACESII_{rocketID}_ECEF_to_FAC.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}{outputPath_modifier}\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
for val in wRocket:
    ECEF_to_Field_Aligned(val)