# --- L2_EFI_to_L3_EFI_proxy.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Use the Low Flyer EFI data and project it up into the High Flyer using the invariant latitude
# coordinates

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2026-03-30"
__version__ = "1.0.0"

from src.ACESII.Science.LangmuirProbes.L3_LP_fit_swept_parameters import data_dict_output

# --- --- --- --- ---


######################
# --- DATA TOGGLES ---
######################
just_print_file_names_bool = True
rocket_str = 'low'
wInstr = 'EFI'
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['L2', [[0],[1]]],
    '':['coordinates/invariant_coordinates',[[0],[0]]]
}
outputData = True


#################
# --- IMPORTS ---
#################
from scipy.integrate import simpson
from src.ACESII.data_tools.my_imports import *
start_time = time.time()
import datetime as dt


def L2_EFI_to_L3_EFI_proxy(data_dicts):

    #######################
    # --- LOAD THE DATA ---
    #######################
    stl.prgMsg(f'Loading data')
    data_dict_EFI = deepcopy(data_dicts[0])
    data_dict_invarCoords = deepcopy(data_dicts[1])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {**{}, **deepcopy(data_dict_invarCoords)}

    # --- determine the coordinate system used ---

    def wCoordUsed(dict_keys):
        coord_trio = [['N', 'T', 'p'], ['E', 'N', 'U'], ['r', 'e', 'p'], ['X', 'Y', 'Z']]
        wCoords = ['auroral', 'ENU', 'FAC', 'ECEF']

        for i, trio in enumerate(coord_trio):
            E_keys = ['E_' + strV for strV in trio]
            if all(target in dict_keys for target in E_keys):
                return trio, wCoords[i]

        raise Exception('Could not determine Coordinate System')

    coord_keys, wCoords = wCoordUsed(data_dict_EFI.keys())
    E_keys = ['E_' + strV for strV in coord_keys]

    # Create a proxy dataset for the High Flyer
    if rocket_str == 'high':
        E_Field_interp = [[],[],[]]
        for i, keyVal in enumerate(E_keys):
            E_Field_interp[i] = np.interp(data_dict_invarCoords['ILat'][0],data_dict_EFI['ILat'][0], data_dict_EFI[keyVal][0])
        E_Field_interp = np.array(E_Field_interp).T

        for i, keyVals in enumerate(coord_keys):
            data_dict_output = {
                **data_dict_output,
                **{
                    f'E_{keyVals}':[E_Field_interp[:,i],{'DEPEND_0':'Epoch','VAR_TYPE':'data','UNITS':'V/m','LABLAXIS':f'E_{keyVals}'}]
                }
            }

        fileoutName = f'ACESII_{ACESII.fliers_dict['high']}_l3_{wInstr}_proxy_{wCoords}.cdf'
    else:
        data_dict_output = deepcopy(data_dict_EFI)
        fileoutName = f'ACESII_{ACESII.fliers_dict['high']}_l3_{wInstr}_{wCoords}.cdf'
    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')
        # write out the data

        outputPath = f'{DataPaths.ACES_data_folder}/L3/{wInstr}/{rocket_str}//{fileoutName}'
        stl.outputDataDict(outputPath, data_dict=data_dict_output)
        stl.Done(start_time)


#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(L2_EFI_to_L3_EFI_proxy,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)
