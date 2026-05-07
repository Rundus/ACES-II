# --- MAG_L2_to_L3.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Take the L2 MAG data and subtract off the CHAOS magnetic field model.
# Then filter the result to show the AC vs DC components

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2026-03-30"
__version__ = "1.0.0"
# --- --- --- --- ---

######################
# --- DATA TOGGLES ---
######################
just_print_file_names_bool = False
rocket_str = 'low'
wInstr = 'MAG'
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['L2', [[3],[3]]],
}
outputData = True

fs = 256 # sample frequency of the data
low_cutoff = 0.05 # 20 seconds (or 1/20 freq) butterworth cutoff
order = 4
filtType = 'lowpass'


#################
# --- IMPORTS ---
#################
from scipy.integrate import simpson
from src.ACESII.data_tools.my_imports import *
start_time = time.time()


def MAG_L2_to_L3(data_dicts):

    #######################
    # --- LOAD THE DATA ---
    #######################
    stl.prgMsg(f'Loading data')
    data_dict_MAG = deepcopy(data_dicts[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
    }
    for key in ['Epoch','ILat','L-Shell','ILong','Lat','Alt','Long']:
        if key in data_dict_MAG.keys():
            data_dict_output = {**data_dict_output,**{f'{key}':deepcopy(data_dict_MAG[key])}}

    # Determine the coordinates of the datafile
    def wCoordUsed(dict_keys):
        coord_trio = [['N', 'T', 'p'], ['E', 'N', 'U'], ['r', 'e', 'p'], ['X', 'Y', 'Z']]
        wCoords = ['auroral','ENU','FAC','ECEF']

        for i,trio in enumerate(coord_trio):
            B_keys = ['B_' + strV for strV in trio]
            if all(target in dict_keys for target in B_keys):
                return trio, wCoords[i]

        raise Exception('Could not determine Coordinate System')

    coord_keys, wCoords = wCoordUsed(data_dict_MAG.keys())
    B_keys = ['B_' + strV for strV in coord_keys]

    # Prepare the output data dict
    data_dict_output = {**data_dict_output,
                        **{
                            f'B_{coord_keys[i]}_residual': [[], {'DEPEND_0': 'Epoch', 'UNITS': None,'VAR_TYPE':'data', 'LABLAXIS': f'&delta; B{coord_keys[i]}'}] for i,key in enumerate(coord_keys)
                        }}


    # Calculate the residuals

    # calculate the difference from model: shows the electrodynamics
    for i, key in enumerate(B_keys):
        data_dict_output[f'{key}_residual'][0] = stl.butterFilter().butter_filter(
            data=data_dict_MAG[f'B_model_{coord_keys[i]}'][0] - data_dict_MAG[f'{key}'][0],
            lowcutoff=low_cutoff,
            highcutoff=low_cutoff,
            order=order,
            fs=fs,
            filtertype=filtType)



    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')
        # write out the data
        fileoutName = f'ACESII_{ACESII.fliers_dict[rocket_str]}_l3_RingCore_{wCoords}_residuals.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/L3/{wInstr}/{rocket_str}//{fileoutName}'
        stl.outputDataDict(outputPath, data_dict=data_dict_output)
        stl.Done(start_time)


#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(MAG_L2_to_L3,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)