# --- L2_MAG_to_Ampere_Current.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Calculate the Parallel Current using Ampere's Law

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2026-03-30"
__version__ = "1.0.0"
import numpy as np
import spaceToolsLib
# --- --- --- --- ---

######################
# --- DATA TOGGLES ---
######################
just_print_file_names_bool = False
rocket_str = 'high'
wInstr = 'MAG'

# Auroral
# dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
#     f'{wInstr}':['L3', [[1],[2]]],
#     'trajectories':['',[[3],[2]]],
# }

# FAC
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['L3', [[2],[1]]],
    'trajectories':['',[[1],[0]]],
}
outputData = True

#################
# --- IMPORTS ---
#################
from scipy.integrate import simpson
import datetime as dt
from scipy.interpolate import CubicSpline
from src.ACESII.data_tools.my_imports import *
start_time = time.time()


def L3_MAG_to_Ampere_Current(data_dicts):

    #######################
    # --- LOAD THE DATA ---
    #######################
    stl.prgMsg(f'Loading data')
    data_dict_MAG = deepcopy(data_dicts[0])
    data_dict_traj = deepcopy(data_dicts[1])
    stl.Done(start_time)


    # --- [00] Determine the coordinate system used ---
    def wCoordUsed(dict_keys):
        import re
        coord_axes = [['N', 'T', 'p'], ['E', 'N', 'U'], ['r', 'e', 'p'], ['X', 'Y', 'Z']]
        coordinate_name = ['auroral', 'ENU', 'FAC', 'ECEF']

        # use the dictonary keys to determine the coordinates
        pat = re.compile(r'(?<=_)(E|N|U|D|X|Y|Z|T|p|r|e)(?=_)')
        axes = set()
        for s in dict_keys:
            m = pat.search(s)
            if m:
                axes.add(m.group(1))

        match = next((lst for lst in coord_axes if set(lst) == axes), None)
        cordNames = coordinate_name[coord_axes.index(match)]

        if len(match)>=1:
            return match,cordNames
        else:
            raise Exception('Could not determine Coordinate System')

    coord_keys,coordinate_name = wCoordUsed(data_dict_MAG.keys())
    print(coordinate_name,coord_keys)

    # --- prepare the output ---
    data_dict_output = {
        f'J_{coord_keys[2]}':[[],{'DEPEND_0':'Epoch','UNITS':'&mu;A/m!A-2!N','LABLAXIS':'Field-aligned Current','VAR_TYPE':'data'}],
        f'DeltaT':[[], {'DEPEND_0':'Epoch','UNITS':'A/m!A-2!N','LABLAXIS':'Parallel Current'}],
        f'dB_{coord_keys[0]}': [[], {'DEPEND_0': 'Epoch', 'UNITS': None, 'LABLAXIS': f'&Delta; B{coord_keys[0]}'}],
        f'dB_{coord_keys[1]}': [[], {'DEPEND_0': 'Epoch', 'UNITS': None, 'LABLAXIS': f'&Delta; B{coord_keys[1]}'}],
        f'dB_{coord_keys[2]}': [[], {'DEPEND_0': 'Epoch', 'UNITS': None, 'LABLAXIS': f'&Delta; B{coord_keys[2]}'}],
        f'dB{coord_keys[0]}d{coord_keys[1]}': [[], {'DEPEND_0': 'Epoch', 'UNITS': None, 'LABLAXIS': f'dB{coord_keys[0]}d{coord_keys[1]}'}],
        f'dB{coord_keys[1]}d{coord_keys[0]}': [[], {'DEPEND_0': 'Epoch', 'UNITS': None, 'LABLAXIS': f'dB{coord_keys[1]}d{coord_keys[0]}'}],
        'Epoch':deepcopy(data_dict_MAG['Epoch']),
        'ILat':deepcopy(data_dict_MAG['ILat']),
        'L-Shell': deepcopy(data_dict_MAG['L-Shell']),
                       }
    data_dict_output['L-Shell'][1]['VAR_TYPE']='support_data'

    # --- [0] Interpolate Trajectory data onto MAG timebase ---
    # Get the Trajectory on the B-Field timebase
    T0 = dt.datetime(2022,11,20,17,20)
    T0_bmed = stl.EpochTo_T0_Rocket(data_dict_MAG['Epoch'][0],T0=T0)
    T0_traj = spaceToolsLib.EpochTo_T0_Rocket(data_dict_traj['Epoch'][0],T0=T0)
    for key in [f'{coord_keys[i].upper()}_VEL' for i in range(len(coord_keys))]:
        cs=CubicSpline(T0_traj, data_dict_traj[key][0])
        data_dict_traj[key][0] = cs(T0_bmed)

    # --- [1] Calculate the deltaT variable ---
    DeltaT = np.diff(T0_bmed,prepend=T0_bmed[0])

    # --- [2] Calculate the dB data then Filter/Smooth ---

    dB= [[],[],[]]
    for i,key in enumerate(coord_keys):
        keyVal = 'B_'+key +'_residual'
        dB[i] = np.array((1E-9)*np.diff(data_dict_MAG[keyVal][0], prepend=data_dict_MAG[keyVal][0][-1]))

    # Calculate J_parallel and Store
    data_dict_output[f'J_{coord_keys[2]}'][0] = (1/1E-6)*(1/(stl.u0*DeltaT)) * (dB[1]/data_dict_traj[f'{coord_keys[0].upper()}_VEL'][0] - dB[0]/data_dict_traj[f'{coord_keys[1].upper()}_VEL'][0])

    data_dict_output[f'dB{coord_keys[0]}d{coord_keys[1]}'][0] =  dB[0] / data_dict_traj[f'{coord_keys[1].upper()}_VEL'][0]
    data_dict_output[f'dB{coord_keys[1]}d{coord_keys[0]}'][0] = dB[1] / data_dict_traj[f'{coord_keys[0].upper()}_VEL'][0]

    data_dict_output['DeltaT'][0] = DeltaT
    data_dict_output[f'dB_{coord_keys[0]}'][0] = dB[0]
    data_dict_output[f'dB_{coord_keys[1]}'][0] = dB[1]
    data_dict_output[f'dB_{coord_keys[2]}'][0] = dB[2]

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')
        # write out the data
        fileoutName = f'ACESII_{ACESII.fliers_dict[rocket_str]}_Ampere_Currents_{coordinate_name}.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/science/Currents/{rocket_str}/{fileoutName}'
        stl.outputDataDict(outputPath, data_dict=data_dict_output)
        stl.Done(start_time)


#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(L3_MAG_to_Ampere_Current,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)
