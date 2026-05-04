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
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['L2', [[0],[0]]],
    'trajectories':['',[[2],[1]]]
}
outputData = True


fs = 256 # sample frequency of the data
low_cutoff = 0.05 # 20 seconds (or 1/20 freq) butterworth cutoff
order = 4

#################
# --- IMPORTS ---
#################
from scipy.integrate import simpson
import datetime as dt
from scipy.interpolate import CubicSpline
from src.ACESII.data_tools.my_imports import *
start_time = time.time()


def L2_MAG_to_Ampere_Current(data_dicts):

    #######################
    # --- LOAD THE DATA ---
    #######################
    stl.prgMsg(f'Loading data')
    data_dict_MAG = deepcopy(data_dicts[0])
    data_dict_traj = deepcopy(data_dicts[1])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
        'J_para':[[],{'DEPEND_0':'Epoch','UNITS':'&mu;A/m!A-2!N','LABLAXIS':'Parallel Current','VAR_TYPE':'data'}],
        'DeltaT':[[], {'DEPEND_0':'Epoch','UNITS':'A/m!A-2!N','LABLAXIS':'Parallel Current'}],
        'DeltaBN': [[], {'DEPEND_0': 'Epoch', 'UNITS': None, 'LABLAXIS': '&Delta; BN'}],
        'DeltaBT': [[], {'DEPEND_0': 'Epoch', 'UNITS': None, 'LABLAXIS': '&Delta; BT'}],
        'DeltaBp': [[], {'DEPEND_0': 'Epoch', 'UNITS': None, 'LABLAXIS': '&Delta; BT'}],
        'DeltaVN': [[], {'DEPEND_0': 'Epoch', 'UNITS': None, 'LABLAXIS': '&Delta; VN'}],
        'DeltaVT': [[], {'DEPEND_0': 'Epoch', 'UNITS': None, 'LABLAXIS': '&Delta; BT'}],
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
    for key in ['N_VEL','T_VEL','P_VEL']:
        cs=CubicSpline(T0_traj, data_dict_traj[key][0])
        data_dict_traj[key][0] = cs(T0_bmed)

    # Calculate the deltas
    DeltaT = np.diff(T0_bmed,prepend=T0_bmed[0])
    Bkeys = ['B_N','B_T','B_p']
    dB_filtered = [[],[],[]]
    for i,key in enumerate(Bkeys):
        dB = (1E-9)*np.diff(data_dict_MAG[key][0], prepend=data_dict_MAG[key][0][-1])

        # Filter the deltas
        dB_filtered[i] = stl.butterFilter().butter_filter(data=dB,
                                                    lowcutoff=low_cutoff,
                                                    highcutoff=low_cutoff,
                                                    order=order,
                                                    fs=fs,
                                                    filtertype='lowPass'
                                                    )
    dB_filtered = np.array(dB_filtered)


    # Calculate J_parallel and Store
    data_dict_output['J_para'][0] = (1/1E-6)*(1/(stl.u0*DeltaT)) * (dB_filtered[1]/data_dict_traj['N_VEL'][0] - dB_filtered[0]/data_dict_traj['T_VEL'][0])
    data_dict_output['DeltaT'][0] = DeltaT
    data_dict_output['DeltaBN'][0] = dB_filtered[0]
    data_dict_output['DeltaBT'][0] = dB_filtered[1]
    data_dict_output['DeltaBp'][0] = dB_filtered[2]

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')
        # write out the data
        fileoutName = f'ACESII_{ACESII.fliers_dict[rocket_str]}_jpara.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/science/Currents/{rocket_str}/{fileoutName}'
        stl.outputDataDict(outputPath, data_dict=data_dict_output)
        stl.Done(start_time)


#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(L2_MAG_to_Ampere_Current,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)
