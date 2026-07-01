# --- MAG_L2_to_L3_residuals.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Take the L2 MAG data and subtract off the CHAOS magnetic field model.
# Then filter the result to show the AC vs DC components

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2026-03-30"
__version__ = "1.0.0"

import matplotlib.pyplot as plt

# --- --- --- --- ---

######################
# --- DATA TOGGLES ---
######################
just_print_file_names_bool = False
rocket_str = 'low'
wInstr = 'MAG'
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['L2', [[0],[1]]],
}

# if you wish to detrend the data after a certain point
detrend_data = False
target_ILat = {'high':70.45,
               'low': 70.9}

plot_smoothing = False
outputData = False

fs = 256 # sample frequency of the data
low_cutoff = {'high':0.05,'low':0.1} # [Best value HF/LF: 0.05/0.1] 20 seconds (or 1/20 freq) butterworth cutoff
order = 4
filtType = 'lowpass'


#################
# --- IMPORTS ---
#################
from scipy.signal import detrend,savgol_filter
from src.ACESII.data_tools.my_imports import *
from scipy.interpolate import UnivariateSpline
import datetime as dt
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

        import re
        coord_axes = [['N', 'T', 'p'], ['E', 'N', 'U'], ['r', 'e', 'p'], ['X', 'Y', 'Z']]
        coordinate_name = ['auroral', 'ENU', 'FAC', 'ECEF']

        # use the dictonary keys to determine the coordinates
        pat = re.compile(r'(?<=_)(E|N|U|D|X|Y|Z|T|p|r|e)')
        axes = set()
        for s in dict_keys:
            m = pat.search(s)
            if m:
                axes.add(m.group(1))

        match = next((lst for lst in coord_axes if set(lst) == axes), None)
        cordNames = coordinate_name[coord_axes.index(match)]

        if len(match) >= 1:
            return match, cordNames
        else:
            raise Exception('Could not determine Coordinate System')



    coord_keys, wCoords = wCoordUsed(data_dict_MAG.keys())
    B_keys = ['B_' + strV for strV in coord_keys]

    # Prepare the output data dict
    data_dict_output = {**data_dict_output,
                        **{
                            f'B_{coord_keys[i]}_residual': [[], {'DEPEND_0': 'Epoch', 'UNITS': None,'VAR_TYPE':'data', 'LABLAXIS': f'&delta; B{coord_keys[i]}'}] for i,key in enumerate(coord_keys)
                        }}

    data_dict_output = {**data_dict_output,
                        **{
                            f'B_{coord_keys[i]}_smooth': [[], {'DEPEND_0': 'Epoch', 'UNITS': None, 'VAR_TYPE': 'data', 'LABLAXIS': f'B{coord_keys[i]} (smooth)'}] for i, key in enumerate(coord_keys)
                        }}

    ##########################################
    # --- Calculate the local-smooth model ---
    ##########################################
    stl.prgMsg('Smoothing Data')
    # reduce the data to only the relevant regions
    low_idx = np.abs(data_dict_MAG['ILat'][0] - target_ILat[rocket_str]).argmin()

    # smooth the data
    for i, key in enumerate(B_keys):

        # === First Fill the smooth data container with raw data ===
        data_dict_output[f'{key}_smooth'][0] = deepcopy(data_dict_MAG[f'{key}'][0])

        # === Savitz-golay a portion of the data ===
        filter_length =len(data_dict_MAG[f'{key}'][0][low_idx:])
        data_dict_output[f'{key}_smooth'][0][low_idx:] = savgol_filter(x=deepcopy(data_dict_MAG[f'{key}'][0][low_idx:]),
                                                               polyorder=3,
                                                               window_length=filter_length)
    stl.Done(start_time)


    if plot_smoothing:
        fig, ax = plt.subplots(3)

        for i,key in enumerate(B_keys):
            ax[i].plot(data_dict_MAG['Epoch'][0],data_dict_MAG[f'{key}'][0])
            ax[i].plot(data_dict_MAG['Epoch'][0], data_dict_output[f'{key}_smooth'][0],color='red')

            if i == 0:
                lims = [-300,300]
            elif i==1:
                lims = [-500, 500]
            elif i==2:
                lims = [45000, 52000]
            ax[i].set_ylim(*lims)

        plt.show()

    #################################
    # --- Calculate the residuals ---
    #################################
    stl.prgMsg('Filtering Data')

    # calculate the difference from model: shows the electrodynamics. Should be B_measured - B_model
    for i, key in enumerate(B_keys):

        data = data_dict_MAG[f'{key}'][0][low_idx:] - data_dict_output[f'{key}_smooth'][0][low_idx:]

        # SSA Filter the spin/coning
        # wL=2000
        # SSAobj = stl.SSA(data,L=wL)
        # SSAobj.components_to_df().plot()
        # SSAobj.orig_TS.plot(alpha=0.4)
        # plt.xlabel('$t$')
        # plt.ylabel('F')
        # plt.show()


        # butterworth Filter away the spin/coning
        data_dict_output[f'{key}_residual'][0] = stl.butterFilter().butter_filter(
            # data=data_dict_MAG[f'{key}'][0] - data_dict_MAG[f'B_model_{coord_keys[i]}'][0],
            data=data - data_dict_output[f'{key}_smooth'][0],
            # data=data_dict_MAG[f'{key}'][0],
            lowcutoff=low_cutoff[rocket_str],
            highcutoff=low_cutoff[rocket_str],
            order=order,
            fs=fs,
            filtertype=filtType)
    stl.Done(start_time)

    if detrend_data:

        # === [0] LINEAR DETREND ===
        # find the relevant index
        target_idx = np.abs(data_dict_MAG['ILat'][0] - target_ILat[rocket_str]).argmin()

        for i, key in enumerate(B_keys):
            data = data_dict_output[f'{key}_residual'][0][target_idx::]
            data_dict_output[f'{key}_residual'][0][target_idx::] = detrend(data,type='linear')

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