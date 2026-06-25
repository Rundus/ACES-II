# --- L3_calculate_Poynting_Flux.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Determine the PoyntingFLux of the data using E-Field and B-Field fullcal Measurements
# for the low flyer


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
# --- --- --- --- ---
#
######################
# --- DATA TOGGLES ---
######################
just_print_file_names_bool = False
rocket_str = 'high'

# auroral
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    'EFI':['L3', [[2],[2]]],
    'MAG':['L3',[[1],[2]]],
}

# FAC
# dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
#     f'EFI':['L3', [[1],[1]]],
#     'MAG':['L3',[[2],[1]]]
# }

# ENU
# dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
#     f'EFI':['L3', [[0],[0]]],
#     'MAG':['L3',[[0],[3]]]
# }

outputData = True


#################
# --- IMPORTS ---
#################
from src.ACESII.data_tools.my_imports import *
start_time = time.time()
import datetime as dt
from scipy.interpolate import CubicSpline


def L3_calculate_Poynting_Flux(data_dicts):

    #######################
    # --- LOAD THE DATA ---
    #######################
    stl.prgMsg(f'Loading data')
    data_dict_EFI = deepcopy(data_dicts[0])
    print(data_dict_EFI.keys())
    data_dict_MAG = deepcopy(data_dicts[1])
    print(data_dict_MAG.keys())
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
                       }

    # include the ephemeris data
    for key in ['Epoch','ILat','L-Shell','ILong','Lat','Alt','Long']:
        if key in data_dict_EFI.keys():
            data_dict_output = {**data_dict_output,**{f'{key}':deepcopy(data_dict_EFI[key])}}

    # Interpolate the MAG data onto the EFI timebase
    stl.prgMsg('Interpolating B Residuals onto EFI')
    T0 = dt.datetime(2022,11,20,17,20)
    T0_MAG = stl.EpochTo_T0_Rocket(data_dict_MAG['Epoch'][0],T0=T0)
    T0_EFI = stl.EpochTo_T0_Rocket(data_dict_EFI['Epoch'][0],T0=T0)

    def wCoordUsed(dict_keys):
        coord_trio = [['N', 'T', 'p'], ['E', 'N', 'U'], ['r', 'e', 'p'], ['X', 'Y', 'Z']]
        wCoords = ['auroral','ENU','FAC','ECEF']

        for i,trio in enumerate(coord_trio):
            E_keys = ['E_' + strV for strV in trio]
            if all(target in dict_keys for target in E_keys):
                return trio, wCoords[i]

        raise Exception('Could not determine Coordinate System')

    coord_keys, wCoords = wCoordUsed(data_dict_EFI.keys())
    B_keys = ['B_' + strV+'_residual' for strV in coord_keys]
    E_keys = ['E_' + strV for strV in coord_keys]

    B_residual_interp = [[],[],[]]

    for i,keyVal in enumerate(B_keys):
        cs = CubicSpline(T0_MAG, data_dict_MAG[keyVal][0])
        B_residual_interp[i] = cs(T0_EFI)

    B_residual_interp = np.array(B_residual_interp).T
    E_Field = np.array([data_dict_EFI[keyVal][0] for keyVal in E_keys]).T
    stl.Done(start_time)

    # --- Calculate the Poynting Flux ---
    stl.prgMsg('Calculating Poynting Flux')
    S = np.cross(E_Field, (1E-9)*B_residual_interp,axis=1)/stl.u0

    for i, keyVals in enumerate(coord_keys):
        data_dict_output = {
            **data_dict_output,
            **{
                f'S_{keyVals}':[S[:,i],{'DEPEND_0':'Epoch','VAR_TYPE':'data','UNITS':'W/m!A2!N','LABLAXIS':f'&delta;S_{keyVals}'}]
            }
        }
    stl.Done(start_time)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')
        # write out the data
        fileoutName = f'ACESII_{ACESII.fliers_dict[rocket_str]}_poynting_flux_{wCoords}_DC.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/science/PoyntingFlux/{rocket_str}//{fileoutName}'
        stl.outputDataDict(outputPath, data_dict=data_dict_output)
        stl.Done(start_time)


#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(L3_calculate_Poynting_Flux,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)