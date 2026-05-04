# --- File_Template_ACESII.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2026-03-30"
__version__ = "1.0.0"

import numpy as np

# --- --- --- --- ---


######################
# --- DATA TOGGLES ---
######################
just_print_file_names_bool = False
rocket_str = 'low'
wInstr = 'EFI'
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['L2', [[0],[2]]],
}
outputData = True


#################
# --- IMPORTS ---
#################
from scipy.integrate import simpson
from src.ACESII.data_tools.my_imports import *
import ppigrf
start_time = time.time()


def file_template(data_dicts):

    #######################
    # --- LOAD THE DATA ---
    #######################
    stl.prgMsg(f'Loading data')
    data_dict_EFI = deepcopy(data_dicts[0])
    stl.Done(start_time)


    # --- Calculate the CHAOS Model --
    stl.prgMsg('Calculating IGRF model')
    lats = data_dict_EFI['Lat'][0]
    lons = data_dict_EFI['Long'][0]
    alts = data_dict_EFI['Alt'][0]
    date = data_dict_EFI['Epoch'][0][0]
    Be,Bn,Bu = ppigrf.igrf(lons,lats,alts,date)
    B_MODEL = 1E-9*np.array([Be[0],Bn[0],Bu[0]]).T
    stl.Done(start_time)

    # --- Calculate ExB ----
    def wCoordUsed(dict_keys):
        coord_trio = [['N', 'T', 'p'], ['E', 'N', 'Up'], ['r', 'e', 'p'], ['X', 'Y', 'Z']]

        for trio in coord_trio:
            coord_keys = ['E_' + strV for strV in trio]
            if all(target in dict_keys for target in coord_keys):
                return trio

        raise Exception('Could not determine Coordinate System')
    trio = wCoordUsed(data_dict_EFI.keys())
    E_keys = ['E_'+val for val in trio]
    EField = np.array([data_dict_EFI[f'{key}'][0] for key in E_keys]).T
    ExB = np.cross(EField,B_MODEL)
    Bmag = np.square(np.array([np.linalg.norm(B_MODEL[i]) for i in range(len(data_dict_EFI['Epoch'][0]))]))
    ExB_drift = np.array([ExB[:,i]/Bmag for i in range(3)]).T
    data_dict_output = {
        f'ExB_{val}' : [ExB_drift[:,i],{'UNITS':'m/s','VAR_TYPE':'data','LABLAXIS':f'ExB_{val}'}] for i,val in enumerate(trio)
    }
    data_dict_output={**data_dict_output,**{'Epoch':deepcopy(data_dict_EFI['Epoch'])}}

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')
        # write out the data
        fileoutName = f'ACESII_{ACESII.fliers_dict[rocket_str]}_ExB.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/science/ExB/{rocket_str}//{fileoutName}'
        stl.outputDataDict(outputPath, data_dict=data_dict_output)
        stl.Done(start_time)


#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(file_template,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)
