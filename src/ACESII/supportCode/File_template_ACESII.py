# --- File_Template_ACESII.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:

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
wInstr = 'EEPAA'
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['L2', [[0],[0]]],
    'Lshell':['coordinates',[[0],[0]]]
}
outputData = True


#################
# --- IMPORTS ---
#################
from scipy.integrate import simpson
from src.ACESII.data_tools.my_imports import *
start_time = time.time()


def file_template(data_dicts):

    #######################
    # --- LOAD THE DATA ---
    #######################
    stl.prgMsg(f'Loading data')
    data_dict_ESA = deepcopy(data_dicts[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
                       }


    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')
        # write out the data
        fileoutName = f'ACESII_{ACESII.fliers_dict[rocket_str]}_l3_{wInstr.lower()}_energy_flux.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/L3/{wInstr}/{rocket_str}//{fileoutName}'
        stl.outputDataDict(outputPath, data_dict=data_dict_output)
        stl.Done(start_time)


#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(file_template,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)
