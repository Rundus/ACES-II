# --- L0_to_L1.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert electrostatic analyzer data from diffNFlux to Distribution Function



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import itertools
# --- --- --- --- ---
import time
import spaceToolsLib as stl
start_time = time.time()
# --- --- --- --- ---




######################
# --- DATA TOGGLES ---
######################
just_print_file_names_bool = False
rocket_str = 'high'
# wInstr = 'EEPAA'
wInstr = 'IEPAA'
# wInstr = 'LEESA'
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['L2', [[0],[0]]],
    f'':['science/ion_mass_effective/',[[0],[0]]]
}
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.ACESII.data_tools.my_imports import *
from scipy.interpolate import CubicSpline

def L2_to_DistFunc(data_dicts):

    # --- get the data ---
    stl.prgMsg(f'Loading data')
    data_dict_esa = deepcopy(data_dicts[0])
    stl.Done(start_time)

    stl.prgMsg('Loading data from L2Files')
    data_dict_meff_i = deepcopy(data_dicts[1])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}

    # --- --- --- --- --- --- --- --- ---
    # --- Calculate Instrument Data ---
    # --- --- --- --- --- --- --- --- ---
    stl.prgMsg('Calculating the Distribution Function')

    # --- CALCULATE DISTRIBUTION FUNCTION ---
    diffNFlux = data_dict_esa['Differential_Number_Flux'][0]
    # oneCountLevel = data_dict_esa['oneCountLevel'][0]
    Energies = data_dict_esa['Energy'][0]

    # define empty numpy array
    sizes = [len(diffNFlux),len(diffNFlux[0]), len(diffNFlux[0][0])]
    ranges = [range(sizes[0]), range(sizes[1]), range(sizes[2])]
    distFunc = np.zeros(shape=(sizes[0], sizes[1], sizes[2]))
    # distFunc_oneCount = np.zeros(shape=(sizes[0], sizes[1], sizes[2]))

    print(f'Num. of iterations: {sizes[0]*sizes[1]*sizes[2]}',end='\n')

    # --- Calculate DistFunc in SI units ---
    if wInstr[1] in ['eepaa','leesa']:
        mass = [stl.m_e for i in range(len(diffNFlux))]
    else:
        Epoch_meff_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_meff_i['Epoch'][0]])
        Epoch_esa_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_esa['Epoch'][0]])
        cs = CubicSpline(Epoch_meff_tt2000, data_dict_meff_i['m_eff_i'][0])
        mass = cs(Epoch_esa_tt2000)

    for tme, ptch, engy in tqdm(itertools.product(*ranges)):

        if diffNFlux[tme][ptch][engy] <= ACESII.epoch_fillVal:
            distFunc[tme][ptch][engy] = ACESII.epoch_fillVal
        else:
            distVal = (stl.cm_to_m*stl.cm_to_m/(stl.q0*stl.q0 ))*(((mass[tme]**2)*diffNFlux[tme][ptch][engy]) / (2 * Energies[engy]))
            if distVal < 0:
                distFunc[tme][ptch][engy] = 0
            else:
                distFunc[tme][ptch][engy] = distVal

        # distFunc_oneCount[tme][ptch][engy] = (stl.cm_to_m*stl.cm_to_m/(stl.q0*stl.q0 ))*(((mass[tme]**2)*oneCountLevel[tme][ptch][engy]) / (2 * Energies[engy]))

    stl.Done(start_time)


    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        stl.prgMsg('Creating output file')
        outputPath = f'{DataPaths.ACES_data_folder}/L3/{wInstr}/{rocket_str}//ACESII_{ACESII.fliers_dict[rocket_str]}_distFunc_{wInstr}.cdf'

        data_dict_output = {**data_dict_output, **{'Distribution_Function':
                                         [distFunc, {'LABLAXIS': 'Distribution_Function',
                                                   'DEPEND_0': 'Epoch',
                                                   'DEPEND_1': 'Pitch_Angle',
                                                   'DEPEND_2': 'Energy',
                                                   'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2',
                                                   'UNITS': 'm!A-6!Ns!A3!N',
                                                   'VALIDMIN': distFunc.min(), 'VALIDMAX': distFunc.max(),
                                                   'VAR_TYPE': 'data', 'SCALETYP': 'log'}]}}

        for key in ['Epoch','Energy','Pitch_Angle','Alt','Lat','Long']:
            data_dict_output = {**data_dict_output, **{f'{key}':deepcopy(data_dict_esa[f'{key}'])}}

        # data_dict = {**data_dict_output, **{'oneCountLevel':
        #                                  [distFunc_oneCount, {'LABLAXIS': 'Distribution_Function',
        #                                              'DEPEND_0': 'Epoch',
        #                                              'DEPEND_1': 'Pitch_Angle',
        #                                              'DEPEND_2': 'Energy',
        #                                              'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2',
        #                                              'UNITS': 'm!A-6!Ns!A3!N',
        #                                              'VALIDMIN': distFunc.min(), 'VALIDMAX': distFunc.max(),
        #                                              'VAR_TYPE': 'support_data', 'SCALETYP': 'log'}]}}

        stl.outputDataDict(outputPath=outputPath, data_dict=data_dict_output)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
DataClasses().ACEII_file_executor(L2_to_DistFunc,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)
