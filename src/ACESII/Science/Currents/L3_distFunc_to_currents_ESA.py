# --- distFunc_to_currents_ESA_old.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert electrostatic analyzer data from Distribution Function to Current using
# the first moment of the distribution funciton


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
# --- --- --- --- ---

######################
# --- DATA TOGGLES ---
######################
just_print_file_names_bool = False
rocket_str = 'low'
# wInstr = 'EEPAA'
wInstr = 'IEPAA'
# wInstr = 'LEESA'
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['L3', [[1],[0]]],
    'Lshell':['coordinates',[[0],[0]]]
}
outputData = True


#################
# --- IMPORTS ---
#################
from scipy.integrate import simpson
from src.ACESII.data_tools.my_imports import *
start_time = time.time()


def DistFunc_to_ESAcurrents(data_dicts):


    # --- Load the Data ---
    stl.prgMsg('Loading data from distribution function Files')
    data_dict_ESA = deepcopy(data_dicts[0])
    data_dict_LShell = deepcopy(data_dicts[1])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}

    # --- --- --- --- --- --- --- --- --- --- -
    # --- ELECTRON - CALCULATE ESA CURRENTS ---
    # --- --- --- --- --- --- --- --- --- --- -
    stl.prgMsg('Calculating Parallel Current/n')
    distFunc = deepcopy(data_dict_ESA['Distribution_Function'][0])
    instr_nam = wInstr

    # [0] Define the charge
    charge = -1*stl.q0 if instr_nam in ['eepaa','leesa'] else stl.q0

    # [1] Define the mass
    mass = stl.m_e if instr_nam in ['eepaa','leesa'] else stl.ion_dict['O+']

    # [2] Convert the energies to joules
    Energy_Joules = np.array(deepcopy(data_dict_ESA['Energy'][0])*stl.q0)

    # [3] Prepare the integrand

    # [3a] Multiply by the energy
    integrand = np.multiply(distFunc,Energy_Joules)

    # [3b] Multiply the integrand by pitch angle dependence
    for idx, ptch in enumerate(data_dict_ESA['Pitch_Angle'][0]):
        integrand[:,idx,:] = integrand[:, idx, :]*np.cos(np.radians(ptch))*np.cos(np.radians(ptch))

    # [4] Integrate over energy
    elec_pitches = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180]
    ion_pitches = [0, 30, 60, 90, 120, 150, 180]
    if rocket_str=='high':
        if wInstr in ['EEPAA','LEESA']:
            pitch_integrates = elec_pitches[1:17 + 1]
            integrand = integrand[:, 2:2 + len(pitch_integrates), :]
        else:
            pitch_integrates = ion_pitches[1:]
            integrand = integrand[:, 1:1 + len(pitch_integrates), :]
    elif rocket_str=='low':
        if wInstr in ['EEPAA','LEESA']:
            pitch_integrates = elec_pitches
            integrand = integrand[:, 1:1 + len(pitch_integrates), :]
        else:
            pitch_integrates = ion_pitches

    J_para_ptch = np.zeros(shape=(len(data_dict_ESA['Epoch'][0]),len(pitch_integrates)))

    for tme in tqdm(range(len(J_para_ptch))):
        for ptch in range(len(pitch_integrates)):
            J_para_ptch[tme][ptch] = simpson(integrand[tme][ptch], x=Energy_Joules)

    # [4] Integrate over Pitch Angle - NOTE: ONLY integrate over 0 to 180deg
    J_para = np.zeros(shape=len(data_dict_ESA['Epoch'][0]))
    for tme in tqdm(range(len(J_para_ptch))):
        J_para[tme] = simpson(J_para_ptch[tme], x=np.radians(pitch_integrates))

    # [5] Multiply by the constant term
    # TODO: Adjust mass for IEPAA
    J_para = (4*np.pi*charge/np.power(mass,2)) * J_para

    # [6] Remove significant outliers
    J_para[np.abs(J_para)>200E-6] = 0

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:
        stl.prgMsg('Creating output file')

        fileOutPath = f'ACESII_{ACESII.fliers_dict[rocket_str]}_Jpara_{wInstr}.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/science/ESA_currents/{wInstr}/{rocket_str}//{fileOutPath}'

        # output the main data
        data_dict_output = {**data_dict_output, **{'j_para':
                                         [J_para, {'LABLAXIS': 'J_parallel',
                                                   'DEPEND_0': 'Epoch',
                                                   'DEPEND_1': None,
                                                   'DEPEND_2': None,
                                                   'UNITS': '!N A!N m!U-2!N',
                                                   'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}

        # Add the other data
        data_dict_output = {**data_dict_output,
                            **{f'{key}':deepcopy(data_dict_ESA[f'{key}']) for key in ['Epoch','Alt','Lat','Long']}}

        # Include the L-Shell
        data_dict_output = {**data_dict_output,
                            **{f'{key}': deepcopy(data_dict_LShell[f'{key}']) for key in ['L-Shell']}}

        data_dict_output['L-Shell'][1]['VAR_TYPE'] = 'support_data'

        stl.outputDataDict(data_dict=data_dict_output,outputPath=outputPath)

        stl.Done(start_time)


#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(DistFunc_to_ESAcurrents,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)
