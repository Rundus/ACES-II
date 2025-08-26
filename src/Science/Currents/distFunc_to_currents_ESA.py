# --- distFunc_to_currents_ESA_old.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert electrostatic analyzer data from Distribution Function to Current using
# the first moment of the distribution funciton


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

# --- --- --- --- ---
from src.my_imports import *

start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = []

outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.integrate import simpson
from tqdm import tqdm

def DistFunc_to_ESAcurrents(wRocket, justPrintFileNames,wFile):
    # Set the paths for the file names
    rocketFolderPath = DataPaths.ACES_data_folder
    inputFiles = [f for f in glob(f'{rocketFolderPath}\\L3\\DistFunc\\{ACESII.fliers[wRocket - 4]}\*.cdf*') if "eepaa" in f or "iepaa" in f or "leesa" in f]
    names = ['eepaa','iepaa','leesa']
    inputFiles_Lshell = glob(f'{rocketFolderPath}\\coordinates\\Lshell\\{ACESII.fliers[wRocket - 4]}\*.cdf*')

    if justPrintFileNames:
        for idx,file in enumerate(inputFiles):
            print(f'[{idx}] ', file)
        return

    print('\n')
    print(stl.color.UNDERLINE + f'Calculating J_parallel from {names[wFile]} data' + stl.color.END)

    # --- Load the Data ---
    stl.prgMsg('Loading data from distribution function Files')
    data_dict_ESA = stl.loadDictFromFile(inputFiles[wFile])
    data_dict_LShell = stl.loadDictFromFile(inputFiles_Lshell[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}

    # --- --- --- --- --- --- --- --- --- --- -
    # --- ELECTRON - CALCULATE ESA CURRENTS ---
    # --- --- --- --- --- --- --- --- --- --- -
    stl.prgMsg('Calculating Parallel Current\n')
    distFunc = deepcopy(data_dict_ESA['Distribution_Function'][0])
    instr_nam = names[wFile]

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
    # TODO: Adjust pitch_integrates for IEPAA and LEESA
    elec_pitches = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180]
    ion_pitches = [0, 30, 60, 90, 120, 150, 180]
    if wRocket-4 == 0:
        if names[wFile] != 'iepaa':
            pitch_integrates = elec_pitches[1:17 + 1]
            integrand = integrand[:, 2:2 + len(pitch_integrates), :]
        else:
            pitch_integrates = ion_pitches[1:]
            integrand = integrand[:, 1:1 + len(pitch_integrates), :]
    elif wRocket-4 == 1:
        if names[wFile] != 'iepaa':
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

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:
        stl.prgMsg('Creating output file')

        fileOutPath = f'ACESII_{ACESII.payload_IDs[wRocket-4]}_Jpara_{names[wFile]}.cdf'
        outputPath = f'C:\Data\ACESII\science\ESA_currents\\{ACESII.fliers[wRocket-4]}\\{fileOutPath}'

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

        stl.outputCDFdata(data_dict=data_dict_output,outputPath=outputPath)

        stl.Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

if justPrintFileNames:
    DistFunc_to_ESAcurrents(wRocket, justPrintFileNames, 0)
else:
    if wFiles == []:
        for fileNo in range(len(glob(f'{DataPaths.ACES_data_folder}\\L3\\DistFunc\\{ACESII.fliers[wRocket - 4]}\*.cdf*'))):
            DistFunc_to_ESAcurrents(wRocket, justPrintFileNames, fileNo)
    else:
        for fileNo in wFiles:
            DistFunc_to_ESAcurrents(wRocket, justPrintFileNames, fileNo)
