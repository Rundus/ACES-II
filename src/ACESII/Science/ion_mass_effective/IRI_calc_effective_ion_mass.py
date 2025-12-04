# --- IRI_calc_effective_ion_mass.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Use the IRI model (or my ionosphere simluation) to calculate the
# likely effective ion mass vs. altitude for the two payloads on the EEPAA timebase

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
from src.ACESII.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False # Just print the names of files

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4

# --- OutputData ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---

# --- OutputData ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.interpolate import CubicSpline

#######################
# --- MAIN FUNCTION ---
#######################
def IRI_calc_effective_ion_mass(wRocket, justPrintFileNames):

    # --- FILE I/O ---
    input_files_alt = glob(rf'C:\Data\ACESII\attitude\{ACESII.fliers[wRocket-4]}\\*.cdf*')
    input_files_simIRI = glob(rf'C:\Data\physicsModels\ionosphere\plasma_environment\\*.cdf*')
    input_file_ESA = glob(rf'C:\Data\ACESII\L2\{ACESII.fliers[wRocket-4]}\\*eepaa_fullCal.cdf*')

    if justPrintFileNames:
        for i, file in enumerate(input_files_alt):
            print(file)
        for i, file in enumerate(input_files_simIRI):
            print(file)
        return

    # --- get the data from the L2 file ---
    stl.prgMsg(f'Loading data')
    data_dict_ESA = stl.loadDictFromFile(input_file_ESA[0])
    data_dict_simIRI = stl.loadDictFromFile(input_files_simIRI[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {'Epoch':deepcopy(data_dict_ESA['Epoch'])}

    ######################
    # --- INTERPOLATE  ---
    ######################
    stl.prgMsg('Interpolating IRI m_eff_i')
    for key in ['m_eff_i', 'C_O+', 'C_H+', 'C_NO+', 'Te', 'Ti', 'n_NO+', 'n_O2+', 'n_O+']:
        cs = CubicSpline(data_dict_simIRI['simAlt'][0], data_dict_simIRI[f'{key}'][0][0])
        data_dict_output = {**data_dict_output,
                            **{f'{key}':[np.array(cs(data_dict_ESA['Alt'][0])), deepcopy(data_dict_simIRI[f'{key}'][1])]}
                            }
        data_dict_output[f'{key}'][1]['DEPEND_0'] = 'Epoch'

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:
        stl.prgMsg('Creating output file')
        fileoutName = f'ACESII_{ACESII.payload_IDs[wRocket-4]}_effective_ion_mass.cdf'
        outputPath = f'C:\Data\ACESII\science\ion_mass_effective\\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(data_dict=data_dict_output, outputPath=outputPath)
        stl.Done(start_time)


IRI_calc_effective_ion_mass(wRocket, justPrintFileNames)
