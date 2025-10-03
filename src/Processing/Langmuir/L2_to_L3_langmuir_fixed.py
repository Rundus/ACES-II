# --- L2_to_L3_langmuir_fixed.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Takes in the L2 Langmuir Data and
# (1) Converts fixed Fixed LP current to ni
# (2) performs the Chi Square analysis to get the Temperature and Density from the
# characteristic curves of the swept probe


# NOTE: to convert to plasma density, the average ion mass and temperature was assumed to be
# Ti_assumed = 1 eV OR uses the smoothed IRI model to estimate
# m_i = average mass from the smoothed IRI model


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
from src.my_imports import *
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
wRocket = 5

# --- OutputData ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.Processing.Langmuir.toggles import L2toL3_FixedLPToggles as fToggles

#######################
# --- MAIN FUNCTION ---
#######################
def L2_to_L3_langmuir_fixed(wRocket, justPrintFileNames):

    # --- FILE I/O ---
    rocket_folder_path = DataPaths.ACES_data_folder
    rocketID = ACESII.payload_IDs[wRocket-4]
    globalAttrsMod = deepcopy(ACESII.global_attributes[wRocket-4])
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'LangmuirData'

    inputFiles = glob(f'{rocket_folder_path}{fToggles.inputPath_modifier}\{ACESII.fliers[wRocket-4]}{fToggles.modifier}\*langmuir_fixed*')
    input_names = [ifile.replace(f'{rocket_folder_path}{fToggles.inputPath_modifier}\{ACESII.fliers[wRocket-4]}{fToggles.modifier}\\', '') for ifile in inputFiles]
    input_names_searchable = [ifile.replace(fToggles.inputPath_modifier.lower() + '_', '').replace('_v00', '') for ifile in input_names]
    dataFile_name = inputFiles[0].replace(f'{rocket_folder_path}{fToggles.inputPath_modifier}\{ACESII.fliers[wRocket-4]}{fToggles.modifier}\\','')

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(os.path.getsize(file) / (10 ** 6), 1)))
        return

    print('\n')
    print(stl.color.UNDERLINE + f'Converting to {fToggles.outputPath_modifier} data for {dataFile_name}' + stl.color.END)
    print('[' + str(0) + ']   ' + str(round(os.path.getsize(inputFiles[0]) / (10 ** 6), 1)) + 'MiB')

    # --- get the data from the L2 file ---
    stl.prgMsg(f'Loading data')
    data_dict = stl.loadDictFromFile(inputFiles[0])

    # load the Langmuir Probe Post-flight Calibration data
    data_dict_LPpostFlight = stl.loadDictFromFile(glob(rf'C:\Data\ACESII\calibration\LP_postFlight_calibration\\{ACESII.fliers[wRocket-4]}\\*_postFlight_cal*')[0])
    stl.Done(start_time)
    ##################################################################
    # --- Calculate the plasma density from Ion Saturation Current ---
    ##################################################################

    stl.prgMsg('Calculating Fixed ni')
    fixed_current = data_dict['fixed_current'][0]  # note: current is in Nano Amps
    Te_eV = data_dict_LPpostFlight['Te_DERPA2'][0]
    Ti = Te_eV/data_dict_LPpostFlight['Tr'][0]
    vth_i = np.array(np.sqrt((stl.q0*Ti)/ (data_dict_LPpostFlight['m_eff_i'][0]*2*np.pi)))
    Phi_plas_minus_absphi_floating_assumed = 1 # ASSUMES a plasma potential 1V above the floating potential
    Phi_probe_minus_phi_plas = -1*(np.abs(data_dict_LPpostFlight['floating_potential'][0]) + Phi_plas_minus_absphi_floating_assumed + np.abs(ACESII.LP_fixed_probe_bias[wRocket-4]))

    ni_density_m3 = (-1*fixed_current*(1E-9)) / (stl.q0*ACESII.LP_probe_areas[wRocket-4][0]*vth_i*(1 + Phi_probe_minus_phi_plas/Ti))

    # create an output data_dict for the fixed data
    data_dict_fixed = {'ni': [ni_density_m3, {'LABLAXIS': 'plasma density', 'DEPEND_0': 'Epoch','UNITS': '!Nm!A-3!N', 'VAR_TYPE': 'data'}],
                       'Epoch': deepcopy(data_dict['Epoch']),
                       'fixed_boom_monitor': deepcopy(data_dict['fixed_boom_monitor']),
                       'Epoch_boom_monitor': deepcopy(data_dict['Epoch_boom_monitor'])
                       }

    if outputData:

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        stl.prgMsg('Creating output file')
        fileoutName_fixed = f'ACESII_{rocketID}_l3_langmuir_fixed.cdf'
        outputPath = f'{rocket_folder_path}{fToggles.outputPath_modifier}\{ACESII.fliers[wRocket-4]}\\{fileoutName_fixed}'
        stl.outputCDFdata(outputPath, data_dict_fixed, instrNam= 'Langmuir Probe',globalAttrsMod=globalAttrsMod)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
data_directory = f'{DataPaths.ACES_data_folder}{fToggles.inputPath_modifier}\{ACESII.fliers[wRocket-4]}{fToggles.modifier}\*.cdf'

if len(glob(data_directory)) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    L2_to_L3_langmuir_fixed(wRocket, justPrintFileNames)

