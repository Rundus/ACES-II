# --- diffFlux_to_Energy_Flux.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: using the specs of the ACESII ESAs, convert from differential Energy Flux
# to just energy flux as described in EEPAA_Flux_Conversion.pdf document in Overleaf



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

from src.my_imports import *
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = True
wRocket = 4
wFiles = [0]
modifier = ''
inputPath_modifier = 'L2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'L3\Energy_Flux' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# ---------------------------
outputData = False
# ---------------------------
erg_to_eV = 6.242E11 # eV per erg




# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.integrate import simpson
import spaceToolsLib as stl
from tqdm import tqdm


def diffFlux_to_Energy_Flux(wRocket, rocketFolderPath, justPrintFileNames, wflyer, wfile):


    # --- FILE I/O ---
    rocket_folder_path = DataPaths.ACES_data_folder
    rocketID = ACESII.payload_IDs[wflyer]
    globalAttrsMod = deepcopy(ACESII.global_attributes[wflyer])
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L3'

    input_files = glob(f'{rocketFolderPath}{inputPath_modifier}\{DataPaths.fliers[wflyer]}{modifier}\*.cdf')
    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{DataPaths.fliers[wflyer]}{modifier}\\', '') for ifile in input_files]
    input_names_searchable = [ifile.replace(inputPath_modifier.lower() +'_', '') for ifile in input_names]

    if justPrintFileNames:
        for i, file in enumerate(input_files):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(os.path.getsize(file) / (10 ** 6), 1)))
        return

    # --- --- --- --- --- -
    # --- LOAD THE DATA ---
    # --- --- --- --- --- -
    # --- get the data from the file ---
    stl.prgMsg(f'Loading data from ESA Files')
    data_dict, globalAttrs = stl.loadDictFromFile(inputFilePath=input_files[wfile], getGlobalAttrs=True)
    stl.Done(start_time)

    # --- --- --- --- -
    # --- INTEGRATE ---
    # --- --- --- --- -
    stl.prgMsg('Calculating Fluxes')
    Epoch = data_dict['Epoch'][0]
    Pitch = data_dict['Pitch_Angle'][0][1:20] # only get pitch angles 0deg to 180deg
    Energy = data_dict['Energy'][0]
    diffEFlux = data_dict['Differential_Energy_Flux'][0][:,1:20,:] # ONLY get 0deg to 180deg
    diffNFlux = data_dict['Differential_Number_Flux'][0][:,1:20,:] # ONLY get 0deg to 180deg

    downwardPitchRange = [0, 9 + 1]  # what pitch indicies to consider when calculating parallel (downward) # 0-90deg (IDX: 6)
    upwardPitchRange = [11, 19 + 1]  # 90-180deg

    Phi_N = np.zeros(shape=(len(Epoch)))
    Phi_N_antiParallel = np.zeros(shape=(len(Epoch)))
    Phi_N_Parallel = np.zeros(shape=(len(Epoch)))
    varPhi_N = np.zeros(shape=(len(Epoch), len(Energy)))
    varPhi_N_antiParallel = np.zeros(shape=(len(Epoch), len(Energy)))
    varPhi_N_Parallel = np.zeros(shape=(len(Epoch), len(Energy)))

    # --- perform the integration ---
    steradian = 2 * np.pi * np.array([np.sin(np.radians(ptch)) for ptch in Pitch])

    for tme in tqdm(range(len(Epoch))):

        # calculate the partial differential flux
        prepared = np.transpose(diffNFlux[tme]*steradian)
        varPhi_N[tme] = np.array([simpson(y=prepared[idx], x=Pitch) for idx, engy in enumerate(Energy)])
        varPhi_N_antiParallel[tme] = np.array([simpson(y=prepared[idx][0:10], x=Pitch[0:10]) for idx, engy in enumerate(Energy)])
        varPhi_N_Parallel[tme] = np.array([simpson(y=prepared[idx][10:20], x=Pitch[10:20]) for idx, engy in enumerate(Energy)])

        Phi_N[tme] = simpson(y=varPhi_N[tme], x=Energy)
        Phi_N_antiParallel[tme] = simpson(y=varPhi_N_antiParallel[tme], x=Energy)
        Phi_N_Parallel[tme] = simpson(y=varPhi_N_Parallel[tme], x=Energy)

    stl.Done(start_time)



    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')

        data_dict_output = {'Phi_N': [Phi_N, deepcopy(data_dict['Differential_Energy_Flux'][1])],
                            'Phi_N_antiParallel': [Phi_N_antiParallel, deepcopy(data_dict['Differential_Energy_Flux'][1])],
                            'Phi_N_Parallel': [Phi_N_Parallel, deepcopy(data_dict['Differential_Energy_Flux'][1])],
                            'Pitch_Angle': data_dict['Pitch_Angle'],
                            'Energy': data_dict['Energy'],
                            'Epoch': data_dict['Epoch']
                            }

        data_dict_output['Phi_N'][1]['LABLAXIS'] = 'Number_Flux'
        data_dict_output['Phi_N'][1]['UNITS'] = 'cm!A-2!N s!A-1!N'

        data_dict_output['Phi_N_antiParallel'][1]['LABLAXIS'] = 'Anti_Parallel_Number_Flux'
        data_dict_output['Phi_N_antiParallel'][1]['UNITS'] = 'cm!A-2!N s!A-1!N'

        data_dict_output['Phi_N_Parallel'][1]['LABLAXIS'] = 'Parallel_Number_Flux'
        data_dict_output['Phi_N_Parallel'][1]['UNITS'] = 'cm!A-2!N s!A-1!N'

        # write out the data
        fileoutName = f'ACESII_{rocketID}_eepaa_Energy_Flux.cdf'
        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output, globalAttrsMod=globalAttrs, instrNam=wInstr[1])
        stl.Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocket_folder_path = DataPaths.ACES_data_folder
if len(glob(f'{rocket_folder_path}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}\*.cdf')) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    if justPrintFileNames:
        diffFlux_to_Energy_Flux(wRocket, rocket_folder_path, justPrintFileNames, wRocket-4, 0)
    elif wFiles == []:
        for i, wfile in enumerate(glob(f'{rocket_folder_path}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}\*.cdf')):
            diffFlux_to_Energy_Flux(wRocket, rocket_folder_path, justPrintFileNames, wfile)
    else:
        for wfile in wFiles:
            diffFlux_to_Energy_Flux(wRocket, rocket_folder_path, justPrintFileNames, wfile)
