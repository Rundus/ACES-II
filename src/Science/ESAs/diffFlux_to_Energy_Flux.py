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
justPrintFileNames = False
wRocket = 5
wFiles = [[0], [1]]
modifier = ''
inputPath_modifier = 'L2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder

# ---------------------------
outputData = True
# ---------------------------



# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# from scipy.integrate import trapz
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
    data_dict_eepaa, globalAttrs = stl.loadDictFromFile(inputFilePath=input_files[wfile], getGlobalAttrs=True)
    data_dict_Lshell = stl.loadDictFromFile(glob(f'C:\Data\ACESII\coordinates\Lshell\{ACESII.fliers[wRocket-4]}\\*.cdf*')[0])
    stl.Done(start_time)

    # --- --- --- --- -
    # --- INTEGRATE ---
    # --- --- --- --- -
    stl.prgMsg('Calculating Fluxes')
    Epoch = data_dict_eepaa['Epoch'][0]
    Pitch = data_dict_eepaa['Pitch_Angle'][0][1:20] # only get pitch angles 0deg to 180deg
    Energy = data_dict_eepaa['Energy'][0]
    diffNFlux = data_dict_eepaa['Differential_Number_Flux'][0][:, 1:20, :] # ONLY get 0deg to 180deg

    # Number Fluxes
    Phi_N = np.zeros(shape=(len(Epoch)))
    Phi_N_antiParallel = np.zeros(shape=(len(Epoch)))
    Phi_N_Parallel = np.zeros(shape=(len(Epoch)))
    varPhi_N = np.zeros(shape=(len(Epoch), len(Energy)))
    varPhi_N_antiParallel = np.zeros(shape=(len(Epoch), len(Energy)))
    varPhi_N_Parallel = np.zeros(shape=(len(Epoch), len(Energy)))

    # Energy Fluxes
    Phi_E = np.zeros(shape=(len(Epoch)))
    Phi_E_antiParallel = np.zeros(shape=(len(Epoch)))
    Phi_E_Parallel = np.zeros(shape=(len(Epoch)))
    varPhi_E = np.zeros(shape=(len(Epoch), len(Energy)))
    varPhi_E_antiParallel = np.zeros(shape=(len(Epoch), len(Energy)))
    varPhi_E_Parallel = np.zeros(shape=(len(Epoch), len(Energy)))

    # Average Energy
    Energy_avg_Robinson = np.zeros(shape=(len(Epoch)))
    cutoff_idx = np.abs(Energy - 100).argmin()

    # determine the DeltaE to use for the varphi integrations - DeltaE = the half distance to the next energy value
    deltaEs = []
    for idx, engy in enumerate(Energy):
        if idx == len(Energy)-1:
            deltaEs.append(Energy[-2] - Energy[-1])
        elif idx == 0:
            deltaEs.append(Energy[0] - Energy[1])
        else:
            lowerE =(Energy[idx] - Energy[idx+1])/2
            highE = (Energy[idx-1] - Energy[idx])/2
            deltaEs.append(lowerE+highE)


    # --- perform the integration_tad_files ---
    alpha_grid, engy_grid = np.meshgrid(Pitch, data_dict_eepaa['Energy'][0])

    for tme in tqdm(range(len(Epoch))):

        #######################
        # --- NUMBER FLUXES ---
        #######################
        # --- Integrate over pitch angle ---

        # omni directional
        prepared = np.transpose(diffNFlux[tme]*2*np.pi*np.sin(np.radians(alpha_grid.T)))
        JE_N = np.array([simpson(y=prepared[idx], x=np.radians(Pitch)) for idx, engy in enumerate(Energy)]) # J_N(E) diffNFlux without pitch

        # parallel
        prepared = np.transpose(diffNFlux[tme]*2*np.pi*np.sin(np.radians(alpha_grid.T))*np.cos(np.radians(alpha_grid.T)))
        JE_N_Parallel = np.array([simpson(y=prepared[idx][0:10], x=np.radians(Pitch[0:10])) for idx, engy in enumerate(Energy)])

        # anti-parallel
        prepared = np.transpose(diffNFlux[tme]*2*np.pi* np.sin(np.radians(alpha_grid.T)) * np.cos(np.radians(alpha_grid.T)))
        JE_N_antiParallel = np.array([simpson(y=prepared[idx][10:], x=np.radians(Pitch[10:])) for idx, engy in enumerate(Energy)])

        # --- partially integrate over energy. ---
        # Description: To get in units of [cm^-2 s^-1], we assume j(E) doesn't change over the DeltaE interval between samples.
        # The integral between E-DeltaE and E+DeltaE around a central energy E is just: varphi(E) = DeltaE(E) * J(E) where DeltaE(E) depends
        # on the central energy. In our detector, DeltaE(E) is designed to be ~18% always --> DeltaE(E) = (1+gamma)E -(1-gamma)E = 2*gamma*E
        varPhi_N[tme] = np.array([deltaEs[idx]*JE_N[idx] for idx, engy in enumerate(Energy)])
        varPhi_N_Parallel[tme] = np.array([deltaEs[idx]*JE_N_Parallel[idx] for idx, engy in enumerate(Energy)])
        varPhi_N_antiParallel[tme] = np.array([deltaEs[idx]*JE_N_antiParallel[idx] for idx, engy in enumerate(Energy)])

        # Integrate over energy
        Phi_N[tme] = np.array(-1*simpson(y=JE_N, x=Energy)).clip(min=0)
        Phi_N_antiParallel[tme] = np.array(-1*simpson(y=JE_N_antiParallel, x=Energy)).clip(min=0)
        Phi_N_Parallel[tme] = np.array(-1*simpson(y=JE_N_Parallel, x=Energy)).clip(min=0)

        # ---------------------
        # --- ENERGY FLUXES ---
        # ---------------------
        varPhi_E[tme] = Energy*deepcopy(varPhi_N[tme])
        varPhi_E_antiParallel[tme] = Energy * deepcopy(varPhi_N_antiParallel[tme])
        varPhi_E_Parallel[tme] = Energy * deepcopy(varPhi_N_Parallel[tme])

        # Integrate over energy
        Phi_E[tme] = np.array(-1*simpson(y=JE_N*Energy, x=Energy)).clip(min=0)
        Phi_E_antiParallel[tme] = np.array(-1*simpson(y=JE_N_antiParallel*Energy, x=Energy)).clip(min=0)
        Phi_E_Parallel[tme] = np.array(-1*simpson(y=JE_N_Parallel*Energy, x=Energy)).clip(min=0)

        # -------------------------
        # --- ROBINSON FORMULAE ---
        # -------------------------

        # calculate the average energy
        # NOTE: the low-energy electrons DON'T contribute to the height_integrated conductivity very much,
        # thus ONLY use the PRIMARY beam to calculate the average energy, which is ~ 116 eV
        Energy_avg_Robinson[tme] = np.array(-1*simpson(y=JE_N_Parallel[:cutoff_idx]*Energy[:cutoff_idx], x=Energy[:cutoff_idx]))/np.array(-1*simpson(y=JE_N_Parallel[:cutoff_idx], x=Energy[:cutoff_idx]))


        # -------------------------
        # --- KAEPPLER FORMULAE ---
        # -------------------------

    stl.Done(start_time)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')

        data_dict_output = {'Phi_N': [Phi_N, deepcopy(data_dict_eepaa['Differential_Energy_Flux'][1])],
                            'Phi_N_antiParallel': [Phi_N_antiParallel, deepcopy(data_dict_eepaa['Differential_Energy_Flux'][1])],
                            'Phi_N_Parallel': [Phi_N_Parallel, deepcopy(data_dict_eepaa['Differential_Energy_Flux'][1])],
                            'varPhi_N': [varPhi_N, deepcopy(data_dict_eepaa['Differential_Energy_Flux'][1])],
                            'varPhi_N_antiParallel': [varPhi_N_antiParallel, deepcopy(data_dict_eepaa['Differential_Energy_Flux'][1])],
                            'varPhi_N_Parallel': [varPhi_N_Parallel, deepcopy(data_dict_eepaa['Differential_Energy_Flux'][1])],

                            'Phi_E': [Phi_E, deepcopy(data_dict_eepaa['Differential_Energy_Flux'][1])],
                            'Phi_E_antiParallel': [Phi_E_antiParallel, deepcopy(data_dict_eepaa['Differential_Energy_Flux'][1])],
                            'Phi_E_Parallel': [Phi_E_Parallel, deepcopy(data_dict_eepaa['Differential_Energy_Flux'][1])],
                            'varPhi_E': [varPhi_E, deepcopy(data_dict_eepaa['Differential_Energy_Flux'][1])],
                            'varPhi_E_antiParallel': [varPhi_E_antiParallel,deepcopy(data_dict_eepaa['Differential_Energy_Flux'][1])],
                            'varPhi_E_Parallel': [varPhi_E_Parallel,deepcopy(data_dict_eepaa['Differential_Energy_Flux'][1])],

                            'Pitch_Angle': data_dict_eepaa['Pitch_Angle'],
                            'Energy': data_dict_eepaa['Energy'],
                            'Epoch': data_dict_eepaa['Epoch'],
                            'Alt': data_dict_eepaa['Alt'],
                            'L-Shell':data_dict_Lshell['L-Shell'],

                            'Energy_avg': [Energy_avg_Robinson, data_dict_eepaa['Energy'][1]],
                            }

        data_dict_output['Phi_N'][1]['LABLAXIS'] = 'Number_Flux'
        data_dict_output['Phi_N'][1]['UNITS'] = 'cm!A-2!N s!A-1!N'

        data_dict_output['Phi_E'][1]['LABLAXIS'] = 'Energy_Flux'
        data_dict_output['Phi_E'][1]['UNITS'] = 'eV cm!A-2!N s!A-1!N'

        data_dict_output['varPhi_N'][1]['UNITS'] = 'cm!A-2!N s!A-1!N'
        data_dict_output['varPhi_N'][1]['DEPEND_1'] = 'Energy'
        data_dict_output['varPhi_N'][1]['DEPEND_2'] = None

        data_dict_output['varPhi_E'][1]['UNITS'] = 'eV cm!A-2!N s!A-1!N'
        data_dict_output['varPhi_E'][1]['DEPEND_1'] = 'Energy'
        data_dict_output['varPhi_E'][1]['DEPEND_2'] = None


        data_dict_output['Phi_N_antiParallel'][1]['LABLAXIS'] = 'Anti_Parallel_Number_Flux'
        data_dict_output['Phi_N_antiParallel'][1]['UNITS'] = 'cm!A-2!N s!A-1!N'

        data_dict_output['Phi_E_antiParallel'][1]['LABLAXIS'] = 'Anti_Parallel_Energy_Flux'
        data_dict_output['Phi_E_antiParallel'][1]['UNITS'] = 'eV cm!A-2!N s!A-1!N'

        data_dict_output['Phi_N_Parallel'][1]['LABLAXIS'] = 'Parallel_Number_Flux'
        data_dict_output['Phi_N_Parallel'][1]['UNITS'] = 'cm!A-2!N s!A-1!N'

        data_dict_output['Phi_E_Parallel'][1]['LABLAXIS'] = 'Parallel_Energy_Flux'
        data_dict_output['Phi_E_Parallel'][1]['UNITS'] = 'eV cm!A-2!N s!A-1!N'


        data_dict_output['varPhi_N_antiParallel'][1]['UNITS'] = 'cm!A-2!N s!A-1!N'
        data_dict_output['varPhi_N_antiParallel'][1]['DEPEND_1'] = 'Energy'
        data_dict_output['varPhi_N_antiParallel'][1]['DEPEND_2'] = None

        data_dict_output['varPhi_E_antiParallel'][1]['UNITS'] = 'eV cm!A-2!N s!A-1!N'
        data_dict_output['varPhi_E_antiParallel'][1]['DEPEND_1'] = 'Energy'
        data_dict_output['varPhi_E_antiParallel'][1]['DEPEND_2'] = None

        data_dict_output['varPhi_N_Parallel'][1]['UNITS'] = 'cm!A-2!N s!A-1!N'
        data_dict_output['varPhi_N_Parallel'][1]['DEPEND_1'] = 'Energy'
        data_dict_output['varPhi_N_Parallel'][1]['DEPEND_2'] = None

        data_dict_output['varPhi_E_Parallel'][1]['UNITS'] = 'eV cm!A-2!N s!A-1!N'
        data_dict_output['varPhi_E_Parallel'][1]['DEPEND_1'] = 'Energy'
        data_dict_output['varPhi_E_Parallel'][1]['DEPEND_2'] = None

        # write out the data
        fileoutName = f'ACESII_{rocketID}_l3_eepaa_flux.cdf'
        outputPath = f'{rocketFolderPath}\\L3\Energy_Flux\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output, globalAttrsMod=globalAttrs)
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
    elif wFiles[wRocket-4] == []:
        for i, wfile in enumerate(glob(f'{rocket_folder_path}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}\*.cdf')):
            diffFlux_to_Energy_Flux(wRocket, rocket_folder_path, justPrintFileNames, wRocket-4, wfile)
    else:
        for wfile in wFiles[wRocket-4]:
            diffFlux_to_Energy_Flux(wRocket, rocket_folder_path, justPrintFileNames, wRocket-4,wfile)
