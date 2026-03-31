# --- L2_to_L3_energy_flux.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: using the specs of the ACESII ESAs, convert from differential Energy Flux
# to just energy flux as described in EEPAA_Flux_Conversion.pdf document in Overleaf

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
wInstr = 'EEPAA'
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['L2', [[0],[0]]],
    'Lshell':['coordinates',[[0],[0]]]
}
outputData = True


# Robinson Formulae
char_energy_cutoff = 1000 # in [eV]. Should be the point where "the secondary electrons begin to dominate over the primary electrons. Robinson uses 500eV, but for our case that is clearly too low

#################
# --- IMPORTS ---
#################
from scipy.integrate import simpson
from src.ACESII.data_tools.my_imports import *
start_time = time.time()


def L2_to_L3_energy_flux(data_dicts):

    #######################
    # --- LOAD THE DATA ---
    #######################
    stl.prgMsg(f'Loading data')
    data_dict_ESA = deepcopy(data_dicts[0])
    data_dict_Lshell = deepcopy(data_dicts[1])
    stl.Done(start_time)

    ######################################
    # --- PREPARE THE OUTPUT VARIABLES ---
    ######################################
    stl.prgMsg('Calculating Fluxes')
    Epoch = data_dict_ESA['Epoch'][0]
    Pitch = data_dict_ESA['Pitch_Angle'][0] # only get pitch angles 0deg to 180deg
    Energy = data_dict_ESA['Energy'][0]
    diffNFlux = data_dict_ESA['Differential_Number_Flux'][0] # ONLY get 0deg to 180deg

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
    characteristic_energy_R87 = np.zeros(shape=(len(Epoch)))
    cutoff_idx = np.abs(Energy - char_energy_cutoff).argmin()

    data_dict_output = {'Phi_N': [Phi_N, deepcopy(data_dict_ESA['Differential_Number_Flux'][1])],
                       'Phi_N_antiParallel': [Phi_N_antiParallel, deepcopy(data_dict_ESA['Differential_Number_Flux'][1])],
                       'Phi_N_Parallel': [Phi_N_Parallel, deepcopy(data_dict_ESA['Differential_Number_Flux'][1])],
                       'varPhi_N': [varPhi_N, deepcopy(data_dict_ESA['Differential_Number_Flux'][1])],
                       'varPhi_N_antiParallel': [varPhi_N_antiParallel, deepcopy(data_dict_ESA['Differential_Number_Flux'][1])],
                       'varPhi_N_Parallel': [varPhi_N_Parallel, deepcopy(data_dict_ESA['Differential_Number_Flux'][1])],
                       'Phi_E': [Phi_E, deepcopy(data_dict_ESA['Differential_Energy_Flux'][1])],
                       'Phi_E_antiParallel': [Phi_E_antiParallel, deepcopy(data_dict_ESA['Differential_Energy_Flux'][1])],
                       'Phi_E_Parallel': [Phi_E_Parallel, deepcopy(data_dict_ESA['Differential_Energy_Flux'][1])],
                       'varPhi_E': [varPhi_E, deepcopy(data_dict_ESA['Differential_Energy_Flux'][1])],
                       'varPhi_E_antiParallel': [varPhi_E_antiParallel, deepcopy(data_dict_ESA['Differential_Energy_Flux'][1])],
                       'varPhi_E_Parallel': [varPhi_E_Parallel, deepcopy(data_dict_ESA['Differential_Energy_Flux'][1])],
                       'Pitch_Angle': data_dict_ESA['Pitch_Angle'],
                       'Energy': data_dict_ESA['Energy'],
                       'Epoch': data_dict_ESA['Epoch'],
                       'Alt': data_dict_ESA['Alt'],
                       'L-Shell': data_dict_Lshell['L-Shell'],
                       'characteristic_energy_R87': [characteristic_energy_R87, data_dict_ESA['Energy'][1]],
                       }

    data_dict_output['Phi_N'][1]['LABLAXIS'] = 'Number_Flux'
    data_dict_output['Phi_N'][1]['UNITS'] = 'cm!A-2!N s!A-1!N'

    data_dict_output['Phi_E'][1]['LABLAXIS'] = 'Energy_Flux'
    data_dict_output['Phi_E'][1]['UNITS'] = 'eV cm!A-2!N s!A-1!N'

    data_dict_output['varPhi_N'][1]['LABLAXIS'] = 'Number_Flux'
    data_dict_output['varPhi_N'][1]['UNITS'] = 'cm!A-2!N s!A-1!N'
    data_dict_output['varPhi_N'][1]['DEPEND_1'] = 'Energy'
    data_dict_output['varPhi_N'][1]['DEPEND_2'] = None

    data_dict_output['varPhi_N'][1]['LABLAXIS'] = 'Energy_Flux'
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


    ##################################
    # --- PERFORM THE INTEGRATIONS ---
    ##################################

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

    # For pitch angle integrations, determine the deltaPitch (in radians) for each anode
    para_idx = np.abs(Pitch - 0).argmin()
    perp_idx = np.abs(Pitch - 90).argmin()
    antipara_idx = np.abs(Pitch - 180).argmin()
    limits_radians = np.radians(ACESII.ESA_pitch_bin_integration_limits[wInstr][para_idx:antipara_idx+1])
    deltaPtch = [lim[1]-lim[0] for lim in limits_radians]
    f_x_factor = [2*np.sin(lim[0]/2 + lim[1]/2) * np.cos(lim[0]/2 - lim[1]/2) for lim in limits_radians]

    # Calculate the Fluxes using Trapezoidal Integration
    for tmeIdx in tqdm(range(len(Epoch))):

        data_slice = diffNFlux[tmeIdx,para_idx:antipara_idx+1,:]

        # --- NUMBER FLUX ---
        # trapezoidal integrate 0.5*(b-a)[f(a) + f(b)] away pitch angle dimension. Note: f(a) = diffNFlux[tme,ptch,engy]*sin(a), f(x) = diffNFlux[tme,ptch,engy]*sin(a/2 + b/2)*cos(a/2-b/2)

        # omni directional
        JE_N =  np.nansum(2*np.pi*0.5*np.array(data_slice).T*f_x_factor*deltaPtch,axis=1)

        # parallel flux
        JE_N_Parallel = np.nansum(2*np.pi*0.5*np.array(data_slice[:perp_idx+1,:]).T*f_x_factor[:perp_idx+1]*deltaPtch[:perp_idx+1],axis=1)

        # # anti-parallel flux
        JE_N_antiParallel = np.nansum(2*np.pi*0.5*np.array(data_slice[perp_idx:,:]).T*f_x_factor[perp_idx:]*deltaPtch[perp_idx:],axis=1)

        # --- partially integrate over energy. ---
        # Description: To get in units of [cm^-2 s^-1], we assume j(E) doesn't change over the DeltaE interval between samples.
        # The integral between E-DeltaE and E+DeltaE around a central energy E is just: varphi(E) = DeltaE(E) * J(E) where DeltaE(E) depends
        # on the central energy. In our detector, DeltaE(E) is designed to be ~18% always --> DeltaE(E) = (1+gamma)E -(1-gamma)E = 2*gamma*E
        varPhi_N[tmeIdx] = np.array(deltaEs)*JE_N
        varPhi_N_Parallel[tmeIdx] = deltaEs*JE_N_Parallel
        varPhi_N_antiParallel[tmeIdx] = deltaEs*JE_N_antiParallel

        # Properly Integrate over energy
        Phi_N[tmeIdx] = np.array(-1*simpson(y=JE_N, x=Energy)).clip(min=0)
        Phi_N_antiParallel[tmeIdx] = np.array(-1*simpson(y=JE_N_antiParallel, x=Energy)).clip(min=0)
        Phi_N_Parallel[tmeIdx] = np.array(-1*simpson(y=JE_N_Parallel, x=Energy)).clip(min=0)

        # ---------------------
        # --- ENERGY FLUXES ---
        # ---------------------
        varPhi_E[tmeIdx] = Energy*deepcopy(varPhi_N[tmeIdx])
        varPhi_E_antiParallel[tmeIdx] = Energy * deepcopy(varPhi_N_antiParallel[tmeIdx])
        varPhi_E_Parallel[tmeIdx] = Energy * deepcopy(varPhi_N_Parallel[tmeIdx])

        # Integrate over energy
        Phi_E[tmeIdx] = np.array(-1*simpson(y=JE_N*Energy, x=Energy)).clip(min=0)
        Phi_E_antiParallel[tmeIdx] = np.array(-1*simpson(y=JE_N_antiParallel*Energy, x=Energy)).clip(min=0)
        Phi_E_Parallel[tmeIdx] = np.array(-1*simpson(y=JE_N_Parallel*Energy, x=Energy)).clip(min=0)

        # -------------------------
        # --- ROBINSON FORMULAE ---
        # -------------------------

        # calculate the average energy
        # NOTE: the low-energy electrons DON'T contribute to the height_integrated conductivity very much,
        # thus ONLY use the PRIMARY beam to calculate the average energy, which is ~ 500 eV (according to robinson)
        characteristic_energy_R87[tmeIdx] = np.array(-1*simpson(y=JE_N_Parallel[:cutoff_idx]*Energy[:cutoff_idx], x=Energy[:cutoff_idx]))/np.array(-1*simpson(y=JE_N_Parallel[:cutoff_idx], x=Energy[:cutoff_idx]))

        # -------------------------
        # --- KAEPPLER FORMULAE ---
        # -------------------------


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
DataClasses().ACEII_file_executor(L2_to_L3_energy_flux,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)
