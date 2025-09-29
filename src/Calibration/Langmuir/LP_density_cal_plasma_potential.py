# --- LP_density_cal_plasma_potential.py ---
# use the EISCAT data to generate a nominal value of the plasma potential
# that can be used for the fitting routine later to get a more accurate
# value

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"

import numpy as np

from src.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False  # Just print the names of files

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4


showEISCAT_profiles = False


# --- OutputData ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from glob import glob
from scipy.signal import savgol_filter
from scipy.interpolate import CubicSpline
import spaceToolsLib as stl
import datetime as dt


#######################
# --- MAIN FUNCTION ---
#######################
def LP_density_cal_plasma_potential(wRocket):

    # --- FILE I/O ---
    stl.prgMsg('Loading Data')
    # load the rocket ion saturation current data
    data_path = glob(rf'C:\Data\ACESII\L2\{ACESII.fliers[wRocket-4]}\\*langmuir_fixed.cdf*')[0]
    data_dict_LP_current = stl.loadDictFromFile(data_path)

    # load the EISCAT data
    data_path = r'C:\Data\ACESII\science\EISCAT\tromso\UHF\MAD6400_2022-11-20_beata_ant@uhfa.cdf'
    data_dict_EISCAT = stl.loadDictFromFile(data_path)

    # load the attitude data
    data_path = glob(rf'C:\Data\ACESII\attitude\{ACESII.fliers[wRocket-4]}\\*.cdf*')[0]
    data_dict_attitude = stl.loadDictFromFile(data_path)

    # load the trajectory data
    data_path = glob(rf'C:\Data\ACESII\trajectories\{ACESII.fliers[wRocket - 4]}\\*_auroral.cdf*')[0]
    data_dict_traj = stl.loadDictFromFile(data_path)

    # --- prepare the output ---
    data_dict_output = {}
    stl.Done(start_time)


    # --- shorten the dataset ---
    low_idx = np.abs(data_dict_LP_current['Epoch'][0] - dt.datetime(2022,11,20,17,24)).argmin()
    high_idx = np.abs(data_dict_LP_current['Epoch'][0] - dt.datetime(2022, 11, 20, 17, 28)).argmin()
    data_dict_LP_current['fixed_current'][0]= data_dict_LP_current['fixed_current'][0][low_idx:high_idx+1]
    data_dict_LP_current['Epoch'][0] = data_dict_LP_current['Epoch'][0][low_idx:high_idx + 1]



    ###############################################
    # --- Collect the EISCAT Data and Smooth it ---
    ###############################################


    # Collect the ESICAT Profiles
    target_times = [dt.datetime(2022,11,20,17,21,25),
                    dt.datetime(2022,11,20,17,22,15),
                    dt.datetime(2022,11,20,17,25,25)]
    target_idxs = [np.abs(tme - data_dict_EISCAT['Epoch'][0]).argmin() for tme in target_times]
    target_Ti = data_dict_EISCAT['ti'][0][target_idxs]
    target_ni = data_dict_EISCAT['ne'][0][target_idxs]

    # Average the profiles
    # Ti
    Ti_profile = np.nanmean(target_Ti,axis=0) * (stl.kB/stl.q0)
    Ti_profile_idxs = np.isnan(Ti_profile)
    Ti_profile_alts = data_dict_EISCAT['range'][0][np.where(Ti_profile_idxs==False)[0]]
    Ti_profile_vals = Ti_profile[np.where(Ti_profile_idxs==False)[0]]

    # ni
    ni_profile = np.nanmean(target_ni, axis=0)
    ni_profile_idxs = np.isnan(ni_profile)
    ni_profile_alts = data_dict_EISCAT['range'][0][np.where(ni_profile_idxs==False)[0]]
    ni_profile_vals = ni_profile[np.where(ni_profile_idxs==False)[0]]

    # ---Smooth the profiles ---
    window_length = 20
    porder = 3
    Ti_profile_smooth = savgol_filter(x=Ti_profile_vals,window_length=window_length,polyorder=porder)
    Ti_profile_smoothed = savgol_filter(x=Ti_profile_smooth, window_length=window_length+20, polyorder=porder+2)
    ni_profile_smooth = savgol_filter(x=ni_profile_vals, window_length=window_length, polyorder=porder)
    ni_profile_smoothed = savgol_filter(x=ni_profile_smooth, window_length=window_length+20, polyorder=porder + 2)

    if showEISCAT_profiles:
        import matplotlib.pyplot as plt
        fig, ax =plt.subplots(nrows=1,ncols=2)

        ax[0].plot(Ti_profile_alts, Ti_profile_vals,color='tab:blue')
        ax[0].plot(Ti_profile_alts, Ti_profile_smoothed, color='tab:red')
        ax[0].set_xlabel('Alt [km]')
        ax[0].set_ylabel('Ti [eV]')

        ax[1].plot(ni_profile_alts, ni_profile_vals,color='tab:blue')
        ax[1].plot(ni_profile_alts, ni_profile_smoothed, color='tab:red')
        ax[1].set_yscale('log')
        ax[1].set_xlabel('Alt [km]')
        ax[1].set_ylabel('n [m^-3]')
        plt.show()

    # --- Interpolate the Ti and ni data onto the langmuir current timebase ---

    # create the altitude interpolation object
    T0 = dt.datetime(2022,11,20,17,20)
    T0_LP_current = stl.EpochTo_T0_Rocket(data_dict_LP_current['Epoch'][0], T0=T0)
    T0_attitude = stl.EpochTo_T0_Rocket(data_dict_attitude['Epoch'][0],T0=T0)
    cs1 = CubicSpline(T0_attitude, data_dict_attitude['Alt'][0]/stl.m_to_km)
    alt_langmuir = cs1(T0_LP_current)

    # Create the Ti profile for Langmuir
    cs2 = CubicSpline(Ti_profile_alts, Ti_profile_smoothed)
    Ti_langmuir = cs2(alt_langmuir)

    # Create the ni profile for langmuir
    cs3 = CubicSpline(ni_profile_alts, ni_profile_smoothed)
    ni_langmuir = cs3(alt_langmuir)

    # Create the rocket velocity profile for langmuir
    vec_rkt = np.array([data_dict_traj['N_VEL'][0],data_dict_traj['T_VEL'][0],data_dict_traj['P_VEL'][0]]).T
    v_rkt = np.array([np.linalg.norm(vec) for vec in vec_rkt])
    T0_trajectory = stl.EpochTo_T0_Rocket(data_dict_traj['Epoch'][0],T0=T0)
    cs4 = CubicSpline(T0_trajectory,v_rkt)
    vrkt_langmuir = cs4(T0_LP_current)

    # store everything in the output data dict

    data_dict_output = {**data_dict_output,
                        **{
                            'Epoch':deepcopy(data_dict_LP_current['Epoch']),
                            'v_rkt':[vrkt_langmuir,{'DEPEND_0':'Epoch','UNITS':'m/s'}],
                            'Ti': [Ti_langmuir, {'DEPEND_0': 'Epoch', 'UNITS': 'eV'}],
                            'ni': [ni_langmuir, {'DEPEND_0': 'Epoch', 'UNITS': 'm^-3'}],
                            'Alt':[alt_langmuir,{'DEPEND_0': 'Epoch', 'UNITS': 'km'}]
                           }
                        }



    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:
        stl.prgMsg('Creating output file')
        fileoutName = f'ACESII_{ACESII.payload_IDs[wRocket-4]}_plasma_potential_ESICAT_guess.cdf'
        outputPath = f'C:\Data\ACESII\calibration\LP_density_calibration\\{ACESII.fliers[wRocket-4]}\\' + fileoutName
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)

# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
LP_density_cal_plasma_potential(wRocket)

