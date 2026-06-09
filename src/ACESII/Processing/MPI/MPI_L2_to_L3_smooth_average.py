# --- MPI_L2_to_L3_smooth_average.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: collect all the MPI into one timebase and
# average it, then smooth it

import matplotlib.pyplot as plt
# imports
import numpy as np
from copy import deepcopy
import spaceToolsLib as stl
from scipy.interpolate import CubicSpline
from src.ACESII.data_tools.data_paths import DataPaths
import datetime as dt
from matplotlib.widgets import Slider
from scipy.signal import savgol_filter


plot_interactive_despin_MPI = True


def MPI_L2_to_L3_smooth_average():

    # Load the data
    path_to_MPI = rf'{DataPaths.ACES_data_folder}/L2/MPI/low/ACESII_36364_L2_MPI_ENU.cdf'
    data_dict_attitude = stl.loadDictFromFile(rf'{DataPaths.ACES_data_folder}/attitude/low/ACESII_36364_Attitude_Solution.cdf')
    data_dict_MPI = stl.loadDictFromFile(path_to_MPI)

    # --- PROCESS THE DATA ---


    # select 1 MPI as the "prime"
    wMPI = 4
    spin_rate = 0.5474
    deltaT = 3*(1/spin_rate)
    MPI_time_prime = data_dict_MPI[f'time{wMPI}'][0]
    MPI_E = [[val] for val in data_dict_MPI[f'MPI{wMPI}_E'][0]]
    MPI_N = [[val] for val in data_dict_MPI[f'MPI{wMPI}_N'][0]]
    MPI_U = [[val] for val in data_dict_MPI[f'MPI{wMPI}_U'][0]]

    for tmeIdx in range(len(MPI_time_prime)):

        # get the timepoint of the look MPI
        timeP = MPI_time_prime[tmeIdx]

        # search around this point to find the other MPI data that could be averaged together
        wMPI_list = [wMPI]
        for idx,mpiIdx in enumerate(wMPI_list):
            other_MPI_time = data_dict_MPI[f'time{mpiIdx}'][0]
            closest_idx = np.abs(other_MPI_time-timeP).argmin()
            closest_time = other_MPI_time[closest_idx]
            if np.abs(timeP-closest_time) <= deltaT:
                MPI_E[tmeIdx].append(data_dict_MPI[f'MPI{mpiIdx}_E'][0][closest_idx])
                MPI_N[tmeIdx].append(data_dict_MPI[f'MPI{mpiIdx}_N'][0][closest_idx])
                MPI_U[tmeIdx].append(data_dict_MPI[f'MPI{mpiIdx}_U'][0][closest_idx])


    # Average all the MPI data together
    MPI_E = np.array([np.nanmean(MPI_E[i]) for i in range(len(MPI_E))])
    MPI_N = np.array([np.nanmean(MPI_N[i]) for i in range(len(MPI_N))])
    MPI_U = np.array([np.nanmean(MPI_U[i]) for i in range(len(MPI_U))])
    MPI_ENU = np.array([MPI_E,MPI_N,MPI_U]).T

    # Determine the true epoch for this data
    Epoch = np.array([dt.datetime(2022,11,20,17,20) + dt.timedelta(seconds=data_dict_MPI[f'time{wMPI}'][0][i]) for i in range(len(data_dict_MPI[f'time{wMPI}'][0]))])

    # Smooth the data
    for i in range(3):
        MPI_ENU[:,i] = savgol_filter(x=MPI_ENU[:,i],
                                   window_length=8,
                                   polyorder=3)

    # Prepare the output
    data_dict_output = {
        'MPI_E':[MPI_ENU[:,0], deepcopy(data_dict_MPI[f'MPI{wMPI}_E'][1])],
        'MPI_N': [MPI_ENU[:,1], deepcopy(data_dict_MPI[f'MPI{wMPI}_N'][1])],
        'MPI_U': [MPI_ENU[:,2], deepcopy(data_dict_MPI[f'MPI{wMPI}_U'][1])],
        'Epoch': [Epoch,{'LABLAXIS':'Epoch'}]
    }
    for coord in ['E','N','U']:
        data_dict_output[f'MPI_{coord}'][1]['DEPEND_0'] = 'Epoch'
        data_dict_output[f'MPI_{coord}'][1]['UNITS'] = 'm/s'

    # --- OUTPUT THE DATA ----
    file_out_name = 'ACESII_36364_l3_MPI_ENU.cdf'
    outputPath = rf'{DataPaths.ACES_data_folder}/L3/MPI/low/' + file_out_name
    stl.outputDataDict(outputPath=outputPath,
                       data_dict=data_dict_output)






###EXECUTE###
MPI_L2_to_L3_smooth_average()
