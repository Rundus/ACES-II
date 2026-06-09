# --- ACESII_EISCAT_Compare_ion_Velocity.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Use the ExB drift data, the MPI data and the EISCAT ion drift speeds at specific locations to
# try to find a coherent "story" of ion drifts

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2026-03-30"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np

# --- --- --- --- ---


######################
# --- DATA TOGGLES ---
######################
just_print_file_names_bool = False
rocket_str = 'low'
wInstr = 'MPI'
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['L3', [[0],[0]]],
    '':['science/ExB',[[0],[0]]],
    'attitude':['',[[0],[0]]],
}
plot_ExB_MPI_compare = True
outputData = False


#################
# --- IMPORTS ---
#################
from scipy.integrate import simpson
from src.ACESII.data_tools.my_imports import *
start_time = time.time()
import datetime as dt
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import griddata


def ACESII_EISCAT_Compare_ion_velocity(data_dicts):

    #######################
    # --- LOAD THE DATA ---
    #######################
    stl.prgMsg(f'Loading data')
    data_dict_MPI = deepcopy(data_dicts[0])
    data_dict_ExB = deepcopy(data_dicts[1])
    data_dict_attitude = deepcopy(data_dicts[2])
    data_dict_vIon = stl.loadDictFromFile(glob('/home/connor/Data/ROCKETS/ACESII/science/EISCAT/tromso/UHF/ion_drifts/*.cdf*')[0])
    data_dict_EISCAT_highRez = stl.loadDictFromFile('/home/connor/Data/ROCKETS/ACESII/science/EISCAT/tromso/UHF/MAD6301_2022-11-20_beata_ant@uhfa.cdf')
    data_dict_EISCAT_lowRez = stl.loadDictFromFile('/home/connor/Data/ROCKETS/ACESII/science/EISCAT/tromso/UHF/MAD6400_2022-11-20_beata_ant@uhfa.cdf')
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
                       }

    # --- Clean up the ExB ---
    # Replace data before E-Field boom deployment with NaNs
    bad_point = dt.datetime(2022, 11, 20, 17, 23, 40)
    bad_idx = np.abs(data_dict_ExB['Epoch'][0] - bad_point).argmin()
    data_dict_ExB['ExB_E'][0][:bad_idx] = np.nan
    data_dict_ExB['ExB_N'][0][:bad_idx] = np.nan
    data_dict_ExB['ExB_Up'][0][:bad_idx] = np.nan

    # remove outliers in ExB
    bad_idxs = np.where(np.abs(data_dict_ExB['ExB_E'][0])>2000)
    data_dict_ExB['ExB_E'][0][bad_idxs] = np.nan
    bad_idxs = np.where(np.abs(data_dict_ExB['ExB_N'][0]) > 2000)
    data_dict_ExB['ExB_N'][0][bad_idxs] = np.nan

    # --- Isolate a few EISCAT vIon Profiles ---
    target_times = [
        # dt.datetime(2022, 11, 20, 16, 44, 10),
        # dt.datetime(2022, 11, 20, 16, 48, 10),
        # dt.datetime(2022, 11, 20, 16, 52, 10), # "Good"
        # dt.datetime(2022, 11, 20, 16, 56, 10),
        # dt.datetime(2022, 11, 20, 17, 0, 10),
        # dt.datetime(2022, 11, 20, 17, 4, 10),
        # dt.datetime(2022, 11, 20, 17, 8, 10),
        # dt.datetime(2022, 11, 20, 17, 12, 10),
        # dt.datetime(2022, 11, 20, 17, 16, 10),
        dt.datetime(2022, 11, 20, 17, 20, 10),
        # dt.datetime(2022, 11, 20, 17, 24, 10),
        # dt.datetime(2022, 11, 20, 17, 28, 10),
        # dt.datetime(2022, 11, 20, 17, 32, 10),
        # dt.datetime(2022, 11, 20, 17, 36, 10),
        # dt.datetime(2022, 11, 20, 17, 40, 10),
    ]

    Vel_ions_E = []
    Vel_ions_N = []
    Vel_ions_U = []
    Vel_ions_E_err = []
    Vel_ions_N_err = []
    Vel_ions_U_err = []
    alt_ions = data_dict_vIon['gdalt'][0]
    for tmeIdx,tmeVal in enumerate(target_times):
        closest_idx = np.abs(data_dict_vIon['Epoch'][0] - tmeVal).argmin()
        Vel_ions_E.append(deepcopy(data_dict_vIon['vi1'][0][closest_idx, :]))
        Vel_ions_N.append(deepcopy(data_dict_vIon['vi2'][0][closest_idx, :]))
        Vel_ions_U.append(deepcopy(data_dict_vIon['vi3'][0][closest_idx, :]))
        Vel_ions_E_err.append(deepcopy(data_dict_vIon['dvi1'][0][closest_idx, :]))
        Vel_ions_N_err.append(deepcopy(data_dict_vIon['dvi2'][0][closest_idx, :]))
        Vel_ions_U_err.append(deepcopy(data_dict_vIon['dvi3'][0][closest_idx, :]))

    # average all the time slices together
    Vel_ions_ENU = np.array([ np.nanmean(Vel_ions_E,axis=0),
                              np.nanmean(Vel_ions_N,axis=0),
                              np.nanmean(Vel_ions_U,axis=0)])

    Vel_ions_error_ENU = np.array([np.nanmean(Vel_ions_E_err, axis=0),
                             np.nanmean(Vel_ions_N_err, axis=0),
                             np.nanmean(Vel_ions_U_err, axis=0)])

    alt_val = [alt_ions,alt_ions]

    # --- Separate the data into upleg/downleg ---
    apogee_time_idx = np.abs(data_dict_attitude['Alt'][0] - np.max(data_dict_attitude['Alt'][0])).argmin()
    apogee_time = data_dict_attitude['Epoch'][0][apogee_time_idx]

    # MPI
    MPI_ENU = np.array([data_dict_MPI['MPI_E'][0],
                        data_dict_MPI['MPI_N'][0],
                        data_dict_MPI['MPI_U'][0]]).T
    MPI_idx = np.abs(data_dict_MPI['Epoch'][0] - apogee_time).argmin()


    # ExB
    ExB_ENU = np.array([data_dict_ExB['ExB_E'][0],
                        data_dict_ExB['ExB_N'][0],
                        data_dict_ExB['ExB_Up'][0]]).T
    ExB_idx = np.abs(data_dict_ExB['Epoch'][0] - apogee_time).argmin()


    # --- Fill the output data dict ---
    for i in range(3):
        Ekeys = ['E','N','Up']
        keys = ['E', 'N', 'U']
        data_dict_output = {
            **data_dict_output,
            **{
                f'ExB_{keys[i]}_upleg':[deepcopy(data_dict_ExB[f'ExB_{Ekeys[i]}'][0][0:ExB_idx]),deepcopy(data_dict_ExB[f'ExB_{Ekeys[i]}'][1])],
                f'ExB_{keys[i]}_downleg': [deepcopy(data_dict_ExB[f'ExB_{Ekeys[i]}'][0][ExB_idx:]), deepcopy(data_dict_ExB[f'ExB_{Ekeys[i]}'][1])],

                f'MPI_{keys[i]}_upleg': [deepcopy(data_dict_MPI[f'MPI_{keys[i]}'][0][0:MPI_idx]), deepcopy(data_dict_MPI[f'MPI_{keys[i]}'][1])],
                f'MPI_{keys[i]}_downleg': [deepcopy(data_dict_MPI[f'MPI_{keys[i]}'][0][MPI_idx:]), deepcopy(data_dict_MPI[f'MPI_{keys[i]}'][1])],

            }
        }

        if i == 0:
            for ephm in ['Alt','ILat']:
                data_dict_output = {
                    **data_dict_output,
                    **{
                        f'ExB_{ephm}_upleg':[data_dict_ExB[ephm][0][0:ExB_idx],deepcopy(data_dict_ExB[ephm][1])],
                        f'ExB_{ephm}_downleg': [data_dict_ExB[ephm][0][ExB_idx:], deepcopy(data_dict_ExB[ephm][1])],
                        f'MPI_{ephm}_upleg': [data_dict_MPI[ephm][0][0:MPI_idx], deepcopy(data_dict_MPI[ephm][1])],
                        f'MPI_{ephm}_downleg': [data_dict_MPI[ephm][0][MPI_idx:], deepcopy(data_dict_MPI[ephm][1])],
                    }
                }
    for keyVal in data_dict_output.keys():
        data_dict_output[keyVal][1]['DEPEND_0'] = None



    if plot_ExB_MPI_compare:
        #########################
        # --- PLOT EVERYTHING ---
        #########################
        fig, ax = plt.subplots(nrows=2,ncols=2)
        width = 12
        height = 1.618 * width
        fig.set_figwidth(width)
        fig.set_figheight(height)

        # UPLEG
        ax[0, 0].plot(data_dict_ExB['ExB_E'][0][0:ExB_idx],data_dict_ExB['Alt'][0][0:ExB_idx],color='tab:gray')
        ax[0, 0].plot(data_dict_MPI['MPI_E'][0][0:MPI_idx], data_dict_MPI['Alt'][0][0:MPI_idx], color='tab:red')
        ax[0, 0].errorbar(Vel_ions_ENU[0], alt_val[0], xerr=Vel_ions_error_ENU[0],color='tab:blue', capsize=10,elinewidth=1)
        ax[0, 0].scatter(Vel_ions_ENU[0], alt_val[0], s=10, color='tab:blue')
        ax[0, 0].set_title('E_East (Upleg)')

        ax[0, 1].plot(data_dict_ExB['ExB_N'][0][0:ExB_idx], data_dict_ExB['Alt'][0][0:ExB_idx], color='tab:gray')
        ax[0, 1].plot(data_dict_MPI['MPI_N'][0][0:MPI_idx], data_dict_MPI['Alt'][0][0:MPI_idx], color='tab:red')
        ax[0, 1].errorbar(Vel_ions_ENU[1], alt_val[1], xerr=Vel_ions_error_ENU[1], color='tab:blue', capsize=10,elinewidth=1)
        ax[0, 1].scatter(Vel_ions_ENU[1], alt_val[1], s=10, color='tab:blue')
        ax[0, 1].set_title('E_North (Upleg)')

        # DOWNLEG
        ax[1, 0].plot(data_dict_ExB['ExB_E'][0][ExB_idx::], data_dict_ExB['Alt'][0][ExB_idx::], color='tab:gray')
        ax[1, 0].plot(data_dict_MPI['MPI_E'][0][MPI_idx::], data_dict_MPI['Alt'][0][MPI_idx::], color='tab:red')
        ax[1, 0].errorbar(Vel_ions_ENU[0], alt_val[0], xerr=Vel_ions_error_ENU[0], color='tab:blue', capsize=10,elinewidth=1)
        ax[1, 0].scatter(Vel_ions_ENU[0], alt_val[0], s=10, color='tab:blue')
        ax[1, 0].set_title('E_East (Downleg)')

        ax[1, 1].plot(data_dict_ExB['ExB_N'][0][ExB_idx::], data_dict_ExB['Alt'][0][ExB_idx::], color='tab:gray')
        ax[1, 1].plot(data_dict_MPI['MPI_N'][0][MPI_idx::], data_dict_MPI['Alt'][0][MPI_idx::], color='tab:red')
        ax[1, 1].errorbar(Vel_ions_ENU[1], alt_val[1], xerr=Vel_ions_error_ENU[1], color='tab:blue', capsize=10, elinewidth=1)
        ax[1, 1].scatter(Vel_ions_ENU[1], alt_val[1], s=10, color='tab:blue')
        ax[1, 1].set_title('E_North (Downleg)')

        for i in range(2):
            for j in range(2):
                ax[i, j].set_xlim(-1000,1000)
                ax[i, j].set_ylim(90,190)
                ax[i,j].axhline(y=120,color='red',linestyle='--',alpha=0.5)
        plt.tight_layout()
        plt.show()

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')
        # write out the data
        fileoutName = f'ACESII_{ACESII.fliers_dict[rocket_str]}_neutral_wind.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/science/neutral_wind/{fileoutName}'
        stl.outputDataDict(outputPath, data_dict=data_dict_output)
        stl.Done(start_time)


#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(ACESII_EISCAT_Compare_ion_velocity,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)
