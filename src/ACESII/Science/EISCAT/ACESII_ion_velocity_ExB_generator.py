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
import spaceToolsLib

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
plot_ExB_MPI_compare = False
plot_EISCAT_model_compare = True
outputData = False
error_threshold = 1000 # [Percent] any radar points with error/value >= error_threshold are not included


#################
# --- IMPORTS ---
#################
from scipy.integrate import simpson
from src.ACESII.data_tools.my_imports import *
start_time = time.time()
import datetime as dt
from scipy.stats import binned_statistic
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
    data_dict_HWM14 = stl.loadDictFromFile('/home/connor/Data/ROCKETS/ACESII/science/HWM14_horizontal_neutral_wind_model/ACESII_HWM14.cdf')
    data_dict_bulkvel = spaceToolsLib.loadDictFromFile(glob('/home/connor/Data/MODELS/ACESII_ionosphere/bulk_velocity/*.cdf*')[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
                       }

    ##########################
    ##########################
    # --- Clean up the ExB ---
    ##########################
    ##########################
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

    #################################
    #################################
    # --- PREPARE THE EISCAT DATA ---
    #################################
    #################################
    # --- Isolate a few EISCAT vIon Profiles ---
    # Separate by upleg/downleg


    target_times_legs = [
        [
            # dt.datetime(2022, 11, 20, 16, 44, 10),
            # dt.datetime(2022, 11, 20, 16, 48, 10),
            # dt.datetime(2022, 11, 20, 16, 52, 10), # "Good"
            # dt.datetime(2022, 11, 20, 16, 56, 10),
            # dt.datetime(2022, 11, 20, 17, 0, 10),
            # dt.datetime(2022, 11, 20, 17, 4, 10),
            # dt.datetime(2022, 11, 20, 17, 8, 10),
            # dt.datetime(2022, 11, 20, 17, 12, 10),
            dt.datetime(2022, 11, 20, 17, 16, 10),
            dt.datetime(2022, 11, 20, 17, 20, 10),

            # dt.datetime(2022, 11, 20, 17, 24, 10),
            # dt.datetime(2022, 11, 20, 17, 28, 10), # "Good" downleg
            # dt.datetime(2022, 11, 20, 17, 32, 10),  # "good" downleg
            # dt.datetime(2022, 11, 20, 17, 36, 10),
        ], # "Good" on upleg East, poor on North
        [
            dt.datetime(2022, 11, 20, 17, 32, 10),  # "good" downleg
            dt.datetime(2022, 11, 20, 17, 36, 10),
            # dt.datetime(2022, 11, 20, 17, 40, 10),
        ], # "good" downleg (Bad on upleg)

    ]

    Vel_ions_means_upleg = [[],[],[]]
    Vel_ions_ENU_upleg = [[], [], []]
    Vel_ions_ENU_errors_upleg = [[], [], []]
    alts_upleg = [[], [], []]
    alts_mean_upleg = [[],[],[]]

    Vel_ions_means_downleg = [[], [], []]
    Vel_ions_ENU_downleg = [[], [], []]
    Vel_ions_ENU_errors_downleg = [[], [], []]
    alts_downleg = [[], [], []]
    alts_mean_downleg = [[], [], []]

    # For the upleg/downleg journeys, collect the good data, sort it, then error-weighted average it in altitude bins
    for k in range(len(target_times_legs)):
        target_times = target_times_legs[k]

        Vel_ions_raw = [[],[],[]]
        Vel_ions_err_raw = [[],[],[]]
        alt_ions = data_dict_vIon['gdalt'][0]
        for tmeIdx,tmeVal in enumerate(target_times):
            closest_idx = np.abs(data_dict_vIon['Epoch'][0] - tmeVal).argmin()
            for i in range(3):
                Vel_ions_raw[i].append(deepcopy(data_dict_vIon[f'vi{i+1}'][0][closest_idx, :]))
                Vel_ions_err_raw[i].append(deepcopy(data_dict_vIon[f'dvi{i+1}'][0][closest_idx, :]))

        # For each axis, collect only the valid points
        Vel_ions_ENU = [[],[],[]]
        Vel_ions_error_ENU = [[], [], []]
        alts = [[],[],[]]
        for i in range(3):

            data = Vel_ions_raw[i] # will look like [[],[],[],... for each EISCAT time slice]
            errors = Vel_ions_err_raw[i]

            # find only the good data
            for j in range(len(data)):
                good_idxs = np.where(np.isnan(data[j])==False)
                Vel_ions_ENU[i].extend(data[j][good_idxs])
                Vel_ions_error_ENU[i].extend(errors[j][good_idxs])
                alts[i].extend(alt_ions[good_idxs])

            # Sort the good data by altitude
            sorted_data = sorted(zip(alts[i], Vel_ions_ENU[i], Vel_ions_error_ENU[i]))
            alts_sorted, Vels_sorted, errors_sorted = map(list,zip(*sorted_data))
            Vel_ions_ENU[i] = Vels_sorted
            Vel_ions_error_ENU[i] = errors_sorted
            alts[i] = alts_sorted

        # Find the bin-weighted mean
        weighted_means = [[],[],[]]
        bin_altitudes = np.linspace(70, 200, int(1 + (200 - 70) / 2.5))
        for i in range(3):

            data = np.array(Vel_ions_ENU[i])
            data_error = np.array(Vel_ions_error_ENU[i])
            altitude = np.array(alts[i])

            # digitize the wind data based on altitude bins
            bin_idxs_wind = np.digitize(altitude,bin_altitudes)

            # cast the data into a dictionary
            grouped_wind = {
                j: data[bin_idxs_wind == j]
                for j in range(1, len(bin_altitudes))
            }

            # digitize the wind ERROR data based on altitude bins
            bin_idxs_wind_error = np.digitize(data_error, bin_altitudes)

            # cast the data into a dictionary
            grouped_wind_error = {
                j: data_error[bin_idxs_wind == j]
                for j in range(1, len(bin_altitudes))
            }

            # calculate the error-weighted means
            for key in grouped_wind.keys():
                if len(grouped_wind[key]) == 0:
                    weighted_means[i].append(np.nan)
                else:
                    values = grouped_wind[key]
                    weights = 1/grouped_wind_error[key]
                    weighted_means[i].append(np.average(values,weights=weights))

        # weighted_alts = (bin_altitudes[:-1]+bin_altitudes[1:])/2
        weighted_alts = bin_altitudes[:-1]


        # Remove any points with error above a certain threshold
        bad_alts_idxs = []
        for i in range(2): # I DONT care about the up direction
            # error_percent = np.abs((np.array(Vel_ions_error_ENU[i])/np.array(Vel_ions_ENU[i]))*100)
            # bad_idxs = np.where(error_percent >= error_threshold)[0]
            bad_idxs = np.where(np.abs(Vel_ions_error_ENU[i])>150)[0]
            if len(bad_idxs)!= 0:
                bad_alts_idxs.extend(bad_idxs)

        Vel_ions_ENU = np.delete(Vel_ions_ENU,bad_alts_idxs,axis=1)
        Vel_ions_error_ENU = np.delete(Vel_ions_error_ENU,bad_alts_idxs,axis=1)
        alts = np.delete(alts,bad_alts_idxs,axis=1)

        # store the output
        if k == 0:
            Vel_ions_means_upleg = weighted_means
            alts_mean_upleg = weighted_alts
            alts_upleg = alts
            Vel_ions_ENU_upleg = Vel_ions_ENU
            Vel_ions_ENU_errors_upleg = Vel_ions_error_ENU
        elif k ==1:
            Vel_ions_means_downleg = weighted_means
            alts_mean_downleg = weighted_alts
            alts_downleg = alts
            Vel_ions_ENU_downleg = Vel_ions_ENU
            Vel_ions_ENU_errors_downleg = Vel_ions_error_ENU

    ################################
    ################################
    # --- ISOLATE THE HWM14 DATA ---
    ################################
    ################################

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

    ############################
    ############################
    # --- CONSTRUCT V_I Data ---
    ############################
    ############################

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

                f'EISCAT_vi_{keys[i]}_upleg': [np.array(Vel_ions_ENU_upleg[i]),{}],
                f'EISCAT_vi_{keys[i]}_errors_upleg': [np.array(Vel_ions_ENU_errors_upleg[i]), {}],
                f'EISCAT_alts_{keys[i]}_upleg': [np.array(alts_upleg[i]), {}],

                f'EISCAT_vi_{keys[i]}_downleg': [np.array(Vel_ions_ENU_downleg[i]),{}],
                f'EISCAT_vi_{keys[i]}_errors_downleg': [np.array(Vel_ions_ENU_errors_downleg[i]), {}],
                f'EISCAT_alts_{keys[i]}_downleg': [np.array(alts_downleg[i]), {}],
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
        ax[0, 0].errorbar(Vel_ions_ENU_upleg[0], alts_upleg[0], xerr=Vel_ions_ENU_errors_upleg[0],color='tab:blue', capsize=10,elinewidth=1,linewidth=0)
        ax[0, 0].scatter(Vel_ions_ENU_upleg[0], alts_upleg[0], s=10, color='tab:blue')
        ax[0, 0].plot(data_dict_HWM14['zonal_wind'][0], data_dict_HWM14['altitude'][0])
        ax[0, 0].set_title('E_East (Upleg)')

        # ax[0,0].plot(weighted_means[0],weighted_alts,color='black')
        # ax[0, 0].scatter(weighted_means[0], weighted_alts, s=10, color='black')

        ax[0, 1].plot(data_dict_ExB['ExB_N'][0][0:ExB_idx], data_dict_ExB['Alt'][0][0:ExB_idx], color='tab:gray')
        ax[0, 1].plot(data_dict_MPI['MPI_N'][0][0:MPI_idx], data_dict_MPI['Alt'][0][0:MPI_idx], color='tab:red')
        ax[0, 1].errorbar(Vel_ions_ENU_upleg[1], alts_upleg[1], xerr=Vel_ions_ENU_errors_upleg[1], color='tab:blue', capsize=10,elinewidth=1,linewidth=0)
        ax[0, 1].scatter(Vel_ions_ENU_upleg[1], alts_upleg[1], s=10, color='tab:blue')
        ax[0, 1].plot(data_dict_HWM14['meridional_wind'][0], data_dict_HWM14['altitude'][0])
        ax[0, 1].set_title('E_North (Upleg)')

        # ax[0, 1].plot(weighted_means[1], weighted_alts,  color='black')
        # ax[0, 1].scatter(weighted_means[1], weighted_alts, s=10, color='black')


        # DOWNLEG
        ax[1, 0].plot(data_dict_ExB['ExB_E'][0][ExB_idx::], data_dict_ExB['Alt'][0][ExB_idx::], color='tab:gray')
        ax[1, 0].plot(data_dict_MPI['MPI_E'][0][MPI_idx::], data_dict_MPI['Alt'][0][MPI_idx::], color='tab:red')
        ax[1, 0].errorbar(Vel_ions_ENU_downleg[0], alts_downleg[0], xerr=Vel_ions_ENU_errors_downleg[0], color='tab:blue', capsize=10,elinewidth=1,linewidth=0)
        ax[1, 0].scatter(Vel_ions_ENU_downleg[0], alts_downleg[0], s=10, color='tab:blue')
        ax[1, 0].set_title('E_East (Downleg)')
        ax[1, 0].plot(data_dict_HWM14['zonal_wind'][0], data_dict_HWM14['altitude'][0])

        # ax[1, 0].plot(weighted_means[0], weighted_alts, color='black')
        # ax[1, 0].scatter(weighted_means[0], weighted_alts, s=10, color='black')

        ax[1, 1].plot(data_dict_ExB['ExB_N'][0][ExB_idx::], data_dict_ExB['Alt'][0][ExB_idx::], color='tab:gray')
        ax[1, 1].plot(data_dict_MPI['MPI_N'][0][MPI_idx::], data_dict_MPI['Alt'][0][MPI_idx::], color='tab:red')
        ax[1, 1].errorbar(Vel_ions_ENU_downleg[1], alts_downleg[1], xerr=Vel_ions_ENU_errors_downleg[1], color='tab:blue', capsize=10, elinewidth=1,linewidth=0)
        ax[1, 1].scatter(Vel_ions_ENU_downleg[1], alts_downleg[1], s=10, color='tab:blue')
        ax[1, 1].set_title('E_North (Downleg)')
        ax[1, 1].plot(data_dict_HWM14['meridional_wind'][0], data_dict_HWM14['altitude'][0])

        # ax[1, 1].plot(weighted_means[0], weighted_alts, color='black')
        # ax[1, 1].scatter(weighted_means[0], weighted_alts, s=10, color='black')

        for i in range(2):
            for j in range(2):
                ax[i, j].set_xlim(-1000,1000)
                ax[i, j].set_ylim(70,190)
                ax[i, j].axhline(y=120,color='red',linestyle='--',alpha=0.5)
        plt.tight_layout()
        plt.show()


    if plot_EISCAT_model_compare:
        #########################
        # --- PLOT EVERYTHING ---
        #########################
        fig, ax = plt.subplots(nrows=2,ncols=3)
        width = 14
        height = 1.618 * width
        fig.set_figwidth(width)
        fig.set_figheight(height)


        # Pick the location (L shell) of the Model data to choose. Average between that range
        v_i_East = [[],[]]
        v_i_North = [[], []]
        v_i_Up = [[], []]
        target_LShell = [[8.6, 9.2],[8.7,8.9]]
        for i in range(2):
            low_idx, high_idx = np.abs(data_dict_bulkvel['simLShell'][0] - target_LShell[i][0]).argmin(),np.abs(data_dict_bulkvel['simLShell'][0] - target_LShell[i][1]).argmin()
            v_i_East[i] = np.mean(data_dict_bulkvel['v_i_East'][0][low_idx:high_idx+1],axis=0)
            v_i_North[i] = np.mean(data_dict_bulkvel['v_i_North'][0][low_idx:high_idx+1],axis=0)
            v_i_Up[i] = np.mean(data_dict_bulkvel['v_i_Up'][0][low_idx:high_idx + 1], axis=0)
        simAlt = data_dict_bulkvel['simAlt'][0]/stl.m_to_km
        v_i_model = [v_i_East,v_i_North,v_i_Up]

        # UPLEG
        for i in range(3):
            ax[0, i].errorbar(Vel_ions_ENU_upleg[i], alts_upleg[0], xerr=Vel_ions_ENU_errors_upleg[i],color='tab:blue', capsize=10,elinewidth=1,linewidth=0)
            ax[0, i].scatter(Vel_ions_ENU_upleg[i], alts_upleg[0], s=10, color='tab:blue')
            ax[0, i].plot(v_i_model[i][0], simAlt)

            if i == 0:
                ax[0, i].set_title('E_East (Upleg)')
                ax[0, i].set_ylabel('Altitude [km]')
            elif i ==1:
                ax[0, i].set_title('E_North (Upleg)')
            elif i ==2:
                ax[0, i].set_title('E_North (Upleg)')

            # DOWNLEG
            ax[1, i].errorbar(Vel_ions_ENU_downleg[i], alts_downleg[i], xerr=Vel_ions_ENU_errors_downleg[i], color='tab:blue', capsize=10,elinewidth=1,linewidth=0)
            ax[1, i].scatter(Vel_ions_ENU_downleg[i], alts_downleg[i], s=10, color='tab:blue')
            ax[1, i].plot(v_i_model[i][1], simAlt)
            if i == 0:
                ax[0, i].set_title('E_East (Upleg)')
                ax[0, i].set_ylabel('Altitude [km]')
            elif i ==1:
                ax[0, i].set_title('E_North (Upleg)')
            elif i ==2:
                ax[0, i].set_title('E_North (Upleg)')

            for j in range(2):
                ax[i, j].set_xlim(-1000,1000)
                ax[i, j].set_ylim(70,190)
                ax[i, j].axhline(y=120,color='red',linestyle='--',alpha=0.5)
        plt.tight_layout()
        plt.show()


    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')
        # write out the data
        fileoutName = f'ACESII_{ACESII.fliers_dict[rocket_str]}_bulk_ion_flow.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/science/bulk_ion_flow/{fileoutName}'
        stl.outputDataDict(outputPath, data_dict=data_dict_output)
        stl.Done(start_time)


#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(ACESII_EISCAT_Compare_ion_velocity,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)
