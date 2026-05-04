# --- cal2_visualize_current_system.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: adjust the timebase of the ACS DCM to see if a cleaner MAG despin can be derived

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
wInstr = 'MAG'
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['calibration', [[0],[0]]],
    'EEPAA':['L2',[[0],[0]]],
    'trajectories':['',[[0],[0]]],
    '/':['science/ESA_currents/EEPAA',[[0],[0]]],
    '':['science/ESA_currents/IEPAA',[[0],[0]]]
}
outputData = False

# --- Filter Toggles ---
fs = 256 # sample frequency of the data
low_cutoff = 0.05 # 20 seconds (or 1/20 freq) butterworth cutoff
order = 4
filter_window = 30*fs+1 # corresponds to 30 seconds
plot_filtered_dB = True


#################
# --- IMPORTS ---
#################
from src.ACESII.data_tools.my_imports import *
from scipy.signal import medfilt
start_time = time.time()
from scipy.signal import savgol_filter
from scipy.interpolate import CubicSpline
import datetime as dt

def cal2_visualize_current_system(data_dicts):

    #######################
    # --- LOAD THE DATA ---
    #######################
    stl.prgMsg(f'Loading data')
    data_dict_MAG = deepcopy(data_dicts[0])
    data_dict_EEPAA = deepcopy(data_dicts[1])
    data_dict_traj = deepcopy(data_dicts[2])
    data_dict_Jpara_EEPAA = deepcopy(data_dicts[3])
    data_dict_Jpara_IEPAA = deepcopy(data_dicts[4])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = deepcopy(data_dict_MAG)


    # --- [0] Do some initial Filtering ---

    # === Trajectories ===


    # === B - Field ===

    def wCoordUsed(dict_keys):
        coord_trio = [['N', 'T', 'p'], ['E', 'N', 'U'], ['r', 'e', 'p'], ['X', 'Y', 'Z']]

        for trio in coord_trio:
            B_keys = ['B_' + strV for strV in trio]
            if all(target in dict_keys for target in B_keys):
                return trio

        raise Exception('Could not determine Coordinate System')

    coord_keys = wCoordUsed(data_dict_MAG.keys())
    B_keys = ['B_' + strV for strV in coord_keys]


    # Play around with filtering
    N = 10000
    Bdata = np.array([data_dict_output[key][0][N:] for key in B_keys]).T

    if plot_filtered_dB:
        # Filter the data
        fig, ax = plt.subplots(5,2)
        height = 10
        width= height*1.618
        fig.set_size_inches(width,height)
        xLimits = [dt.datetime(2022,11,20,17,23),
                   dt.datetime(2022,11,20,17,28)]

        for i, key in enumerate(['test','test']+B_keys):
            if i == 0:
                for k in range(2):
                    colormap = stl.apl_rainbow_black0_cmap()
                    colormap.set_bad((0,0,0))
                    cmap=ax[i,k].pcolormesh(data_dict_EEPAA['Epoch'][0],
                                       data_dict_EEPAA['Energy'][0],
                                       data_dict_EEPAA['Differential_Energy_Flux'][0][:,2,:].T,
                                       cmap=colormap,
                                       norm='log')
                    ax[i,k].set_yscale('log')
                    ax[i,k].set_ylabel('Energy [eV]')
                    ax[i,k].set_xlim(*xLimits)
                    ax[i,k].legend(loc='upper right')

            elif i == 1:
                for k in range(2):
                    eData = data_dict_Jpara_EEPAA['j_para'][0]/(1E-6)
                    iData =data_dict_Jpara_IEPAA['j_para'][0]/(1E-6)

                    iDta_filtered = stl.butterFilter().butter_filter(data=iData,
                                                                   lowcutoff=0.15,
                                                                   highcutoff=0.15,
                                                                   order=order,
                                                                   fs=20,
                                                                   filtertype='lowPass'
                                                                   )

                    ax[i,k].plot(data_dict_Jpara_EEPAA['Epoch'][0],eData,color='tab:blue',label=r'$J_{\parallel}^{e-}$')
                    ax[i, k].plot(data_dict_Jpara_IEPAA['Epoch'][0], -1*iDta_filtered, color='tab:red',label=r'$J_{\parallel}^{i+}$')
                    ax[i,k].set_ylim(-10,10)
                    ax[i,k].set_ylabel(r'Current [$\mu$A]')
                    ax[i, k].set_xlim(*xLimits)
                    ax[i,k].legend(loc='upper right')
            else:
                adjust = 2
                filtData = Bdata[:,i-adjust]

                # --- butterworth filter the data to remove any spin tone effects ---
                filtData = stl.butterFilter().butter_filter(data=filtData,
                                                            lowcutoff=low_cutoff,
                                                            highcutoff=low_cutoff,
                                                            order=order,
                                                            fs=fs,
                                                            filtertype='lowPass'
                                                            )

                # Median filter the data
                # filtData = medfilt(volume=filtData, kernel_size=filter_window)

                # savitz-golay filter the data
                # filtData = savgol_filter(x=filtData, window_length=filter_window, polyorder=3)

                # plot the filtered data
                data_dict_output[f'{key}_filter'] = [filtData,deepcopy(data_dict_output[key][1])]
                ax[i,0].plot(data_dict_output['Epoch'][0][N:], data_dict_output[f'{key}'][0][N:], label=f'{B_keys[i-adjust]}')
                ax[i,0].plot(data_dict_output['Epoch'][0][N:], filtData, label=f'{B_keys[i-adjust]}_filter')
                ax[i,0].set_ylabel(f'{key} [nT]')
                ax[i,0].legend(loc='upper right')
                ax[i, 0].set_xlim(*xLimits)

                dat = data_dict_output[f'{key}'][0][N:]
                dB_raw = np.diff(dat, prepend=dat[0])
                dB_adjust = np.diff(filtData, prepend=filtData[0])
                dB_filtData = stl.butterFilter().butter_filter(data=dB_raw,
                                                            lowcutoff=low_cutoff,
                                                            highcutoff=low_cutoff,
                                                            order=order,
                                                            fs=fs,
                                                            filtertype='lowPass'
                                                            )

                ax[i, 1].plot(data_dict_output['Epoch'][0][N:], dB_raw, label=f'{B_keys[i-adjust]}')
                ax[i, 1].plot(data_dict_output['Epoch'][0][N:], dB_adjust, label=f'{B_keys[i-adjust]}_filter_first')
                ax[i, 1].plot(data_dict_output['Epoch'][0][N:], dB_filtData, label=f'{B_keys[i-adjust]}_filter_after')
                ax[i,1].axhline(y=0,alpha=0.5,linestyle='--',color='tab:red')
                ax[i,1].legend(loc='upper right')
                ax[i, 1].set_ylabel(f'd{key} [nT]')
                ax[i,1].set_ylim(-0.2,0.2)
                ax[i, 1].set_xlim(*xLimits)
        plt.tight_layout()
        plt.show()

    # === [1] Remove small-scale noise through savgol filter ===
    B_sinwave = [[], [], []]
    fig, ax = plt.subplots(3,2)
    for i in range(3):
        filtData = Bdata[:, i ]

        # --- butterworth filter the data to remove any spin tone effects ---
        filtData = stl.butterFilter().butter_filter(data=filtData,
                                                    lowcutoff=0.15,
                                                    highcutoff=0.15,
                                                    order=4,
                                                    fs=fs,
                                                    filtertype='highpass'
                                                    )
        B_sinwave[i] = filtData
        ax[i,0].plot(data_dict_MAG['Epoch'][0][N:], Bdata[:,i],color='tab:blue')
        ax[i,0].plot(data_dict_MAG['Epoch'][0][N:], B_sinwave[i],color='tab:red')
        ax[i,1].plot(data_dict_MAG['Epoch'][0][N:],Bdata[:,i]-B_sinwave[i],color='tab:purple')

    plt.show()


    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')

        # write out the data
        fileoutName = f'ACESII_{ACESII.fliers_dict[rocket_str]}_l1_RingCore_cal1_DCM_time_adjust.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/calibration/{wInstr}/{rocket_str}//{fileoutName}'
        stl.outputDataDict(outputPath, data_dict=data_dict_output)
        stl.Done(start_time)




#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(cal2_visualize_current_system,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)

