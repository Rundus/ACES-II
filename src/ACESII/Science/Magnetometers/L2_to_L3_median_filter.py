# --- L2_to_L3_median_filter.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Input some ACES-II RingCore magnetometer data, determine the coordinate system, then median filter it


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2026-03-30"
__version__ = "1.0.0"

import matplotlib.pyplot as plt

# --- --- --- --- ---

# --- Filter Toggles ---
fs = 256 # sample frequency of the data
low_cutoff = 1/60 # 20 seconds (or 1/20 freq) butterworth cutoff
order = 4
filter_window = 30*fs+1 # corresponds to 30 seconds

######################
# --- DATA TOGGLES ---
######################
just_print_file_names_bool = False
rocket_str = 'low'
wInstr = 'MAG'
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['L2', [[0],[0]]],
}
outputData = False


#################
# --- IMPORTS ---
#################
from scipy.integrate import simpson
from src.ACESII.data_tools.my_imports import *
from scipy.signal import savgol_filter
start_time = time.time()


def L2_to_L3_median_filter(data_dicts):

    #######################
    # --- LOAD THE DATA ---
    #######################
    stl.prgMsg(f'Loading data')
    data_dict = deepcopy(data_dicts[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = deepcopy(data_dict)

    #########################################
    # --- DETERMINE THE COORDINATE SYSTEM ---
    #########################################

    def wCoordUsed(dict_keys):
        coord_trio = [['N', 'T', 'p'], ['E', 'N', 'U'], ['r', 'e', 'p'], ['X', 'Y', 'Z']]

        for trio in coord_trio:
            B_keys = ['B_' + strV for strV in trio]
            if all(target in dict_keys for target in B_keys):
                return trio

        raise Exception('Could not determine Coordinate System')

    coord_keys = wCoordUsed(data_dict.keys())
    B_keys = ['B_' + strV for strV in coord_keys]

    ############################
    # --- MEDIAN FILTER DATA ---
    ############################
    stl.prgMsg('Median Filtering Data')
    from scipy.signal import medfilt


    # Play around with filtering
    N = 17500
    Bdata = np.array([data_dict_output[key][0][N:] for key in B_keys])

    # --- Plot the data ---
    fig,ax = plt.subplots()

    # raw data
    filter_window = 30*fs + 1
    dataChoice = Bdata[1]
    filtData = savgol_filter(x=dataChoice, window_length=filter_window, polyorder=3)
    ax.plot(data_dict_output['Epoch'][0][N:], dataChoice)
    ax.plot(data_dict_output['Epoch'][0][N:],filtData,color='red')
    ax.legend()
    plt.show()


    # remove the detrend

    # SSA
    wL = 500
    F_ssa = stl.SSA(dataChoice - filtData, wL)
    F_ssa.components_to_df().plot()
    plt.show()




    # # Filter the data
    # for key in B_keys:
    #
    #     filtData = deepcopy(data_dict[f'{key}'][0])
    #
    #     # --- butterworth filter the data to remove any spin tone effects ---
    #     filtData = stl.butterFilter().butter_filter(data=filtData,
    #                                                 lowcutoff=low_cutoff,
    #                                                 highcutoff=low_cutoff,
    #                                                 order=order,
    #                                                 fs=fs,
    #                                                 filtertype='lowPass'
    #                                                 )
    #
    #     # Median filter the data
    #     # filtData = medfilt(volume=filtData, kernel_size=filter_window)
    #
    #     # savitz-golay filter the data
    #     filtData = savgol_filter(x=filtData, window_length=filter_window, polyorder=3)
    #
    #
    #     # Append the Key
    #     data_dict_output[f'{key}_median'] = [filtData,deepcopy(data_dict[key][1])]


    stl.Done(start_time)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:
        stl.prgMsg('Creating output file')
        fileoutName = f'ACESII_{ACESII.fliers_dict[rocket_str]}_L3_RingCore_auroral_median_filter.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/L3/{wInstr}/{rocket_str}//{fileoutName}'
        stl.outputDataDict(outputPath, data_dict_output)
        stl.Done(start_time)



#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(L2_to_L3_median_filter,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)

