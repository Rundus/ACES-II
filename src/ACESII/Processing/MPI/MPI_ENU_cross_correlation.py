# --- MPI_ENU_create_reduced_dataset.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: isolate points of the MPI data via cross-correlation
import matplotlib.pyplot as plt
# imports
import spaceToolsLib as stl
from src.ACESII.data_paths import DataPaths
import numpy as np
from scipy.signal import correlate


def MPI_ENU_cross_correlation():

    # 1. Load MPI ENU data
    path_to_MPI = rf'{DataPaths.ACES_data_folder}/L2/low/ACESII_35364_L2_MPI_ENU.cdf'
    data_dict_MPI = stl.loadDictFromFile(path_to_MPI)

    # prepare the output
    data_dict_output = {}

    # 2. define a regularized time grid to bin all the MPI data onto
    # Note: It MUST be finer resolution than any of the MPI points so that
    # multiple datapoints aren't in the same time bin
    N = 1000
    time_space = np.linspace(data_dict_MPI['time1'][0][0], data_dict_MPI['time1'][0][-1]+5, N)
    MPI_E_time_grided = [np.zeros_like(time_space) for i in range(4)]
    MPI_N_time_grided = [np.zeros_like(time_space) for i in range(4)]

    for idx in range(4):

        time = data_dict_MPI[f'time{idx+1}'][0]
        MPI_East = data_dict_MPI[f'MPI{idx + 1}_E'][0]
        MPI_North = data_dict_MPI[f'MPI{idx + 1}_N'][0]

        for tdx, tme in enumerate(time):
            target_idx = np.abs(time_space-tme).argmin()
            MPI_E_time_grided[idx][target_idx] = MPI_East[tdx]
            MPI_N_time_grided[idx][target_idx] = MPI_North[tdx]


    fig, ax = plt.subplots(3)
    for i in range(4):
        ax[0].plot(time_space, MPI_E_time_grided[i])
        ax[1].plot(time_space, MPI_N_time_grided[i])

    # just plot the multipliction fo everything
    ax[2].plot(time_space, np.abs(MPI_E_time_grided[0]*MPI_E_time_grided[3]*MPI_E_time_grided[1]))

    plt.show()



    #
    # # # 4. put our data in an output
    # #
    # # data_dict_output = {**data_dict_output,
    # #                     **{'time': [np.array(target_times), {}]},
    # #                     **{f'MPI{j}_E_reduced': [np.array(MPI_E_reduced[j]), {'DEPEND_0':'time'}] for j in range(4)},
    # #                     **{f'MPI{k}_N_reduced': [np.array(MPI_N_reduced[k]), {'DEPEND_0': 'time'}] for k in range(4)},
    # #                     }
    # #
    # #
    # # # --- 9. Output new CDF ---
    # # file_out_path = rf'{DataPaths.ACES_data_folder}/L3/MPI/low/ACESII_35364_L3_MPI_ENU_reduced.cdf'
    # # stl.outputDataDict(outputPath=file_out_path, data_dict=data_dict_output)


###EXECUTE###
MPI_ENU_cross_correlation()
