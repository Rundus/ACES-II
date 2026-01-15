# --- MPI_ENU_create_reduced_dataset.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: isolate points of the MPI data

# imports
import spaceToolsLib as stl
from src.ACESII.data_paths import DataPaths
import numpy as np


def MPI_ENU_create_reduced_dataset():

    # 1. Load MPI ENU data
    path_to_MPI = rf'{DataPaths.ACES_data_folder}/L2/low/ACESII_35364_L2_MPI_ENU.cdf'
    data_dict_MPI = stl.loadDictFromFile(path_to_MPI)

    # prepare the output
    data_dict_output = {}

    # 2. define some ideal time points as seconds from launch
    # Note: ONLY pick times where there's good correlation between MPI_East AND MPI_North
    target_times = [
        184.18,
        207.94,
        264.7,
    ]


    # 3. collect all the points into a new dataset and write them out in a data_dict_output.
    # note: write each MPI data individually, not as a single dataset.

    MPI_E_reduced = [[], [], [], []]
    MPI_N_reduced = [[], [], [], []]

    for idx in range(4):
        time = data_dict_MPI[f'time{idx + 1}'][0]
        MPI_East = data_dict_MPI[f'MPI{idx+1}_E'][0]
        MPI_North = data_dict_MPI[f'MPI{idx+1}_N'][0]

        for tme in target_times:

            target_idx = np.abs(time - tme).argmin()

            MPI_E_reduced[idx].append(MPI_East[target_idx])
            MPI_N_reduced[idx].append(MPI_North[target_idx])


    # 4. put our data in an output

    data_dict_output = {**data_dict_output,
                        **{'time': [np.array(target_times), {}]},
                        **{f'MPI{j}_E_reduced': [np.array(MPI_E_reduced[j]), {'DEPEND_0':'time'}] for j in range(4)},
                        **{f'MPI{k}_N_reduced': [np.array(MPI_N_reduced[k]), {'DEPEND_0': 'time'}] for k in range(4)},
                        }


    # --- 9. Output new CDF ---
    file_out_path = rf'{DataPaths.ACES_data_folder}/L3/MPI/low/ACESII_35364_L3_MPI_ENU_reduced.cdf'
    stl.outputDataDict(outputPath=file_out_path, data_dict=data_dict_output)


###EXECUTE###
MPI_ENU_create_reduced_dataset()
