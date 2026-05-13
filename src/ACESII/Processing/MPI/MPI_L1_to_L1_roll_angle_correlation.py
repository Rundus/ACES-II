# imports
import numpy as np
from copy import deepcopy
import spaceToolsLib as stl
import pandas as pd
from src.ACESII.data_tools.data_paths import DataPaths


def MPI_roll_association_check():

    # 1 Load data
    path_to_attitude = rf'{DataPaths.ACES_data_folder}/attitude/low/ACESII_36364_Attitude_Solution.cdf'
    path_to_MPI = fr'{DataPaths.ACES_data_folder}/L1/MPI/low/ACESII_36364_l1_MPI_rktFrm.cdf'

    data_dict_attitude = stl.loadDictFromFile(path_to_attitude)
    data_dict_MPI = stl.loadDictFromFile(path_to_MPI)

    #  prepare output
    data_dict_output = {}

    # 2 Extract MPI time + RollI
    time_epoch_attitude = data_dict_attitude['Epoch'][0]
    time_sec_attitude = stl.EpochTo_T0_Rocket(time_epoch_attitude,T0=time_epoch_attitude[0])
    rollI = np.array(data_dict_attitude['RollI'][0])

    # 3 loop over 4 MPIs

    for idx in range(4):
        # MPI time
        time_epoch_MPI = data_dict_MPI[f'Epoch_MPI{idx + 1}'][0]
        time_sec_MPI = stl.EpochTo_T0_Rocket(time_epoch_MPI,T0=time_epoch_attitude[0])

        N = len(time_sec_MPI)

        # 4 Interpolate roll onto MPI time base
        roll_MPI = np.interp(time_sec_MPI, time_sec_attitude, rollI)
        #
        # # 5 Extract MPI data
        # mpi_values = np.array(data_dict_MPI[f'Vx_rkt_MPI{idx + 1}'][0])
        #
        # # 6 panda data frame
        # df = pd.DataFrame({
        #     "roll": roll_MPI,
        #     "mpi": mpi_values
        # })
        #
        # # 7 Define roll bins
        # bin_width = 15  # degrees
        # df["roll_bin"] = (df["roll"] // bin_width) * bin_width
        #
        # # 8 Bin MPI values
        # binned = df.groupby("roll_bin")["mpi"].mean()
        #
        # roll_bins = binned.index.to_numpy()
        # mpi_binned = binned.to_numpy()

        # 9 Store into output dictionary (same structure you already use)
        data_dict_output = {
            **data_dict_output,
            **{
                f'Vx_rkt_MPI{idx + 1}': deepcopy(data_dict_MPI[f'Vx_rkt_MPI{idx + 1}']),
                f'Vy_rkt_MPI{idx + 1}': deepcopy(data_dict_MPI[f'Vy_rkt_MPI{idx + 1}']),
                f'Vz_rkt_MPI{idx + 1}': deepcopy(data_dict_MPI[f'Vz_rkt_MPI{idx + 1}']),
                #f'MPI{idx+1}_roll_binned': [np.array(mpi_binned), {'DEPEND_0': 'rollI'}],
                f'roll_MPI{idx+1}': [np.array(roll_MPI), {'DEPEND_0': None} ],
                f'Epoch_MPI{idx + 1}': deepcopy(data_dict_MPI[f'Epoch_MPI{idx + 1}']),
            }
        }

        # --- 10. Output new CDF ---
        outputPath = rf'{DataPaths.ACES_data_folder}/L2/MPI/low/'
        file_out_name = 'ACESII_36364_l2_MPI_roll_bin_analysis.cdf'
        stl.outputDataDict(outputPath=outputPath + file_out_name, data_dict=data_dict_output)

MPI_roll_association_check()