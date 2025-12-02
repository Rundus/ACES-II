# --- MPI_ENU_create_reduced_dataset.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: isolate points of the MPI data

# imports
import numpy as np
from copy import deepcopy
import spaceToolsLib as stl
from scipy.interpolate import CubicSpline
from src.data_paths import DataPaths


def MPI_ENU_create_reduced_dataset():

    # 1. Load MPI ENU data
    path_to_MPI = rf'{DataPaths.ACES_data_folder}/L2/low/ACESII_35364_L2_MPI_ENU.cdf'
    data_dict_MPI = stl.loadDictFromFile(path_to_MPI)

    # prepare the output
    data_dict_output = {}

    # 2. define some ideal time points as seconds from launch
    # Note: ONLY pick times where there's good correlation between MPI_East AND MPI_North
    target_times = [
        187.721,
        187.95,
        193.2,
        .... # and so on
    ]


    # --- 9. Output new CDF ---
    file_out_path = rf'{DataPaths.ACES_data_folder}/L3/MPI/low/ACESII_35364_L3_MPI_ENU_reduced.cdf'
    stl.outputDataDict(outputPath=file_out_path, data_dict=data_dict_output)


###EXECUTE###
MPI_ENU_create_reduced_dataset()
