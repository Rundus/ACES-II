# --- MPI_ENU_create_reduced_dataset.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: isolate points of the MPI data

# imports
import spaceToolsLib as stl
from src.ACESII.data_paths import DataPaths


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

    # 3. define a window of time where all points from the MPIs must fall within
    # # Idea: find the point in each MPI dataset which is nearest to the target_times
    # Then use your window to decide if that point (for each MPI) is within the time window you've chosen
    time_window = 0.05

    # 4. collect all the points into a new dataset and write them out in a data_dict_output.
    # note: write each MPI data individually, not as a single dataset.


    # --- 9. Output new CDF ---
    file_out_path = rf'{DataPaths.ACES_data_folder}/L3/MPI/low/ACESII_35364_L3_MPI_ENU_reduced.cdf'
    stl.outputDataDict(outputPath=file_out_path, data_dict=data_dict_output)


###EXECUTE###
MPI_ENU_create_reduced_dataset()
