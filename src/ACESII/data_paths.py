# --- data_paths.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Place to store the pathing information of the project

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
# -------------------

# --- imports ---
from glob import glob


class DataPaths:

    # --- --- --- --- --- --- ---
    # --- USER SPECIFIC DATA ---
    # --- --- --- --- --- --- ---
    user = 'connor'
    # PATH_TO_DATA_FOLDER = r'C:/Data/'
    PATH_TO_DATA_FOLDER = r'/home/connor/Data/ROCKETS/'
    HOMEDRIVE = 'C:'
    HOMEPATH = 'C:\\'
    fliers = ['high', 'low']

    ACES_data_folder = fr'{PATH_TO_DATA_FOLDER}ACESII/'
    TRICE_data_folder = fr'{PATH_TO_DATA_FOLDER}TRICEII\\'
    Integration_data_folder = fr'{PATH_TO_DATA_FOLDER}Integration\\'
