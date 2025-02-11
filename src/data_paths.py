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
    user = 'cfelt'
    PATH_TO_DATA_FOLDER = r'C:\Data\\'
    # SYSTEM_PATH = f'C:\\Users\\{user}\\PycharmProjects\\UIOWA_CDF_operator\\ACESII_code'
    HOMEDRIVE = 'C:'
    HOMEPATH = 'C:\\'
    fliers = ['high', 'low']


    # --- --- --- --- --- ---
    # --- ACES II PATHS ---
    # --- --- --- --- --- ---
    ACES_data_folder = fr'{PATH_TO_DATA_FOLDER}ACESII\\'

    ACES_tad_files = [glob(f'{ACES_data_folder}tad\{fliers[0]}\*.tad'),glob(f'{ACES_data_folder}tad\{fliers[1]}\*.tad') ]
    ACES_tmCDF_files = [glob(f'{ACES_data_folder}tmCDF\{fliers[0]}\*.cdf'),glob(f'{ACES_data_folder}tmCDF\{fliers[1]}\*.cdf') ]
    ACES_L0_files = [glob(f'{ACES_data_folder}L0\{fliers[0]}\*.cdf'),glob(f'{ACES_data_folder}L0\{fliers[1]}\*.cdf')]
    ACES_L1_files = [glob(f'{ACES_data_folder}L1\{fliers[0]}\*.cdf'),glob(f'{ACES_data_folder}L1\{fliers[1]}\*.cdf')]
    ACES_L2_files = [glob(f'{ACES_data_folder}L2\{fliers[0]}\*.cdf'),glob(f'{ACES_data_folder}L2\{fliers[1]}\*.cdf')]
    ACES_csv_trajectories = [glob(f'{ACES_data_folder}trajectories\{fliers[0]}\*_GPSdata.csv'),glob(f'{ACES_data_folder}trajectories\{fliers[1]}\*_GPSdata.csv')]
    ACES_cdf_trajectories = [glob(f'{ACES_data_folder}trajectories\{fliers[0]}\*_GPSdata.cdf'),glob(f'{ACES_data_folder}trajectories\{fliers[1]}\*_GPSdata.cdf')]

    # --- --- --- --- --- ---
    # --- TRICE II PATHS ---
    # --- --- --- --- --- ---
    TRICE_data_folder = fr'{PATH_TO_DATA_FOLDER}TRICEII\\'

    TRICE_ACS_files = [glob(f'{TRICE_data_folder}attitude\{fliers[0]}\*.cdf'), glob(f'{TRICE_data_folder}attitude\{fliers[1]}\*.cdf')]
    TRICE_tad_files = [glob(f'{TRICE_data_folder}tad\{fliers[0]}\*.tad'),glob(f'{TRICE_data_folder}tad\{fliers[1]}\*.tad')]
    TRICE_tmCDF_files = [glob(f'{TRICE_data_folder}tmCDF\{fliers[0]}\*.cdf'),glob(f'{TRICE_data_folder}tmCDF\{fliers[1]}\*.cdf')]
    TRICE_L0_files = [glob(f'{TRICE_data_folder}L0\{fliers[0]}\*.cdf'),glob(f'{TRICE_data_folder}L0\{fliers[1]}\*.cdf')]
    TRICE_L1_files = [glob(f'{TRICE_data_folder}L1\{fliers[0]}\*.cdf'),glob(f'{TRICE_data_folder}L1\{fliers[1]}\*.cdf')]
    TRICE_L2_files = [glob(f'{TRICE_data_folder}L2\{fliers[0]}\*.cdf'),glob(f'{TRICE_data_folder}L2\{fliers[1]}\*.cdf')]



    # --- --- --- --- --- --- ---
    # --- INTEGRATION FILES ---
    # --- --- --- --- --- --- ---
    Integration_data_folder = fr'{PATH_TO_DATA_FOLDER}Integration\\'

    Integration_tad_files = [glob(fr'{PATH_TO_DATA_FOLDER}tad\{fliers[0]}\*.tad'),glob(fr'{PATH_TO_DATA_FOLDER}tad\{fliers[1]}\*.tad')]
    Integration_tmCDF_files =[glob(fr'{PATH_TO_DATA_FOLDER}tmCDF\{fliers[0]}\*.cdf'),glob(fr'{PATH_TO_DATA_FOLDER}tmCDF\{fliers[1]}\*.cdf')]
    Integration_L0_files = [glob(fr'{PATH_TO_DATA_FOLDER}L0\{fliers[0]}\*.cdf'),glob(fr'{PATH_TO_DATA_FOLDER}L0\{fliers[1]}\*.cdf')]
    Integration_L1_files = [glob(fr'{PATH_TO_DATA_FOLDER}L1\{fliers[0]}\*.cdf'),glob(fr'{PATH_TO_DATA_FOLDER}L1\{fliers[1]}\*.cdf')]
    Integration_L2_files = [glob(fr'{PATH_TO_DATA_FOLDER}L2\{fliers[0]}\*.cdf'),glob(fr'{PATH_TO_DATA_FOLDER}L2\{fliers[1]}\*.cdf')]