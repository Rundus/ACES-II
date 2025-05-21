# --- traj_to_auroral_coordinates.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: input the payload trajectory files and convert the ECEF velocities, positions into auroral coordinates
# as well as determine the gradient in the distances

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
import numpy as np
from src.my_imports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none

def traj_to_auroral_coordinates(wRocket):

    # --- Load the trajectory Data ---
    inputFiles = glob(f'{DataPaths.ACES_data_folder}\\trajectories\{ACESII.fliers[wRocket-4]}\\*Flight_trajectory_GPSdata*')[0]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles[i].replace(f'{DataPaths.ACES_data_folder}\{ACESII.fliers[wRocket-4]},''), round(getsize(file) / (10 ** 6), 1)))
        return


    # --- get the data from the tmCDF file ---
    stl.prgMsg(f'Loading data')
    data_dict_traj = stl.loadDictFromFile(inputFiles)
    data_dict_ECEF_auroral_transform = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\{ACESII.fliers[wRocket-4]}\*ECEF_to_auroral*')[0])
    stl.Done(start_time)

    # --- prepare output ---
    data_dict_output = {}
    
    #############################################
    # --- CONVERT DATA TO AURORAL COORDINATES ---
    #############################################

    # form the ECEF position vector
    rkt_pos_ECEF = np.array([data_dict_traj['ECEFXPOS'][0],data_dict_traj['ECEFYPOS'][0],data_dict_traj['ECEFZPOS'][0] ]).T


    # form the ECEF velocity vector
    rkt_vel_ECEF = np.array([data_dict_traj['ECEFXVEL'][0], data_dict_traj['ECEFYVEL'][0], data_dict_traj['ECEFZVEL'][0]]).T







    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')

        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

        outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod, wInstr[2])




# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
traj_to_auroral_coordinates(wRocket)