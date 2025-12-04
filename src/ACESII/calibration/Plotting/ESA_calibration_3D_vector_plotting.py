# --- ESA_calibration_3D_vector_plotting.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import itertools
# --- --- --- --- ---

import time

import matplotlib.pyplot as plt

from ACESII_code.class_var_func import Done, setupPYCDF

start_time = time.time()
# --- --- --- --- ---



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

modifier = ''
inputPath_modifier = 'mag' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
wInstr = 0 # 0 - eepaa 1- leesa 2- ieepaa


outputData = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import os
from tqdm import tqdm
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, prgMsg, L2_TRICE_Quick
from myspaceToolsLib.CDF_load import outputCDFdata
from myspaceToolsLib.transforms import Rz,Ry,Rx
from glob import glob
from os.path import getsize

setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)


def main(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]

    ############################
    # --- 3D vector plotting ---
    ###########################
    pitch = [-10 + i*10 for i in range(21)]

    # Define the SPHERICAL description of the unit vectors IN THE MAG FRAME
    if wInstr in [0,2]: # EEPAA
        ThetaPolar = np.radians(90)
        unit_vect_temp = [ [np.sin(ThetaPolar)*np.cos(np.radians(-90 - ptch)),np.sin(ThetaPolar)*np.sin(np.radians(-90 - ptch)),np.cos(ThetaPolar)   ] for ptch in pitch]
        unit_vect = [np.matmul(Ry(45),temp_vect) for temp_vect in unit_vect_temp]

    elif wInstr == 1: # LEESA
        ThetaPolar = np.radians(90)
        unit_vect = np.array([[np.sin(ThetaPolar) * np.cos(np.radians(ptch - 90)), np.sin(ThetaPolar) * np.sin(np.radians(ptch - 90)), np.cos(ThetaPolar)] for ptch in pitch])


    # --- PLOTTING ---
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # color the 0deg arrow red
    colors = ['black' for i in range(len(pitch))]
    colors[1] = 'magenta'

    # plot the vectors
    for i in range(len(unit_vect)):
        ax.quiver(0, 0, 0, unit_vect[i][0], unit_vect[i][1], unit_vect[i][2], color=colors[i])

    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim([-1.5, 1.5])
    ax.set_zlim([-1.5, 1.5])

    # plot the integration_tad_files axes
    ax.quiver(0, 0, 0, 1, 0, 0, color='red')
    ax.quiver(0, 0, 0, 0, 1, 0, color='green')
    ax.quiver(0, 0, 0, 0, 0, 1, color='blue')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()








# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

main(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
