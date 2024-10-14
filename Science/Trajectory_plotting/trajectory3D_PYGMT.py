# --- trajectory3D_PYGMT.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: read in NOM_predict trajectories in order to roughly plot the ACESII rockets

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import itertools
import os.path
import time
from ACESII_code.class_var_func import Done, prgMsg,setupPYCDF
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files in
justPrintFileNames = False


# --- Select the Rocket ---
# 0 -> ACES II High Flier
# 1 -> ACES II Low Flier
wRocket = 0

# select which files to convert
# [] --> all files
# [#1,#2,...etc] --> only specific files
wFiles = []


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import csv
import numpy as np
from warnings import filterwarnings # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.
filterwarnings("ignore")
from tqdm import tqdm
from ACESII_code.missionAttributes import ACES_mission_dicts,TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder,fliers
from ACESII_code.class_var_func import color, L0_TRICE_Quick
from glob import glob
from os.path import getsize
import pygmt
import matplotlib.pyplot as plt
import pandas as pd

data_dicts = []
rocketFolderPath = ACES_data_folder

trajCSVFiles_low = glob(f'{rocketFolderPath}trajectories\{fliers[1]}\*.csv')
trajCSVFiles_high = glob(f'{rocketFolderPath}trajectories\{fliers[0]}\*.csv')
trajCSVFiles = [trajCSVFiles_high[0],trajCSVFiles_low[0]]


for wflyer in [0,1]:
    data_dict = {}
    vars = []

    with open(trajCSVFiles[wflyer]) as openedFile:
        csvFile = csv.reader(openedFile, delimiter=',')
        counter = 0
        for row in csvFile:
            if counter > 3:
                for i in range(len(row)):
                    data_dict[vars[i]][0].append(float(row[i]))

            elif counter == 2:
                vars = row
                for key in row:
                    data_dict = { **data_dict, **{key:[ [], {'VarName':key,'UNIT':'nan'}] } }

            elif counter == 3:
                for i in range(len(row)):
                    data_dict[vars[i]][1]['UNIT'] = row[i]

            counter +=1

    # Convert to numpy arrays
    for key, val in data_dict.items():
        data_dict[key][0] = np.array(data_dict[key][0])

    data_dicts.append(data_dict)



# --- --- --- --- ---
# --- PLOTTING ---
# --- --- --- --- ---


# --- Toggles ---
fig = pygmt.Figure()
region_rocket = [8, 20, 65, 76, 0, 400] #
perspective = [-120, 35]  # azimuth, elevation (in deg)
resolution = "15m"
# projection = "G16.020833/69.294167/12c+a0+t45+v60/80+w0+z400"
projection = "M15c"
styleR1 = "3p,red,-"
styleR2 = "3p,blue,-"
zscale = 0.05
frame =['+t"Andoya, Norway"','xa2g','ya2g']
cmap = "geo"
registration = "gridline"
registration = None



grid = pygmt.datasets.load_earth_relief(
    resolution = resolution,
    region = region_rocket,
    registration=registration)
fig.grdview(
    grid = grid,
    perspective = perspective,
    frame = frame,
    projection = projection,
    zsize = "1.5c",
    surftype = "s",
    cmap =cmap,
)


fig.plot3d(
    perspective = perspective,
    x=np.array(data_dicts[0]['Longitude'][0]),
    y=np.array(data_dicts[0]['Latitude'][0]),
    z=np.array(data_dicts[0]['Alt'][0]),
    region = region_rocket,
    zscale = zscale,
    pen = styleR1,
    frame = frame,
    label = 'High Flyer'
)

fig.plot3d(
    perspective = perspective,
    x=np.array(data_dicts[1]['Longitude'][0]),
    y=np.array(data_dicts[1]['Latitude'][0]),
    z=np.array(data_dicts[1]['Alt'][0]),
    region = region_rocket,
    zscale = zscale,
    pen = styleR2,
    frame = frame,
    label = 'Low Flyer'
)

plt.show()
# fig.savefig(r"D:\Data\ACESII\trajectories\trajectory_plots\testfig.png")














