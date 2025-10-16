# --- Plot_ASI_auroral_coordinates_rotation_vec.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:
# [1] Develop an algorthim to center the E-East around 0 mV/m
# [2] develop a method to rotate E-Field, field-aligned coordinate data in order to minimize
# the perpendicular E-Field. This is the "slab" approximation. Then, plot over All Sky imager to compare

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

from src.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False # Just print the names of files

# --- OutputData ---
outputData = False

# --- Plots ---

######################
# --- PLOT TOGGLES ---
######################
figure_height = (20)
figure_width = (20)
# ---------------BigAllSky-----------------
wavelength = ['5570', '6300']
trajColors = ['tab:red', 'tab:orange']
lonW = 8
lonE = 30
latS = 68.5
latN = 74
res = '50m'
cbarVmin,cbarVmax = 0, 12 # in kRayleigh
BigAllSky_textSize = 20
BigAllSky_tickLabelSize = 30
BigAllSky_lineThickness = 8
BigAllSky_GridSize = 5
BigAllSky_TitleSize = 40
BigAllSky_costLineSize = 3

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.my_imports import *
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import spaceToolsLib as stl
faceColorChoice = (156 / 255, 156 / 255, 156 / 255, 0.5)  # in normalize RGBA
cmapColor = stl.matlab_parula_cmap()


def Plot_ASI_auroral_coordinates_rotation_vec():

    # find the point (time, lat, long, alt) for when the Flyer's first hit the aurora vs when they left
    wpoints = 0  # 0 for HF 1 for LF
    image_set = [8, 9]
    wImage = image_set[wpoints]


    #######################
    # --- LOAD THE DATA ---
    #######################
    stl.prgMsg(f'Loading data')
    data_dict_allSky5577 = stl.loadDictFromFile(glob(r'C:\Data\ACESII\all_sky\skibotn\5577\\*.cdf')[0])
    data_dict_allSky6300 = stl.loadDictFromFile(glob(r'C:\Data\ACESII\all_sky\skibotn\6300\\*.cdf')[0])
    data_dicts_allsky = [data_dict_allSky5577,data_dict_allSky6300]
    data_dict_EFI_auroral = stl.loadDictFromFile(glob(r'C:\Data\ACESII\science\auroral_coordinates\low\*.cdf*')[0])
    attitudeFolderPath = f'{DataPaths.ACES_data_folder}\\attitude'
    inputFilesTraj = [glob(f'{attitudeFolderPath}\\{DataPaths.fliers[0]}\\*.cdf*')[0],
                      glob(f'{attitudeFolderPath}\\{DataPaths.fliers[1]}\\*.cdf*')[0]]
    data_dicts_attitude = [stl.loadDictFromFile(inputFilesTraj[0]), stl.loadDictFromFile(inputFilesTraj[1])]
    stl.Done(start_time)

    for i in range(2):
        # --- PLOT ---
        fig, ax = plt.subplots(constrained_layout=True)
        fig.set_figwidth(figure_width)
        fig.set_figheight(figure_height)

        # all sky image
        ax.set_title(f'Skibotn {wavelength[i]}$\AA$\n' + data_dicts_allsky[i]['Epoch'][0][wImage].strftime("%Y-%B-%d %H:%M:%S") + ' UTC',weight='bold')
        ax.imshow(data_dicts_allsky[i]['AllSkyImages'][0][wImage], cmap=cmapColor, vmin=cbarVmin, vmax=cbarVmax)

        # all center center coordinates
        Lat_origin = len(data_dict_allSky5577['AllSkyImages'][0][wImage]) / 2
        Long_origin = len(data_dict_allSky5577['AllSkyImages'][0][wImage]) / 2
        offset = 20
        ax.set_ylim(-1*offset, Lat_origin*2 + offset)
        ax.set_xlim(-1*offset, Lat_origin * 2 + offset)
        plt.gca().invert_yaxis()

        # Units vectors
        mag = 300
        East_uvec = mag * np.array([1, 0])
        North_uvec = mag * np.array([0, -1])

        def rotate(angle, vec):
            rad = np.radians(angle)  # the 90deg is to account for the imshow rotation
            matrix = np.array([[np.cos(rad), -1 * np.sin(rad)], [np.sin(rad), np.cos(rad)]])
            return np.matmul(matrix, vec)

        rotation_angle = -1 * data_dict_EFI_auroral['rotation_Angle'][0][0]
        East_uvec_rot = rotate(rotation_angle, East_uvec)
        North_uvec_rot = rotate(rotation_angle, North_uvec)
        east_vector = [[Long_origin, East_uvec_rot[0] + Long_origin], [Lat_origin, East_uvec_rot[1] + Lat_origin]]
        north_Vector = [[Long_origin, North_uvec_rot[0] + Long_origin], [Lat_origin, North_uvec_rot[1] + Lat_origin]]
        ax.plot(north_Vector[0], north_Vector[1], color='purple', linewidth=3)
        ax.scatter([Long_origin], [Lat_origin], label=f'Rotation Angle = {round(-1*rotation_angle, 2)} deg', color='black')

        # draw guiding lines
        ax.axline(xy1=(Lat_origin, 40), slope=(East_uvec_rot[1]/East_uvec_rot[0]), linestyle='--',linewidth=2, color='black', alpha=0.5)
        ax.axline(xy1=(Lat_origin, 115), slope=(East_uvec_rot[1] / East_uvec_rot[0]), linestyle='--', linewidth=2, color='black', alpha=0.5)

        # labels
        ax.axvline(x=Long_origin)
        ax.axhline(y=Long_origin)
        ax.text(Lat_origin, 0, 'Mag N', ha='center', va='center', fontsize=20)
        ax.text(2*Long_origin, Lat_origin, 'Mag E', ha='left', va='center', fontsize=20)
        ax.legend()

        # --- Convert the attitude data into pixel coordinates ---
        stl.prgMsg('Plotting Trajectory Data')

        # [0] Pair up the GLats/Glongs data into datapiars (GLat_val, GLong_val)
        xData = data_dicts_allsky[i]['GLats'][0]
        yData = data_dicts_allsky[i]['GLongs'][0]
        shape = np.shape(xData)
        latlong_pairs = [[] for i in range(shape[0])]
        counter = 0
        for arr1, arr2 in zip(xData, yData):
            latlong_pairs[counter] = [[val1, val2] for (val1, val2) in zip(arr1, arr2)]
            counter += 1

        # [1] For each ephemeris point, find the all sky indicies with the smallest distance between ephemeris lat/long
        high_flyer_test_statistic = np.array([np.sum(np.abs(np.array(latlong_pairs) - np.array([data_dicts_attitude[0]['Lat'][0][i], data_dicts_attitude[0]['Long'][0][i]])), axis=2) for i in range(len(data_dicts_attitude[0]['Epoch'][0]))])
        low_flyer_test_statistic = np.array([np.sum(np.abs(np.array(latlong_pairs) - np.array([data_dicts_attitude[1]['Lat'][0][i], data_dicts_attitude[1]['Long'][0][i]])), axis=2) for i in range(len(data_dicts_attitude[1]['Epoch'][0]))])

        # plot the rocket trajectories
        plt.savefig(rf'C:\Data\ACESII\science\auroral_coordinates\low\ASI_auroral_CoordinateS_check_{wavelength[i]}.png')
        stl.Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
Plot_ASI_auroral_coordinates_rotation_vec()

