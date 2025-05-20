# --- L2_EFI_to_auroral_coordinates.py ---
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

def Plot_allsky_auroral_coordinates_check():

    # Load AllSky data
    stl.prgMsg('Loading Allsky Data')
    data_dict_allSky5577 = stl.loadDictFromFile(glob(r'C:\Data\ACESII\all_sky\skibotn\5577\\*.cdf')[0])
    data_dict_allSky6300 = stl.loadDictFromFile(glob(r'C:\Data\ACESII\all_sky\skibotn\6300\\*.cdf')[0])
    data_dicts_attitude = [stl.loadDictFromFile(glob(f'C:\Data\ACESII\\attitude\\{DataPaths.fliers[0]}\\*.cdf*')[0]),stl.loadDictFromFile(glob(f'C:\Data\ACESII\\attitude\\{DataPaths.fliers[1]}\\*.cdf*')[0])]
    data_dict_auroral = stl.loadDictFromFile(glob(r'C:\Data\ACESII\science\auroral_coordinates\low\\ACESII_36364_E_Field_Auroral_Coordinates.cdf')[0])
    geoAlt = [data_dicts_attitude[0]['Alt'][0], data_dicts_attitude[1]['Alt'][0]]
    geoLat = [data_dicts_attitude[0]['Lat'][0], data_dicts_attitude[1]['Lat'][0]]
    geoLong = [data_dicts_attitude[0]['Long'][0], data_dicts_attitude[1]['Long'][0]]
    geoMagLat = [data_dicts_attitude[0]['Lat_geom'][0], data_dicts_attitude[1]['Lat_geom'][0]]
    stl.Done(start_time)

    # find the point (time, lat, long, alt) for when the Flyer's first hit the aurora vs when they left
    wpoints = 0 # 0 for HF 1 for LF

    image_set = [8,9]
    wImage = image_set[wpoints]

    enter_idx = [np.abs(data_dicts_attitude[0]['Epoch'][0]- dt.datetime(2022,11,20,17,24,10)).argmin(), np.abs(data_dicts_attitude[1]['Epoch'][0]- dt.datetime(2022, 11, 20, 17, 24, 48)).argmin()]
    leave_idx = [np.abs(data_dicts_attitude[0]['Epoch'][0]- dt.datetime(2022,11,20,17,28,50)).argmin(),np.abs(data_dicts_attitude[1]['Epoch'][0]- dt.datetime(2022, 11, 20, 17, 26, 24)).argmin()]

    # choose an entrance origin
    origin_enter = [data_dicts_attitude[wpoints]['Long'][0][enter_idx[wpoints]],data_dicts_attitude[wpoints]['Lat'][0][enter_idx[wpoints]]]
    origin_leave = [data_dicts_attitude[wpoints]['Long'][0][leave_idx[wpoints]]],data_dicts_attitude[wpoints]['Lat'][0][leave_idx[wpoints]]

    # --- --- --- --- --- --
    # --- BigAllSky plot ---
    # --- --- --- --- --- --
    projProjection = ccrs.Orthographic(central_longitude=15, central_latitude=70)
    projTransform = ccrs.PlateCarree()

    for i in range(2):
        # --- PLOT MAP OF NORWAY ---
        fig, axBigAllSky = plt.subplots(1,subplot_kw=dict(projection=projProjection))
        fig.set_figwidth(figure_width)
        fig.set_figheight(figure_height)

        # gridlines
        gl = axBigAllSky.gridlines(draw_labels=True, linewidth=BigAllSky_GridSize,
                                   alpha=0.4,
                                   linestyle='--',
                                   color='black')
        gl.xlabel_style = {'size': BigAllSky_tickLabelSize, 'color': 'black', 'weight': 'bold'}
        gl.ylabel_style = {'size': BigAllSky_tickLabelSize, 'color': 'black', 'weight': 'bold'}
        gl.top_labels = False

        # extent of map
        axBigAllSky.set_extent([lonW, lonE, latS, latN])  # controls lat/long axes display

        # coastlines
        axBigAllSky.coastlines(resolution=res, color='black',  alpha=1,linewidth=BigAllSky_costLineSize)  # adds coastlines with resolution

        if i == 0:
            #--- Plot the Big AllSky image ---
            cmapBigAllSky = axBigAllSky.pcolormesh(data_dict_allSky5577['GLongs'][0], data_dict_allSky5577['GLats'][0], data_dict_allSky5577['AllSkyImages'][0][wImage],
                                                   cmap=cmapColor,
                                                   transform=projTransform,
                                                   vmin=cbarVmin,
                                                   vmax=cbarVmax,
                                                   alpha=1)
            BigAllSky_outputPath = r'C:\Data\ACESII\science\auroral_coordinates\low\5570.png'
            fig.suptitle('Skibotn 5577$\AA$ - 150 km\n' + data_dict_allSky5577['Epoch'][0][wImage].strftime("%Y-%B-%d %H:%M:%S") + ' UTC',fontsize=BigAllSky_TitleSize,weight='bold')
        elif i == 1:
            cmapBigAllSky = axBigAllSky.pcolormesh(data_dict_allSky6300['GLongs'][0], data_dict_allSky6300['GLats'][0],
                                                   data_dict_allSky6300['AllSkyImages'][0][wImage],
                                                   cmap=cmapColor,
                                                   transform=projTransform,
                                                   vmin=cbarVmin,
                                                   vmax=cbarVmax,
                                                   alpha=1)
            BigAllSky_outputPath = r'C:\Data\ACESII\science\auroral_coordinates\low\6300.png'
            fig.suptitle('Skibotn 6300$\AA$ - 250 km\n' + data_dict_allSky6300['Epoch'][0][wImage].strftime("%Y-%B-%d %H:%M:%S")+ ' UTC',fontsize=BigAllSky_TitleSize,weight='bold')
        axBigAllSky.set_facecolor(faceColorChoice)

        # --- plot the rocket trajectory data on the large AllSky plot ---
        axBigAllSky.plot(geoLong[0], geoLat[0], color=trajColors[0], transform=projTransform,linewidth=BigAllSky_lineThickness) # High
        axBigAllSky.plot(geoLong[1], geoLat[1], color=trajColors[1], transform=projTransform,linewidth=BigAllSky_lineThickness) # Low



        ##########################################################################
        # --- Get the transformation matricies from ENU to Auroral Coordiantes ---
        ##########################################################################

        # [1] Get ENU to ECEF
        print(origin_enter)
        matrix_ENUtoECEF = stl.ENUtoECEF(Lat=origin_enter[1], Long=origin_enter[0])

        # [2] Get ECEF to Field-Aligned
        # Get the Data
        B_model = stl.CHAOS(lat=[origin_enter[1]],
                            long=[origin_enter[0]],
                            alt=[150], # in km
                            times=[dt.datetime(2022, 11, 20, 17, 24, 48)])  # CHAOS in ENU coordinates

        # --- Convert B-Data to GEO (ECEF) XYZ coordinates ---
        B_CHAOS_ECEF = np.matmul(matrix_ENUtoECEF, B_model[0])

        # --- determine the Payload's Position Vector in GEO (ECEF) coordinate XYZ ---
        R_REF = 6371.2  # earth Radius in km
        Radius = 150 + R_REF
        coLatRad = np.radians(90 - origin_enter[1])
        LongRad = np.radians(origin_enter[0])
        Rsc = np.array(
            [Radius * np.sin(coLatRad) * np.cos(LongRad),
             Radius * np.sin(coLatRad) * np.sin(LongRad),
             Radius * np.cos(coLatRad)])
        stl.Done(start_time)

        # --- calculate Field Aligned unit vectors over the duration of the flight ---
        stl.prgMsg('Converting to Field Aligned Coordinates')

        # pHat comes from the CHAOS model direction of B in GEO
        pHat = B_CHAOS_ECEF/ np.linalg.norm(B_CHAOS_ECEF)

        # e-hat comes from the cross of pHat and the Rocket's radius vector (in geomagnetic coordinates)
        eHat = np.cross(pHat, Rsc) / np.linalg.norm(np.cross(pHat, Rsc))

        # rHat comes from the cross of eHat and pHat
        rHat = np.cross(eHat, pHat)

        # form the transformation matrix FROM ECEF TO FIELD ALIGNED
        FAC_transform = np.array([eHat, pHat, rHat])

        # form the final transformation matrix
        R_1 = matrix_ENUtoECEF
        R_2 = FAC_transform
        R_3 = stl.Rz(-9)
        # transform_total = np.matmul(R_3,np.matmul(R_2,R_1))
        transform_total = R_3

        E_rotated = np.matmul(transform_total,np.array([1, 0, 0]))
        N_rotated = np.matmul(transform_total, np.array([0, 1, 0]))

        # get the angle between the East direction and the new vector direction
        angle_E = np.degrees(np.arccos(np.dot(E_rotated[0:2], [1, 0])/(np.linalg.norm(E_rotated[0:2]))))
        angle_N = np.degrees(np.arccos(np.dot(N_rotated[0:2], [0, 1]) / (np.linalg.norm(N_rotated[0:2]))))

        # --- Plot the auroral coordinate transform ---
        geographic_angles = [0,0+90]
        # EASTWARD
        X1 = origin_enter[0]
        Y1 = origin_enter[1]
        U1 = 1  # directional vector
        V1 = 0
        X = np.array([X1])
        Y = np.array([Y1])
        U = np.array([U1])
        V = np.array([V1])
        # q1=axBigAllSky.quiver(X, Y, U, V, color='black', transform=projTransform, scale=10, angles=geographic_angles[0])

        # NORTHWARD
        X2 = origin_enter[0]
        Y2 = origin_enter[1]
        U2 = 0 # directional vector
        V2 = 1
        X = np.array([X2])
        Y = np.array([Y2])
        U = np.array([U2])
        V = np.array([V2])
        # q2 = axBigAllSky.quiver(X, Y, U, V, color='black', transform=projTransform, scale=5,angles=geographic_angles[1])

        # make a fake rotation to test the orthognality of the projection on small distances
        auroral_Angle = data_dict_auroral['rotation_Angle'][0]
        auroral_angles = [-auroral_Angle,-auroral_Angle+90]
        east_new = np.matmul(stl.Rz(45),[1,0,0])
        X3 = origin_enter[0]
        Y3 = origin_enter[1]
        U3 = east_new[0]   # directional vector
        V3 = east_new[1]
        X = np.array([X3])
        Y = np.array([Y3])
        U = np.array([U3])
        V = np.array([V3])
        axBigAllSky.quiver(X, Y, U, V, color='red', transform=projTransform, scale=6, angles=auroral_angles[0], label=f'{auroral_Angle} deg')
        axBigAllSky.legend()

        # make a fake rotation to test the orthognality of the projection on small distances
        north_new = np.matmul(stl.Rz(45), [0, 1, 0])
        X4 = origin_enter[0]
        Y4 = origin_enter[1]
        U4 = north_new[0] # directional vector
        V4 = north_new[1]

        X = np.array([X4])
        Y = np.array([Y4])
        U = np.array([U4])
        V = np.array([V4])
        axBigAllSky.quiver(X, Y, U, V, color='green', transform=projTransform, scale=10,
                           angles=auroral_angles[1])
        plt.tight_layout()
        plt.savefig(BigAllSky_outputPath)
        stl.Done(start_time)

#
# def Plot_allsky_auroral_coordinates_check():
#
#     #######################
#     # --- LOAD THE DATA ---
#     #######################
#     stl.prgMsg(f'Loading data')
#     data_dict_allSky5577 = stl.loadDictFromFile(glob(r'C:\Data\ACESII\all_sky\skibotn\5577\\*.cdf')[0])
#     data_dict_allSky6300 = stl.loadDictFromFile(glob(r'C:\Data\ACESII\all_sky\skibotn\6300\\*.cdf')[0])
#     data_dicts_allsky = [data_dict_allSky5577,data_dict_allSky6300]
#     data_dict_EFI_auroral = stl.loadDictFromFile(glob(r'C:\Data\ACESII\science\auroral_coordinates\low\*.cdf*')[0])
#     attitudeFolderPath = f'{DataPaths.ACES_data_folder}\\attitude'
#     inputFilesTraj = [glob(f'{attitudeFolderPath}\\{DataPaths.fliers[0]}\\*.cdf*')[0],
#                       glob(f'{attitudeFolderPath}\\{DataPaths.fliers[1]}\\*.cdf*')[0]]
#     data_dicts_attitude = [stl.loadDictFromFile(inputFilesTraj[0]), stl.loadDictFromFile(inputFilesTraj[1])]
#     stl.Done(start_time)
#
#     for i in range(2):
#         # --- PLOT ---
#         fig, ax = plt.subplots(constrained_layout=True)
#         fig.set_figwidth(figure_width)
#         fig.set_figheight(figure_height)
#
#         # all sky image
#         ax.set_title(f'Skibotn {wavelength[i]}$\AA$\n' + data_dicts_allsky[i]['Epoch'][0][wImage].strftime("%Y-%B-%d %H:%M:%S") + ' UTC',weight='bold')
#         ax.imshow(data_dicts_allsky[i]['AllSkyImages'][0][wImage], cmap=cmapColor, vmin=cbarVmin, vmax=cbarVmax)
#
#         # all center center coordinates
#         Lat_origin = len(data_dict_allSky5577['AllSkyImages'][0][wImage]) / 2
#         Long_origin = len(data_dict_allSky5577['AllSkyImages'][0][wImage]) / 2
#         offset = 20
#         ax.set_ylim(-1*offset, Lat_origin*2 + offset)
#         ax.set_xlim(-1*offset, Lat_origin * 2 + offset)
#         plt.gca().invert_yaxis()
#
#         # Units vectors
#         mag = 300
#         East_uvec = mag * np.array([1, 0])
#         North_uvec = mag * np.array([0, -1])
#
#         def rotate(angle, vec):
#             rad = np.radians(angle)  # the 90deg is to account for the imshow rotation
#             matrix = np.array([[np.cos(rad), -1 * np.sin(rad)], [np.sin(rad), np.cos(rad)]])
#             return np.matmul(matrix, vec)
#
#         rotation_angle = -1 * data_dict_EFI_auroral['rotation_Angle'][0][0]
#         East_uvec_rot = rotate(rotation_angle, East_uvec)
#         North_uvec_rot = rotate(rotation_angle, North_uvec)
#         east_vector = [[Long_origin, East_uvec_rot[0] + Long_origin], [Lat_origin, East_uvec_rot[1] + Lat_origin]]
#         north_Vector = [[Long_origin, North_uvec_rot[0] + Long_origin], [Lat_origin, North_uvec_rot[1] + Lat_origin]]
#         ax.plot(north_Vector[0], north_Vector[1], color='purple', linewidth=3)
#         ax.scatter([Long_origin], [Lat_origin], label=f'Rotation Angle = {round(-1*rotation_angle, 2)} deg', color='black')
#
#         # draw guiding lines
#         ax.axline(xy1=(Lat_origin, 40), slope=(East_uvec_rot[1]/East_uvec_rot[0]), linestyle='--',linewidth=2, color='black', alpha=0.5)
#         ax.axline(xy1=(Lat_origin, 115), slope=(East_uvec_rot[1] / East_uvec_rot[0]), linestyle='--', linewidth=2, color='black', alpha=0.5)
#
#         # labels
#         ax.axvline(x=Long_origin)
#         ax.axhline(y=Long_origin)
#         ax.text(Lat_origin, 0, 'Mag N', ha='center', va='center', fontsize=20)
#         ax.text(2*Long_origin, Lat_origin, 'Mag E', ha='left', va='center', fontsize=20)
#         ax.legend()
#
#         # --- Convert the attitude data into pixel coordinates ---
#         stl.prgMsg('Plotting Trajectory Data')
#
#         # [0] Pair up the GLats/Glongs data into datapiars (GLat_val, GLong_val)
#         xData = data_dicts_allsky[i]['GLats'][0]
#         yData = data_dicts_allsky[i]['GLongs'][0]
#         shape = np.shape(xData)
#         latlong_pairs = [[] for i in range(shape[0])]
#         counter = 0
#         for arr1, arr2 in zip(xData, yData):
#             latlong_pairs[counter] = [[val1, val2] for (val1, val2) in zip(arr1, arr2)]
#             counter += 1
#
#         # [1] For each ephemeris point, find the all sky indicies with the smallest distance between ephemeris lat/long
#         high_flyer_test_statistic = np.array([np.sum(np.abs(np.array(latlong_pairs) - np.array([data_dicts_attitude[0]['Lat'][0][i], data_dicts_attitude[0]['Long'][0][i]])), axis=2) for i in range(len(data_dicts_attitude[0]['Epoch'][0]))])
#         low_flyer_test_statistic = np.array([np.sum(np.abs(np.array(latlong_pairs) - np.array([data_dicts_attitude[1]['Lat'][0][i], data_dicts_attitude[1]['Long'][0][i]])), axis=2) for i in range(len(data_dicts_attitude[1]['Epoch'][0]))])
#
#         # plot the rocket trajectories
#         plt.savefig(rf'C:\Data\ACESII\science\auroral_coordinates\low\ASI_auroral_CoordinateS_check_{wavelength[i]}.png')
#         stl.Done(start_time)
#

# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
Plot_allsky_auroral_coordinates_check()

