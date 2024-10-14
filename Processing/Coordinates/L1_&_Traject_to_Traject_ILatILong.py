# --- L1_&_Traject_to_Traject_ILatILong.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Takes in the ACESII trajectory data (in geodetic and geomagnetic) as well
# as the ESA data (individual Instruments). Interpolates the Geomag and Alt data then produces the Ionospheric Projected Lattitude and
# Longitude


# TODO: This file wont work if there's multple ESA datafiles with "eepaa" or "iepaa" etc in the L1 folder! Must only be one.Should fix this

# --- --- --- --- ---
import time
from ACESII_code.class_var_func import Done, setupPYCDF,prgMsg
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

modifier = '_ILatILong'

targetProjectionAltitude = 100 # altitude you want to project B (in km). Should be ~100km

# select which files to convert
wInstr = 'eepaa' # valid inputs are strings of eepaa, iepaa and leesa

# For the Kilometers conversion, sets the reference point to Andoya
useAndoya = True

# --- output data? ---
outputDataFile = True


# --- --- --- ---
# --- import ---
# --- --- --- ---
import numpy as np
import pyIGRF
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from os import remove,path
from glob import glob
from ACESII_code.data_paths import fliers, ACES_data_folder, ACES_csv_trajectories, TRICE_data_folder, ACES_L0_files, TRICE_L0_files
from ACESII_code.class_var_func import color, L1_TRICE_Quick
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
setupPYCDF()
from spacepy import pycdf
from spacepy import coordinates as coord
coord.DEFAULTS.set_values(use_irbem=False, itol=5)  # maximum separation, in seconds, for which the coordinate transformations will not be recalculated. To force all transformations to use an exact transform for the time, set ``itol`` to zero.
from spacepy.time import Ticktock #used to determine the time I'm choosing the reference geomagentic field

print(color.BOLD + color.CYAN + 'L1_&_Traject_to_Traject_ILatILong.py' + color.END + color.END)
def Trajectory_to_ESA_ILatILong(wInstr, rocketFolderPath):

    if wInstr == 'leesa':
        rangelen = 1
    else:
        rangelen = 2

    # --- ACES II Flight/Integration Data ---
    rocketAttrs,b,c = ACES_mission_dicts()
    globalAttrsMod = rocketAttrs.globalAttributes[0]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L1'

    # Set the paths for the file names
    TrajectoryFiles = []
    L1Files = []
    ILatILongFiles = []
    L1_names = []
    ILatILong_names = []
    dataFile_name = []
    fileoutName = []

    for i in range(rangelen):
        TrajectoryFiles.append(glob(f"{rocketFolderPath}trajectories\{fliers[i]}\*GPSData.cdf*"))
        L1Files.append(glob(f'{rocketFolderPath}L1\{fliers[i]}\*{wInstr}*'))
        ILatILongFiles.append(glob(f'{rocketFolderPath}trajectories\{fliers[i]}\*{wInstr}*'))

        L1_names.append([ifile.replace(f'{rocketFolderPath}L1\{fliers[i]}\\', '') for ifile in L1Files[i]])
        ILatILong_names.append([ofile.replace(f'{rocketFolderPath}trajectories\{fliers[i]}\\', '') for ofile in ILatILongFiles[i]])

        dataFile_name.append(L1Files[i][0].replace(f'{rocketFolderPath}trajectories\{fliers[i]}\\', ''))

        fileoutName.append(TrajectoryFiles[i][0].replace('GPSdata', 'ILat_ILong').replace(f'{rocketFolderPath}trajectories\{fliers[i]}\\',''))

    if wInstr == 'lp':
        print(color.RED + 'Cannot processes Langmuir Probe File' + color.END)
    else:
        print(color.UNDERLINE + f'Processing data for {wInstr.upper()} instrument' + color.END)

        ######################
        # --- LOAD IN DATA ---
        ######################
        prgMsg('Loading data from L1 and GPS Data')
        data_dicts = []
        data_dicts_traj = []

        for i in range(rangelen):
            data_dict = {}
            with pycdf.CDF(L1Files[i][0]) as L1DataFile:
                for key, val in L1DataFile.items():
                    data_dict = {**data_dict, **{key: [L1DataFile[key][...] , {key:val for key, val in L1DataFile[key].attrs.items()  }  ]  }  }

            data_dict['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_esa'][0][i]) for i in range(len(data_dict['Epoch_esa'][0]))])

            data_dicts.append(data_dict)

            data_dict_traj = {}

            with pycdf.CDF(TrajectoryFiles[i][0]) as TrajDataFile:
                for key, val in TrajDataFile.items():
                    data_dict_traj = {**data_dict_traj, **{key: [TrajDataFile[key][...], {key: val for key, val in TrajDataFile[key].attrs.items()}]}}

            data_dict_traj['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_traj['Epoch'][0][i]) for i in range(len(data_dict_traj['Epoch'][0]))])

            data_dicts_traj.append(data_dict_traj)

        Done(start_time)


        ##############################################
        # --- DownSample GPS Alt, Lat, Long, Epoch ---
        ##############################################

        prgMsg('Downsampling GPS data')
        # Using the Trajectory data's Epoch, downsample the traj data to align with the ESA data
        AltDS = [[],[]]
        LatDS = [[],[]]
        LongDS = [[],[]]
        EpochDS = [[],[]]

        for j in range(rangelen):

            for i in range(len(data_dicts[j]['Epoch_esa'][0])):
                targetIndex = np.abs(data_dicts_traj[j]['Epoch'][0] - data_dicts[j]['Epoch_esa'][0][i]).argmin()
                AltDS[j].append(data_dicts_traj[j]['Alt'][0][targetIndex])
                LatDS[j].append(data_dicts_traj[j]['Lat'][0][targetIndex])
                LongDS[j].append(data_dicts_traj[j]['Long'][0][targetIndex])
                EpochDS[j].append(data_dicts_traj[j]['Epoch'][0][targetIndex])

            AltDS[j] = np.array(AltDS[j])
            LatDS[j] = np.array(LatDS[j])
            LongDS[j] = np.array(LongDS[j])
            EpochDS[j] = np.array(EpochDS[j])

        Done(start_time)

        ################################
        # --- Calculate IGRF B-Field ---
        ################################
        prgMsg('Getting IGRF Field')
        # -- Output order forpyIGRF.igrf_value ---
        # [0] Declination (+ E | - W)
        # [1] Inclination (+ D | - U), should be ~78deg for what we're doing
        # [2] Horizontal Intensity
        # [3] North Comp (+ N | - S)
        # [4] East Comp (+ E | - W)
        # [5] Vertical Comp (+ D | - U)
        # [6] Total Field

        IGRF = [[],[]]
        date = 2022 + 323 / 365  # Corresponds to 11/20/2022
        for i in range(rangelen):
            for j in range(len(data_dicts[i]['Epoch_esa'][0])):
                IGRF[i].append(pyIGRF.igrf_value(LatDS[i][j], LongDS[i][j], AltDS[i][j], date))

        Done(start_time)

        #################################
        # --- I-LAT I-LONG PROJECTION ---
        #################################

        prgMsg('Projecting B-Fields')

        lat_to_meter = 111.319488 # 1 deg latitude to kilometers on Earth
        def long_to_meter(lat):
            return 111.319488 * math.cos(lat*math.pi/180)

        intersectionPoints = [[],[]]
        BFieldDirNorm = [[],[]]
        BFieldLoc = [[],[]]

        # Perform Triangulation Projection
        for j in range(rangelen):

            for i in range(len(data_dicts[j]['Epoch_esa'][0])):
                #Coordiantes reported in (Long (x) , Lat (y), Alt (z))
                vLoc = np.array([LongDS[j][i], LatDS[j][i],AltDS[j][i]]) # IGRF vector Location, should be rocket coordinates
                vDir = np.array([IGRF[j][i][4], IGRF[j][i][3],-1*IGRF[j][i][5]]) # IGRF vector Direction. -1 added in third elemental due to down being positive in IGRF given
                vDirNorm = vDir / np.linalg.norm(vDir) # Normalize IGRF to get its direction only. This will make t larger, but that's fine
                BFieldDirNorm[j].append(vDirNorm)
                BFieldLoc[j].append(vLoc)

                # Determine the Delta-Latitude and Longitutde
                Theta_dec = IGRF[j][i][0]
                Theta_in = IGRF[j][i][1]
                h = vLoc[2] - targetProjectionAltitude
                deltaLat = (1/lat_to_meter) * h*math.tan((math.pi/180) *(90 - Theta_in))
                deltaLong = (1/long_to_meter(LatDS[j][i])) * h*math.tan((math.pi/180) * (Theta_dec))
                intersectionPoints[j].append([vLoc[0]+deltaLong, vLoc[1] + deltaLat, 0])

        Done(start_time)

        ############################################
        # --- CONVERT TO GEOMAGNETIC COORDINATES ---
        ############################################
        prgMsg('Converting Coordinates')

        # Trajectory Data
        geomagAlt = []
        geomagLat = []
        geomagLong = []


        for i in range(rangelen):
            geodeticPos = np.array([AltDS[i],LatDS[i],LongDS[i]]).transpose()
            ISOtime = [pycdf.lib.tt2000_to_datetime(EpochDS[i][j]).isoformat() for j in range(len(EpochDS[i]))]
            cvals_GDZ = coord.Coords(geodeticPos, 'GDZ', 'sph')
            cvals_GDZ.ticks = Ticktock(ISOtime, 'ISO')
            cvals_GDZ_MAG = cvals_GDZ.convert('MAG', 'sph')
            geomagAlt.append(cvals_GDZ_MAG.radi)
            geomagLat.append(cvals_GDZ_MAG.lati)
            geomagLong.append(cvals_GDZ_MAG.long)

        # Intersections
        geodeticLongIntersects = [[],[]]
        geodeticLatIntersects = [[],[]]
        geodeticAltIntersects = [[],[]]
        geoMagFootPrint_lat = [[],[]]
        geoMagFootPrint_long = [[],[]]

        for j in range(rangelen):
            geodeticLongIntersects[j] = np.array([intersectionPoints[j][i][0] for i in range(len(intersectionPoints[j]))]) # Long
            geodeticLatIntersects[j] = np.array([intersectionPoints[j][i][1] for i in range(len(intersectionPoints[j]))])  # Lat
            geodeticAltIntersects[j] = np.array([intersectionPoints[j][i][2] for i in range(len(intersectionPoints[j]))])  # Alt
            geodetic = np.array([geodeticAltIntersects[j], geodeticLatIntersects[j], geodeticLongIntersects[j]]).transpose()

            ISOtime = [ pycdf.lib.tt2000_to_datetime(EpochDS[j][i]).isoformat() for i in range(len(EpochDS[j]))]
            cvals_GDZ = coord.Coords(geodetic, 'GDZ', 'sph')
            cvals_GDZ.ticks = Ticktock(ISOtime, 'ISO')
            cvals_GDZ_MAG_intersects = cvals_GDZ.convert('MAG', 'sph')

            geoMagFootPrint_lat[j] = cvals_GDZ_MAG_intersects.lati
            geoMagFootPrint_long[j] = cvals_GDZ_MAG_intersects.long

        Done(start_time)

        #######################
        # --- CONVERT TO KM ---
        #######################

        # Determine the distance in km from some lat/long reference point
        if useAndoya:
            refLat = rocketAttrs.Andoya_Space_Lat_Long[0]
            refLong = rocketAttrs.Andoya_Space_Lat_Long[1]
        else:
            refLat = 0
            refLong = 0

        geodeticLatIntersects_km = [[], []]
        geodeticLongIntersects_km = [[], []]
        LatDS_km = [[], []]
        LongDS_km = [[], []]


        for i in range(rangelen):
            # Convert lat/long data to KM
            for j in range(len(geodeticLatIntersects[i])):
                geodeticLatIntersects_km[i].append(lat_to_meter * (geodeticLatIntersects[i][j] - refLat))
                geodeticLongIntersects_km[i].append(long_to_meter(refLat) * (geodeticLongIntersects[i][j] - refLong))
                LatDS_km[i].append(lat_to_meter * (LatDS[i][j] - refLat))
                LongDS_km[i].append(long_to_meter(refLat) * LongDS[i][j] - long_to_meter(refLat) * refLong)

        # Determine the magntitude distance from the launch point
        distanceFromLaunchPoint = [
            [np.sqrt((LatDS_km[0][i])**2 + (LongDS_km[0][i])**2) for i in range(len(LatDS_km[0]))],
            [np.sqrt((LatDS_km[1][i]) ** 2 + (LongDS_km[1][i]) ** 2) for i in range(len(LatDS_km[1]))]
        ]


        ###############################
        # --- CONSTRUCT OUTPUT DATA ---
        ###############################

        # Remove the particular instrument data. This is to avoid using this ESA data which may be obsolete later...although it may not matter the more I think about it

        for i in range(rangelen):

            popThese = []
            for key, val in data_dicts[i].items():
                if key not in [wInstr,'Pitch_Angle','Energy','Epoch_esa']:
                    popThese.append(key)

            for key in popThese:
                data_dicts[i].pop(key)




            dataNEW = {
                'geoAlt':[np.array(AltDS[i]),'geodetic_Altitude'],
                'geoLat':[np.array(LatDS[i]),'geodetic_Latitude'],
                'geomagAlt':[np.array(geomagAlt[i]),'geomagnetic_Altitude'],
                'geomagLat':[np.array(geomagLat[i]),'geomagnetic_Latitude'],
                'geomagLong':[np.array(geomagLong[i]),'geomagnetic_Longitude'],
                'geoLong':[np.array(LongDS[i]),'geodetic_Longitude'],
                'geoLat_km': [np.array(LatDS_km[i]), 'geodetic_Latitude_km'],
                'geoLong_km': [np.array(LongDS_km[i]), 'geodetic_Longitude_km'],
                'geomagILat': [np.array(geoMagFootPrint_lat[i]),'mapped_Ionospheric_geomagetic_latitude'],
                'geomagILong': [np.array(geoMagFootPrint_long[i]),'mapped_Ionospheric_geomagnetic_longitude'],
                'geoILong': [np.array(geodeticLongIntersects[i]),'mapped_Ionospheric_geodetic_longitude'],
                'geoILat': [np.array(geodeticLatIntersects[i]),'mapped_Ionospheric_geodetic_latitude'],
                'geoILat_km': [np.array(geodeticLatIntersects_km[i]),'mapped_Ionospheric_geodetic_latitude_km'],
                'geoILong_km': [np.array(geodeticLatIntersects_km[i]),'mapped_Ionospheric_geodetic_longitude_km'],
                'distanceFromLaunchPoint_km': [np.array(distanceFromLaunchPoint[i]), 'distance_From_Launch_Point_km']
            }

            # Add the new data to the data_dict
            for key, val in dataNEW.items():
                if '_km' in key or 'Alt' in key:
                    if 'geomagAlt' in key:
                        unitLabel = 'Re'
                    else:
                        unitLabel = 'km'
                else:
                    unitLabel = 'deg'

                data_dicts[i] = {**data_dicts[i], **{key: [val[0], {'LABLAXIS': val[1],
                                                                 'DEPEND_0': 'Epoch_esa', 'DEPEND_1': None,
                                                                 'DEPEND_2': None,
                                                                 'FILLVAL': -1e30, 'FORMAT': 'E12.2',
                                                                 'UNITS': unitLabel,
                                                                 'VALIDMIN': val[0].min(), 'VALIDMAX': val[0].max(),
                                                                 'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}


        if outputDataFile:
            prgMsg('Writing out Data')

            for i in range(rangelen):

                #####################
                # --- Output Data ---
                #####################

                outputPath = f'{rocketFolderPath}Trajectories\{fliers[i]}\\{fileoutName[i]}'

                # --- delete output file if it already exists ---
                if path.exists(outputPath):
                    remove(outputPath)
                pycdf.lib.set_backward(False)

                # --- open the output file ---
                with pycdf.CDF(outputPath, '') as cdfFile:
                    cdfFile.readonly(False)

                    # --- write out global attributes ---
                    globalAttrsMod = rocketAttrs.globalAttributes[i]
                    ModelData = L1_TRICE_Quick(i)
                    inputGlobDic = ModelData.cdfFile.globalattsget()

                    for key, val in inputGlobDic.items():
                        if key == 'Descriptor':
                            globalAttrsMod[key] = 'None'
                        if key in globalAttrsMod:
                            cdfFile.attrs[key] = globalAttrsMod[key]
                        else:
                            cdfFile.attrs[key] = val

                    # --- WRITE OUT DATA ---
                    for varKey, varVal in data_dicts[i].items():
                        if 'Epoch' in varKey:  # epoch data
                            cdfFile.new(varKey, data=varVal[0], type=33)
                        else:  # other data
                            cdfFile.new(varKey, data=varVal[0])

                        # --- Write out the attributes and variable info ---
                        for attrKey, attrVal in data_dicts[i][varKey][1].items():
                            if attrKey == 'VALIDMIN':
                                cdfFile[varKey].attrs[attrKey] = varVal[0].min()
                            elif attrKey == 'VALIDMAX':
                                cdfFile[varKey].attrs[attrKey] = varVal[0].max()
                            elif attrVal != None:
                                cdfFile[varKey].attrs[attrKey] = attrVal

            Done(start_time)






# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder

if len(glob(f'{rocketFolderPath}L1\{fliers[0]}\*.cdf')) == 0 :
    print(color.RED + 'There are no .cdf files in the specified directory:' + color.END)
    print(f'{rocketFolderPath}L1\{fliers[0]}\*.cdf')
elif len(glob(f'{rocketFolderPath}L1\{fliers[1]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory:' + color.END)
    print(f'{rocketFolderPath}L1\{fliers[1]}\*.cdf')
else:
    Trajectory_to_ESA_ILatILong(wInstr, rocketFolderPath)

