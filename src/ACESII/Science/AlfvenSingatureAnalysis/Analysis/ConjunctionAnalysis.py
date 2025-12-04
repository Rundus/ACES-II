# --- Analysis.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Do an analysis between the high flyer and low flyer to determine the
# correlation between the HF wave-events and LF wave-events. Additionally do this for
# the particle data

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import datetime

import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *

start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
modifier = ''
inputPath_modifier_B = 'science\deltaB' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_E = 'science\deltaE'
inputPath_modifier_attitude = 'attitude'
inputPath_modifier_EEPAA = 'l2'

timeWindow = [
        datetime.datetime(2022,11,20,17,22,30,00), datetime.datetime(2022,11,20,17,27,30,00)
    ]
# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.class_var_func import CHAOS,InterpolateDataDict,lat_to_meter,calculateLong_to_meter
from matplotlib.gridspec import GridSpec

def ConjunctionAnalysis(rocketFolderPath):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()

    inputFile_attitude_low = glob(f'{rocketFolderPath}{inputPath_modifier_attitude}\{fliers[1]}{modifier}\*.cdf')[0]
    inputFile_attitude_high = glob(f'{rocketFolderPath}{inputPath_modifier_attitude}\{fliers[0]}{modifier}\*.cdf')[0]

    inputFile_deltaE_low = glob(f'{rocketFolderPath}{inputPath_modifier_E}\{fliers[1]}{modifier}\*.cdf')[0]

    inputFile_deltaB_high = glob(f'{rocketFolderPath}{inputPath_modifier_B}\{fliers[0]}{modifier}\*.cdf')[0]
    inputFile_deltaB_low = glob(f'{rocketFolderPath}{inputPath_modifier_B}\{fliers[1]}{modifier}\*.cdf')[0]

    # --- get the data from the Attitude files ---
    prgMsg(f'Loading data from {inputPath_modifier_attitude} Files')
    data_dict_attitude_low = loadDictFromFile(inputFile_attitude_low, {},reduceData=True,targetTimes=timeWindow)
    data_dict_attitude_low['Alt'][0] = data_dict_attitude_low['Alt'][0]/1000
    data_dict_attitude_high = loadDictFromFile(inputFile_attitude_high, {},reduceData=True,targetTimes=timeWindow)
    data_dict_attitude_high['Alt'][0] = data_dict_attitude_high['Alt'][0]/1000
    Done(start_time)

    # --- get the data from the E-Field files ---
    prgMsg(f'Loading data from {inputPath_modifier_E} Files')
    data_dict_deltaE_low = loadDictFromFile(inputFile_deltaE_low, {},reduceData=True,targetTimes=timeWindow)
    Done(start_time)

    # --- get the data from the B-Field files ---
    prgMsg(f'Loading data from {inputPath_modifier_B} Files')
    data_dict_deltaB_low = loadDictFromFile(inputFile_deltaB_low, {},reduceData=True,targetTimes=timeWindow)
    data_dict_deltaB_high = loadDictFromFile(inputFile_deltaB_high, {},reduceData=True,targetTimes=timeWindow)
    Done(start_time)

    #######################################################
    # --- USE CHAOS TO FIND LF'S PROJECTION POINT AT HF ---
    #######################################################
    prgMsg('Interpolating LF onto HF')
    # interpolate LF attitude data onto HF timebase
    data_dict_attitude_low = InterpolateDataDict(InputDataDict=data_dict_attitude_low,
                                                 InputEpochArray=np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_attitude_low['Epoch'][0]]),
                                                 wKeys=[],
                                                 targetEpochArray=np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_attitude_high['Epoch'][0]]))

    Done(start_time)

    prgMsg('Calculating CHAOS7 for LF')
    # Using CHAOS determine the Direction of B Along the LF's path

    LF_B = CHAOS(lat=data_dict_attitude_low['Latgd'][0],
          long=data_dict_attitude_low['Long'][0],
          alt=(data_dict_attitude_low['Alt'][0]),
          times=data_dict_attitude_low['Epoch'][0],
          )

    # Find the direction of B for all points in the LF's Flight
    B_unit = np.array([ vec/np.linalg.norm(vec) for vec in LF_B])
    Done(start_time)

    prgMsg('Triangulating Projection')

    # using triangulation, determine the latitude/longitude of the B-Field projection from LF


    # DELTA LAT
    # deltaLat_km = (ALT_HF - ALT_LF) * (B_north/|B_up|) * The sign of B_east determines the deltaLong --> need to SUBTRACT deltaLong from LF_long
    deltaLat_km = [(data_dict_attitude_high['Alt'][0][i] - data_dict_attitude_low['Alt'][0][i])*B_unit[i][1]/np.abs(B_unit[i][2]) for i in range(len(B_unit))]
    deltaLat_deg = [val/lat_to_meter for val in deltaLat_km]
    projection_Lat = np.array([data_dict_attitude_low['Latgd'][0][i] - deltaLat_deg[i] for i in range(len(B_unit))])

    # DELTA LONG
    # deltaLong_km = (ALT_HF - ALT_LF) * (B_east/|B_up|) * The sign of B_east determines the deltaLong --> need to SUBTRACT deltaLong from LF_long
    deltaLong_km = [(data_dict_attitude_high['Alt'][0][i] - data_dict_attitude_low['Alt'][0][i]) * B_unit[i][0] / np.abs(B_unit[i][2]) for i in range(len(B_unit))]
    deltaLong_deg = [deltaLong_km[i]/calculateLong_to_meter(data_dict_attitude_low['Latgd'][0][i]) for i in range(len(B_unit))]
    projection_Long = np.array([data_dict_attitude_low['Long'][0][i] - deltaLong_deg[i] for i in range(len(B_unit))])


    # Determine the absolute distance between Projection and HF
    from ACESII_code.class_var_func import GreatCircleDistance
    RadialDistance = np.array([GreatCircleDistance(lat1=data_dict_attitude_high['Latgd'][0][i],
                                                    lat2=projection_Lat[i],
                                                    long1=data_dict_attitude_high['Long'][0][i],
                                                    long2=projection_Long[i]) for i in range(len(B_unit))])
    Done(start_time)


    # Convert coordinates to geomagnetic coordinates


    # --- Make a Trajectory Comparison Plot ---

    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5)
    title = fig.suptitle('ACESII Conjunction')

    gs0 = GridSpec(nrows=2,ncols=1, height_ratios=[0.5, 0.5])


    # --- top Plot ---
    gsTop = gs0[0].subgridspec(nrows=1, ncols=2)

    # LatLong Plot
    axLatLong = fig.add_subplot(gsTop[:, 0])
    axLatLong.plot(data_dict_attitude_high['Long'][0], data_dict_attitude_high['Latgd'][0],color='red',label='High Flyer Trajectory')
    axLatLong.plot(projection_Long,projection_Lat,color='blue', label='LF $B_{geo}$ Projection')
    axLatLong.set_ylabel('Geographic Lat')
    axLatLong.set_xlabel('Geographic Long')
    axLatLong.legend()

    # Alt Lat Plot
    axAltLat = fig.add_subplot(gsTop[:, 1])
    axAltLat.plot(data_dict_attitude_high['Latgd'][0],data_dict_attitude_high['Alt'][0], color='red', label='High Flyer Trajectory')
    axAltLat.plot(projection_Lat,data_dict_attitude_high['Alt'][0], color='blue', label='LF $B_{geo}$ Projection')
    axAltLat.set_ylabel('Altitude [km]')
    axAltLat.set_xlabel('Geographic Lat')
    axAltLat.legend()

    # --- Bottom Plot ---
    gsBot = gs0[1].subgridspec(nrows=1, ncols=1)
    axDist = fig.add_subplot(gsBot[:, :])
    axDist.plot(data_dict_attitude_high['Epoch'][0],RadialDistance)

    plt.show()
















# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
ConjunctionAnalysis(rocketFolderPath)

