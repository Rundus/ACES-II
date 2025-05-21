# --- csv_to_cdf_trajectories.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Turn the .cdf files of the ACESII GPS data into cdf files


# --- --- --- --- ---
from src.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# 4 --> ACESII High Flyer
# 5 --> ACESII Low Flyer
wRocket = 4
outputData = True
wRow = 2 # row corresponding to where all the names of the variables are


# --- --- --- ---
# --- import ---
# --- --- --- ---
import numpy as np
from os import remove,path
from csv import reader
from copy import deepcopy
from scipy.interpolate import CubicSpline
setupPYCDF()
from spacepy import pycdf
from spacepy import coordinates as coord
coord.DEFAULTS.set_values(use_irbem=False, itol=5)  # maximum separation, in seconds, for which the coordinate transformations will not be recalculated. To force all transformations to use an exact transform for the time, set ``itol`` to zero.
from spacepy.time import Ticktock #used to determine the time I'm choosing the reference geomagentic field


def csv_to_cdf_trajectories(wRocket):


    print(f'Converting csv file\n')

    # --- get the path of the input file ---
    input_file = glob(rf'{DataPaths.ACES_data_folder}\\trajectories\\{ACESII.fliers[wRocket-4]}\\*Flight_trajectory_GPSdata.csv*')[0]

    # collect the csv data
    with open(input_file) as csvFile:
        csvAllData = [row for row in reader(csvFile)]

    varNames = csvAllData[wRow]
    csvData = np.array(csvAllData[wRow + 1:],dtype='float64').transpose()

    # --- Create output Data Dict ---
    exampleAttrs = {'DEPEND_0': None,'DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -9223372036854775808,'FORMAT': 'I5','UNITS': None,'VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'data','SCALETYP':'linear',}
    data_dict_output = {varNames[i] : [ csvData[i] ,  deepcopy(exampleAttrs)  ] for i in range(len(varNames))}


    # --- special modifications to the variables already in the csv file ---
    #  in a dictionary format [ {string name of var: {string name of attributeKey1: attributeval1, attributeKey2: attributeval2  }},{...},...  ]
    # Example: specialMods = [time: {'Units':'ns','FORMAT':'I5'} , Alt: {'Units':'km','FORMAT':'I5'} ]
    # Defaults to attributes = None if specialMods == [], except for VALIDMIN/VALIDMAX

    specialMods = {
        'Alt': {'UNITS': 'km', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'Alt'},
        'ECEFXPOS': {'UNITS': 'km', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'ECEFXPOS'},
        'ECEFXVEL': {'UNITS': 'm/s', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'ECEFXVEL'},
        'ECEFYPOS': {'UNITS': 'km', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'ECEFYPOS'},
        'ECEFYVEL': {'UNITS': 'm/s', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'ECEFYVEL'},
        'ECEFZPOS': {'UNITS': 'km', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'ECEFZPOS'},
        'ECEFZVEL': {'UNITS': 'm/s', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'ECEFZVEL'},
        'FlightTime': {'UNITS': 'Seconds from Launch', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'FlightTime'},
        'GPSTimemSecofweek': {'UNITS': 'Seconds', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'GPSTimemSecofweek'},
        'GPSWeek': {'UNITS': 'Weeks', 'LABLAXIS': 'GPSWeek'},
        'Lat': {'UNITS': 'deg', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'Lat'},
        'Long': {'UNITS': 'deg', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'Long'},
        'Epoch': {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -9223372036854775808,
                  'FORMAT': 'I5', 'UNITS': 'ns', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data',
                  'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000', 'TIME_SCALE': 'Terrestrial Time',
                  'REFERENCE_POSITION': 'Rotating Earth Geoid', 'SCALETYP': 'linear', 'LABLAXIS': 'Epoch'}
    }

    # --- Modify Variable Attributes in data dict output to the  ---
    for key, val in specialMods.items():
        if key in data_dict_output:
            for keyAttr, valAttr in specialMods[key].items():
                if keyAttr in data_dict_output[key][1]: # modify an existing attribute
                    data_dict_output[key][1][keyAttr] = valAttr
                else: # append a new attribute
                    data_dict_output[key][1] = {**data_dict_output[key][1], **{keyAttr:valAttr}}
        else: # Create a whole new variable
            data_dict_output = {**data_dict_output, **{key: [[], val ]} }


    # --- Load attitude data ---
    data_dict_attitude = stl.loadDictFromFile(glob(rf'{DataPaths.ACES_data_folder}\\attitude\\{ACESII.fliers[wRocket-4]}\\*.cdf*')[0])

    ###################################################################################
    # --- Create the Epoch Variable and limit it to the launch times of the rockets ---
    ###################################################################################

    # construct the Epoch variable
    Launch_time_tt2000 = ACESII.launch_T0_TT2000[wRocket-4]
    Launch_time_dt = pycdf.lib.tt2000_to_datetime(Launch_time_tt2000)
    Epoch = np.array([ pycdf.lib.tt2000_to_datetime(int(Launch_time_tt2000 + data_dict_output['FlightTime'][0][i]*(10**9))) for i in range(len(data_dict_output['FlightTime'][0])) ])
    data_dict_output['Epoch'][0] = Epoch

    # reduce all variables to just after launch
    launch_idx = np.abs(data_dict_output['Epoch'][0] - Launch_time_dt).argmin()
    for key, val in data_dict_output.items():
        data_dict_output[key][0] = data_dict_output[key][0][launch_idx:]
    stl.Done(start_time)



    ############################################
    # --- Convert to geomagnetic Coordinates ---
    ############################################
    stl.prgMsg('Calculating geomagnetic coordinates')

    geodetic = np.array([data_dict_output['Alt'][0],
                         data_dict_output['Lat'][0],
                         data_dict_output['Long'][0]]).transpose()

    # Get the times that the mission was launched in ISO datetime. Needed for geomagnetic coordinates
    ISOtime = [data_dict_output['Epoch'][0][j].isoformat() for j in range(len(data_dict_output['Epoch'][0]))]

    # --- Convert to geoMAG Coordinates ---

    cvals_GDZ = coord.Coords(geodetic, 'GDZ', 'sph')
    cvals_GDZ.ticks = Ticktock(ISOtime, 'ISO')
    cvals_GDZ_MAG = cvals_GDZ.convert('MAG', 'sph')

    # --- Make all the new variables ---
    geoMagLat = np.array(cvals_GDZ_MAG.lati)
    geoMagLong = np.array(cvals_GDZ_MAG.long)
    geoMagRadi = np.array(cvals_GDZ_MAG.radi)

    data_dict_output = {**data_dict_output, **{'geoMag_lat': [geoMagLat, {'LABLAXIS': 'Geomagnetic Latitude',
                                              'DEPEND_0': 'Epoch',
                                              'DEPEND_1': None,
                                              'DEPEND_2': None,
                                              'FILLVAL': -1e30,
                                              'FORMAT': 'E12.2',
                                              'UNITS': 'deg',
                                              'VALIDMIN': geoMagLat.min(),
                                              'VALIDMAX': geoMagLat.max(),
                                              'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
    data_dict_output = {**data_dict_output, **{'geoMag_long': [geoMagLong, {'LABLAXIS': 'Geomagnetic Longitude',
                                                            'DEPEND_0': 'Epoch',
                                                            'DEPEND_1': None,
                                                            'DEPEND_2': None,
                                                            'FILLVAL': -1e30,
                                                            'FORMAT': 'E12.2',
                                                            'UNITS': 'deg',
                                                            'VALIDMIN': geoMagLong.min(),
                                                            'VALIDMAX': geoMagLong.max(),
                                                            'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
    data_dict_output = {**data_dict_output, **{'geoMag_radi': [geoMagRadi, {'LABLAXIS': 'Geomagnetic Radius',
                                                            'DEPEND_0': 'Epoch',
                                                            'DEPEND_1': None,
                                                            'DEPEND_2': None,
                                                            'FILLVAL': -1e30,
                                                            'FORMAT': 'E12.2',
                                                            'UNITS': 'm',
                                                            'VALIDMIN': geoMagRadi.min(),
                                                            'VALIDMAX': geoMagRadi.max(),
                                                            'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
    stl.Done(start_time)

    ##########################################################
    # --- Interpolate Data onto Attitude Solution timebase ---
    ##########################################################
    stl.prgMsg('Interpolating on Attitude Timebase')
    Epoch_attitude_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_attitude['Epoch'][0]])
    Epoch_traj_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_output['Epoch'][0]])
    for key in data_dict_output.keys():
        cs = CubicSpline(Epoch_traj_tt2000,data_dict_output['key'][0])
        data_dict_output[key][0] = deepcopy(cs(Epoch_attitude_tt2000))
    stl.Done(start_time)

    if outputData:
        stl.prgMsg('outputting data')
        fileoutName = 'ACESII_36364_Flight_trajectory_GPSdata.cdf'
        outputPath = rf'{DataPaths.rocketFolderPath}\trajectories\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output, instrNam='ephemeris')
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
csv_to_cdf_trajectories(wRocket)

