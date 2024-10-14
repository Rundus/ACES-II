# --- csv_to_cdf_trajectories.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Turn the .cdf files of the ACESII GPS data into cdf files


# --- --- --- --- ---
import time
from ACESII_code.class_var_func import Done, setupPYCDF,prgMsg
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

wRocket = 0
# 0 --> ACESII High Flyer
# 1 --> ACESII Low Flyer
# 2 --> TRICEII High Flyer
# 3 --> TRICEII Low Flyer


Directories = [r'D:\Data\ACESII\trajectories\high\\',r'D:\Data\ACESII\trajectories\low\\']
Filenames = [r'Bounds_36359_Flight_trajectory_GPSdata.csv',r'Bounds_36364_Flight_trajectory_GPSdata.csv']
inputDirectory = Directories[wRocket]
inputFileName = Filenames[wRocket]
outputDirectory = inputDirectory
outputFileName = inputFileName.replace('.csv', '.cdf')

wRow = 2 # row corresponding to where all the names of the variables are



# --- special modifications to the variables already in the csv file ---
#  in a dictionary format [ {string name of var: {string name of attributeKey1: attributeval1, attributeKey2: attributeval2  }},{...},...  ]
# Example: specialMods = [time: {'Units':'ns','FORMAT':'I5'} , Alt: {'Units':'km','FORMAT':'I5'} ]
# Defaults to attributes = None if specialMods == [], except for VALIDMIN/VALIDMAX


specialMods = {
    'Alt':               {'UNITS':'km', 'DEPEND_0':'Epoch','LABLAXIS':'Alt'},
    'ECEFXPOS':          {'UNITS':'km', 'DEPEND_0':'Epoch','LABLAXIS':'ECEFXPOS'},
    'ECEFXVEL':          {'UNITS':'m/s','DEPEND_0':'Epoch','LABLAXIS':'ECEFXVEL'},
    'ECEFYPOS':          {'UNITS':'km', 'DEPEND_0':'Epoch','LABLAXIS':'ECEFYPOS'},
    'ECEFYVEL':          {'UNITS':'m/s','DEPEND_0':'Epoch','LABLAXIS':'ECEFYVEL'},
    'ECEFZPOS':          {'UNITS':'km', 'DEPEND_0':'Epoch','LABLAXIS':'ECEFZPOS'},
    'ECEFZVEL':          {'UNITS':'m/s','DEPEND_0':'Epoch','LABLAXIS':'ECEFZVEL'},
    'FlightTime':        {'UNITS':'Seconds from Launch','DEPEND_0':'Epoch','LABLAXIS':'FlightTime'},
    'GPSTimemSecofweek': {'UNITS':'Seconds','DEPEND_0':'Epoch','LABLAXIS':'GPSTimemSecofweek'},
    'GPSWeek':           {'UNITS':'Weeks','LABLAXIS':'GPSWeek'},
    'Lat':               {'UNITS':'deg','DEPEND_0':'Epoch','LABLAXIS':'Lat'},
    'Long':              {'UNITS':'deg','DEPEND_0':'Epoch','LABLAXIS':'Long'},
    'Epoch':             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -9223372036854775808,'FORMAT': 'I5','UNITS': 'ns','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'data','MONOTON':'INCREASE','TIME_BASE':'J2000','TIME_SCALE':'Terrestrial Time','REFERENCE_POSITION':'Rotating Earth Geoid','SCALETYP':'linear','LABLAXIS':'Epoch'}
}



# --- --- --- ---
# --- import ---
# --- --- --- ---
import numpy as np
from os import remove,path
from csv import reader
from ACESII_code.data_paths import fliers, ACES_data_folder,ACES_csv_trajectories,TRICE_data_folder,ACES_L0_files,TRICE_L0_files
from ACESII_code.class_var_func import color, L2_TRICE_Quick
from ACESII_code.missionAttributes import ACES_mission_dicts,TRICE_mission_dicts
from copy import deepcopy
setupPYCDF()
from spacepy import pycdf
from spacepy import coordinates as coord
coord.DEFAULTS.set_values(use_irbem=False, itol=5)  # maximum separation, in seconds, for which the coordinate transformations will not be recalculated. To force all transformations to use an exact transform for the time, set ``itol`` to zero.
from spacepy.time import Ticktock #used to determine the time I'm choosing the reference geomagentic field

print(color.BOLD + color.CYAN + 'csv_to_cdf.py' + color.END + color.END)

def csv_to_cdf(inputFile,outputFile,missionDicts,specialMods,wRow):

    # collect the csv data
    with open(inputFile) as csvFile:
        csvAllData = [row for row in reader(csvFile)]

    varNames = csvAllData[wRow]
    csvData = np.array(csvAllData[wRow + 1:],dtype='float64').transpose()

    # --- Create Data Dict ---
    rocketAttrs, missionDict, data_dict_temp = missionDicts
    exampleAttrs = {'DEPEND_0': None,'DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -9223372036854775808,'FORMAT': 'I5','UNITS': None,'VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'data','SCALETYP':'linear',}
    data_dict = {varNames[i] : [ csvData[i] ,  deepcopy(exampleAttrs)  ] for i in range(len(varNames))}


    prgMsg(f'Converting csv Files for {inputFileName}')

    # --- Modify Variable Attributes in data dict to the  ---
    for key, val in specialMods.items():

        if key in data_dict:

            for keyAttr, valAttr in specialMods[key].items():

                if keyAttr in data_dict[key][1]: # modify an existing attribute

                    data_dict[key][1][keyAttr] = valAttr

                else: # append a new attribute
                    data_dict[key][1] = {**data_dict[key][1], **{keyAttr:valAttr}}
        else: # Create a whole new variable
            data_dict = {**data_dict, **{key: [[], val ]} }

    # Create the Epoch Variable
    Launch_time = rocketAttrs.Launch_Times[wRocket]

    Epoch = np.array([ pycdf.lib.tt2000_to_datetime(int(Launch_time + data_dict['FlightTime'][0][i]*(10**9))) for i in range(len(data_dict['FlightTime'][0])) ])
    data_dict['Epoch'][0] = Epoch

    Done(start_time)

    ############################################
    # --- Convert to geomagnetic Coordinates ---
    ############################################
    prgMsg('Calculating geomagnetic coordinates')

    geodetic = np.array([data_dict['Alt'][0],
                         data_dict['Lat'][0],
                         data_dict['Long'][0]]).transpose()

    # Get the times that the mission was launched in ISO datetime. Needed for geomagnetic coordinates
    ISOtime = [data_dict['Epoch'][0][j].isoformat() for j in range(len(data_dict['Epoch'][0]))]

    # --- Convert to geoMAG Coordinates ---

    cvals_GDZ = coord.Coords(geodetic, 'GDZ', 'sph')
    cvals_GDZ.ticks = Ticktock(ISOtime, 'ISO')
    cvals_GDZ_MAG = cvals_GDZ.convert('MAG', 'sph')

    # --- Make all the new variables ---
    geoMagLat = np.array(cvals_GDZ_MAG.lati)
    geoMagLong = np.array(cvals_GDZ_MAG.long)
    geoMagRadi = np.array(cvals_GDZ_MAG.radi)

    data_dict = {**data_dict, **{'geoMag_lat': [geoMagLat, {'LABLAXIS': 'Geomagnetic Latitude',
                                              'DEPEND_0': 'Epoch',
                                              'DEPEND_1': None,
                                              'DEPEND_2': None,
                                              'FILLVAL': -1e30,
                                              'FORMAT': 'E12.2',
                                              'UNITS': 'deg',
                                              'VALIDMIN': geoMagLat.min(),
                                              'VALIDMAX': geoMagLat.max(),
                                              'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
    data_dict = {**data_dict, **{'geoMag_long': [geoMagLong, {'LABLAXIS': 'Geomagnetic Longitude',
                                                            'DEPEND_0': 'Epoch',
                                                            'DEPEND_1': None,
                                                            'DEPEND_2': None,
                                                            'FILLVAL': -1e30,
                                                            'FORMAT': 'E12.2',
                                                            'UNITS': 'deg',
                                                            'VALIDMIN': geoMagLong.min(),
                                                            'VALIDMAX': geoMagLong.max(),
                                                            'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
    data_dict = {**data_dict, **{'geoMag_radi': [geoMagRadi, {'LABLAXIS': 'Geomagnetic Radius',
                                                            'DEPEND_0': 'Epoch',
                                                            'DEPEND_1': None,
                                                            'DEPEND_2': None,
                                                            'FILLVAL': -1e30,
                                                            'FORMAT': 'E12.2',
                                                            'UNITS': 'm',
                                                            'VALIDMIN': geoMagRadi.min(),
                                                            'VALIDMAX': geoMagRadi.max(),
                                                            'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}







    # --- --- --- --- ---
    # --- Output Data ---
    # --- --- --- --- ---

    # --- delete output file if it already exists ---
    if path.exists(outputFile):
        remove(outputFile)
    pycdf.lib.set_backward(False)


    # --- open the output file ---
    with pycdf.CDF(outputFile, '') as cdfFile:
        cdfFile.readonly(False)

        # --- write out global attributes ---
        globalAttrsMod = rocketAttrs.globalAttributes[wRocket]
        ModelData = L2_TRICE_Quick(wRocket)
        inputGlobDic = ModelData.cdfFile.globalattsget()

        for key, val in inputGlobDic.items():
            if key == 'Descriptor':
                globalAttrsMod[key] = 'None'
            if key in globalAttrsMod:
                cdfFile.attrs[key] = globalAttrsMod[key]
            else:
                cdfFile.attrs[key] = val

        # --- WRITE OUT DATA ---
        for varKey, varVal in data_dict.items():
            if 'Epoch' in varKey:  # epoch data
                cdfFile.new(varKey, data=varVal[0], type=33)
            else:  # other data
                cdfFile.new(varKey, data=varVal[0])

            # --- Write out the attributes and variable info ---
            for attrKey, attrVal in data_dict[varKey][1].items():
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
inputFile = inputDirectory+inputFileName
outputFile = outputDirectory + outputFileName
if path.exists(inputFile):
    if wRocket == 0 or wRocket == 1:#ACESII
        missionDicts = ACES_mission_dicts()
    elif wRocket == 2 or wRocket == 3: #TRICEII
        missionDicts = TRICE_mission_dicts()

    csv_to_cdf(inputFile,outputFile,missionDicts,specialMods,wRow)

else:
    print(color.RED + f'{inputFile} does not exist' + color.END)
