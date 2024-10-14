# --- csv_to_cdf_attitude.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Turn the .cdf files of the ACESII attitude solution data into cdf files


# --- --- --- --- ---
from ACESII_code.myImports import *

start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---


wRocket = 0
# 0 --> ACESII High Flyer
# 1 --> ACESII Low Flyer


Directories = [r'C:\Data\ACESII\attitude\high\\',r'C:\Data\ACESII\attitude\low\\']
inputFilenames = [r'36359_AttitudeSolution.csv',r'36364_AttitudeSolution.csv']
outputFilenames = [r'ACESII_36359_Attitude_Solution.cdf',r'ACESII_36364_Attitude_Solution.cdf']
inputDirectory = Directories[wRocket]
inputFileName = inputFilenames[wRocket]
outputDirectory = inputDirectory
outputFileName = outputFilenames[wRocket]

wRow_names = 5 -1 # Which row contains the variable names
wRow_units = 6 - 1 # Which row contains the variable units
wRow_startData = 8 - 1 # Which row contains the first row of actual data
# ---------------------------
# --- special modifications to the variables already in the csv file ---
#  in a dictionary format [ {string name of var: {string name of attributeKey1: attributeval1, attributeKey2: attributeval2  }},{...},...  ]
# Example: specialMods = [time: {'Units':'ns','FORMAT':'I5'} , Alt: {'Units':'km','FORMAT':'I5'} ]
# Defaults to attributes = None if specialMods == [], except for VALIDMIN/VALIDMAX
specialMods = {
    'Alt':               {'DEPEND_0':'Epoch','LABLAXIS':'Alt'},
    'Epoch':             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -9223372036854775808,'FORMAT': 'I5','UNITS': 'ns','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'data','MONOTON':'INCREASE','TIME_BASE':'J2000','TIME_SCALE':'Terrestrial Time','REFERENCE_POSITION':'Rotating Earth Geoid','SCALETYP':'linear','LABLAXIS':'Epoch'}
}
# ---------------------------
outputData = True
# ---------------------------


# --- --- --- ---
# --- import ---
# --- --- --- ---
from os import path
from csv import reader
from ACESII_code.class_var_func import color, L2_TRICE_Quick
from spacepy import coordinates as coord
coord.DEFAULTS.set_values(use_irbem=False, itol=5)  # maximum separation, in seconds, for which the coordinate transformations will not be recalculated. To force all transformations to use an exact transform for the time, set ``itol`` to zero.
from spacepy.time import Ticktock #used to determine the time I'm choosing the reference geomagentic field

print(color.BOLD + color.CYAN + 'csv_to_cdf_attitude.py' + color.END + color.END)

def csv_to_cdf_attitude(inputFile,outputFile,missionDicts,specialMods,wRow_names,wRow_units,wRow_startData):

    # collect the entire csv data
    with open(inputFile) as csvFile:
        csvAllData = [row for row in reader(csvFile)]

    varNames = np.array(csvAllData[wRow_names])
    varNames = np.array([nam.replace(' ','') for nam in varNames])


    varUnits = np.array(csvAllData[wRow_units])
    csvData = np.array(csvAllData[wRow_startData:], dtype='float64').transpose()

    # --- Create Data Dict ---
    data_dict = {}
    rocketAttrs, missionDict, data_dict_temp = missionDicts
    exampleAttrs = {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -9223372036854775808,
                    'FORMAT': 'I5', 'UNITS': None, 'LABLAXIS':None, 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data',
                    'SCALETYP': 'linear'}

    for index, var in enumerate(varNames):
        varAttrs = deepcopy(exampleAttrs)
        unit = str(varUnits[index]).replace('(','').replace(')','').replace(' ','')
        if unit == 'm':
            unit = 'meters'

        varAttrs['LABLAXIS'] = var
        varAttrs['UNITS'] = unit
        data_dict = {**data_dict,**{varNames[index] : [csvData[index],  varAttrs ]}}

    prgMsg(f'Converting csv Files for {inputFileName}')

    # --- Modify Variable Attributes in data dict  ---
    for key, val in specialMods.items():

        if key in data_dict:

            for keyAttr, valAttr in specialMods[key].items():

                if keyAttr in data_dict[key][1]: # modify an existing attribute
                    data_dict[key][1][keyAttr] = valAttr

                else: # append a new attribute
                    data_dict[key][1] = {**data_dict[key][1], **{keyAttr:valAttr}}

        else: # Create a whole new variable
            data_dict = {**data_dict, **{key: [[], deepcopy(val) ]} }

    # rename latgd --> lat
    data_dict['Lat'] = data_dict.pop('Latgd')

    # Create the Epoch Variable
    Launch_time = rocketAttrs.Launch_Times[wRocket]
    Epoch = np.array([ pycdf.lib.tt2000_to_datetime(int(Launch_time + data_dict['Time'][0][i]*(10**9))) for i in range(len(data_dict['Time'][0])) ])
    data_dict['Epoch'][0] = Epoch
    Done(start_time)

    # --- Create the Geomagnetic Coordinates ---
    prgMsg('Converting Coordinates')

    # Calculate
    geodeticPos = np.array([data_dict['Alt'][0]/1000, data_dict['Lat'][0], data_dict['Long'][0]]).transpose()
    ISOtime = [data_dict['Epoch'][0][i].isoformat() for i in range(len(data_dict['Epoch'][0]))]
    cvals_GDZ = coord.Coords(geodeticPos, 'GDZ', 'sph')
    cvals_GDZ.ticks = Ticktock(ISOtime, 'ISO')
    cvals_GDZ_MAG = cvals_GDZ.convert('MAG', 'sph')
    geomagAlt = cvals_GDZ_MAG.radi
    geomagLat = cvals_GDZ_MAG.lati
    geomagLong = cvals_GDZ_MAG.long

    # store data
    data_dict = {**data_dict, **{'Lat_geom':[np.array(geomagLat),deepcopy(data_dict['Lat'][1])]}}
    data_dict['Lat_geom'][1]['LABLAXIS'] = 'geomagnetic Lat'
    data_dict = {**data_dict, **{'Long_geom':[np.array(geomagLong),deepcopy(data_dict['Long'][1])]}}
    data_dict['Long_geom'][1]['LABLAXIS'] = 'geomagnetic Long'
    data_dict = {**data_dict, **{'Alt_geom':[np.array(geomagAlt),deepcopy(data_dict['Alt'][1])]}}
    data_dict['Alt_geom'][1]['LABLAXIS'] = 'geomagnetic Alt'
    data_dict['Alt_geom'][1]['UNITS'] = 'Re'

    Done(start_time)

    # --- Calculate Invariant Lattitude ---
    invariantLat = np.array([np.degrees(np.arccos(np.cos(np.radians(data_dict['Lat_geom'][0][i])) / np.sqrt(data_dict['Alt_geom'][0][i]))) for i in range(len(data_dict['Epoch'][0]))])
    data_dict = {**data_dict,**{'invarLat':[invariantLat,deepcopy(data_dict['Lat_geom'][1])]}}
    data_dict['invarLat'][1]['LABLAXIS'] = 'Invariant Lat'

    if outputData:
        prgMsg('Creating output file')
        globalAttrsMod = rocketAttrs.globalAttributes[wRocket]
        outputCDFdata(outputFile, data_dict,globalAttrsMod=globalAttrsMod,instrNam= 'Attitude')
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

    csv_to_cdf_attitude(inputFile,outputFile,missionDicts,specialMods,wRow_names,wRow_units,wRow_startData)

else:
    print(color.RED + f'{inputFile} does not exist' + color.END)
