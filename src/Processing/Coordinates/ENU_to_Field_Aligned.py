# --- ENU_to_Field_Aligned.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: input a integration_tad_files or EFI file and use the respective attitude solution to convert
# to field aligned coordinates. Additionally, use the data from the attitude solution to
# determine the geomagnetic coordinates of the deltaB/deltaE data on the mag epoch
# and store them in the output file



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

from ACESII_code.myImports import *

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
wRocket = 5

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [7]

modifier = ''
inputPath_modifier_attitude = 'attitude' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier = r'\l2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = r'\l2' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.class_var_func import CHAOS
from spacepy import coordinates as coord
coord.DEFAULTS.set_values(use_irbem=False, itol=5)  # maximum separation, in seconds, for which the coordinate transformations will not be recalculated. To force all transformations to use an exact transform for the time, set ``itol`` to zero.
from spacepy.time import Ticktock #used to determine the time I'm choosing the reference geomagentic field

def ENU_to_Field_Aligned(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf')
    inputFiles_attitude = glob(f'{rocketFolderPath}{inputPath_modifier_attitude}\{fliers[wflyer]}{modifier}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]

    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    fileoutName = input_names[wFile].replace('ENU','Field_Aligned')


    if 'RingCore' in fileoutName:
        wInstr = 'RingCore'
    elif 'E_Field' in fileoutName:
        wInstr = 'E_Field'


    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return


    # --- get the data from the E-Field or B-Field file ---
    prgMsg(f'Loading data from {inputPath_modifier} Files')
    data_dict = loadDictFromFile(inputFiles[wFile])
    Done(start_time)

    # --- get the data from the attitude file ---
    prgMsg(f'Loading data from {inputPath_modifier_attitude} Files')
    data_dict_attitude = loadDictFromFile(inputFiles_attitude[0],targetVar=[[data_dict['Epoch'][0][0],data_dict['Epoch'][0][-1]],'Epoch'])
    Done(start_time)

    ###############################################################
    # --- interpolate attitude data up to B-Field/E-Field epoch ---
    ###############################################################
    prgMsg('Interpolating Attitude Data')
    data_dict_attitudeInterp = InterpolateDataDict(InputDataDict=data_dict_attitude,
                                                   InputEpochArray=data_dict_attitude['Epoch'][0],
                                                   wKeys=['Alt','Lat','Long','Epoch'],
                                                   targetEpochArray=data_dict['Epoch'][0])

    # convert altitude to km
    data_dict_attitudeInterp['Alt'][0] = data_dict_attitudeInterp['Alt'][0] / 1000
    Done(start_time)

    # #############################################################
    # --- Convert ENU coordinates to Field Aligned coordinates ---
    # #############################################################
    prgMsg(f'Loading CHAOS model for {len(data_dict["Epoch"][0])} points')

    # get the key components
    compNames = [key for key, val in data_dict.items() if key.lower() not in ['de_mag', 'db_mag', 'epoch','e_mag','ilat','ilong','bmag','b_mag','alt']]

    # Get the Data
    B_rkt_ENU = np.array([[data_dict[compNames[0]][0][i], data_dict[compNames[1]][0][i], data_dict[compNames[2]][0][i]] for i in range(len(data_dict['Epoch'][0]))])
    B_model = CHAOS(lat=data_dict_attitudeInterp['Lat'][0],
                    long=data_dict_attitudeInterp['Long'][0],
                    alt=data_dict_attitudeInterp['Alt'][0],
                    times=data_dict['Epoch'][0])  # CHAOS in ENU coordinates

    # --- Convert B-Data to GEO (ECEF) XYZ coordinates ---
    from ACESII_code.class_var_func import ENUtoECEF
    ENUtoGEOmatrix = np.array([ENUtoECEF(Lat=data_dict_attitudeInterp['Lat'][0][i],
                                         Long=data_dict_attitudeInterp['Long'][0][i]) for i in range(len(data_dict['Epoch'][0]))])

    B_rkt_GEO = np.array([np.matmul(ENUtoGEOmatrix[i], B_rkt_ENU[i]) for i in range(len(data_dict['Epoch'][0]))])
    B_CHAOS_GEO = np.array([np.matmul(ENUtoGEOmatrix[i], B_model[i]) for i in range(len(data_dict['Epoch'][0]))])

    # --- determine the Payload's Position Vector in GEO (ECEF) coordinate XYZ ---
    R_REF = 6371.2  # earth Radius in km
    Radius = data_dict_attitudeInterp['Alt'][0] + R_REF
    coLatRad = [np.radians(90 - lat) for lat in data_dict_attitudeInterp['Lat'][0]]
    LongRad = [np.radians(long) for long in data_dict_attitudeInterp['Long'][0]]
    Rsc = np.array([
        [Radius[i] * np.sin(coLatRad[i]) * np.cos(LongRad[i]),
         Radius[i] * np.sin(coLatRad[i]) * np.sin(LongRad[i]),
         Radius[i] * np.cos(coLatRad[i])] for i in range(len(data_dict['Epoch'][0]))])

    Done(start_time)

    prgMsg('Converting to Field Aligned Coordinates')
    # --- calculate Field Aligned unit vectors over the duration of the flight ---

    # pHat comes from the CHAOS model direction of B in GEO
    pHat = np.array([B_CHAOS_GEO[i] / np.linalg.norm(B_CHAOS_GEO[i]) for i in range(len(data_dict['Epoch'][0]))])

    # e-hat comes from the cross of pHat and the Rocket's radius vector (in geomagnetic coordinates)
    eHat = np.array([np.cross(pHat[i], Rsc[i]) / np.linalg.norm(np.cross(pHat[i], Rsc[i])) for i in range(len(data_dict['Epoch'][0]))])

    # rHat comes from the cross of eHat and pHat
    rHat = np.array([np.cross(eHat[i], pHat[i]) for i in range(len(data_dict['Epoch'][0]))])

    # form the transformation matrix FROM GEO TO FIELD ALIGNED
    DCM_FA = np.array([[eHat[i], pHat[i], rHat[i]] for i in range(len(data_dict['Epoch'][0]))])

    # --- Transform B_ENU to B_EPR ---
    B_rkt_EPR = np.array([np.matmul(DCM_FA[i], B_rkt_GEO[i]) for i in range(len(data_dict['Epoch'][0]))])

    data_for_output = deepcopy(B_rkt_EPR)
    Done(start_time)

    ############################################
    # --- CALCULATE GEOMAGNETIC COORDINATES ---
    ############################################
    prgMsg('Calculating Geomagnetic Coordinates')

    # Trajectory Data
    geodeticPos = np.array([data_dict_attitudeInterp['Alt'][0], data_dict_attitudeInterp['Lat'][0], data_dict_attitudeInterp['Long'][0]]).transpose()
    ISOtime = np.array([data_dict['Epoch'][0][i].isoformat() for i in range(len(data_dict['Epoch'][0]))])
    cvals_GDZ = coord.Coords(geodeticPos, 'GDZ', 'sph')
    cvals_GDZ.ticks = Ticktock(ISOtime, 'ISO')
    cvals_GDZ_MAG = cvals_GDZ.convert('MAG', 'sph')
    geomagAlt = np.array(cvals_GDZ_MAG.radi)
    geomagLat = np.array(cvals_GDZ_MAG.lati)
    geomagLong = np.array(cvals_GDZ_MAG.long)

    Done(start_time)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        prgMsg('Creating FieldAligned output file')

        # create the output data_dict # NOTE: We inheret variables like Bmag from the original data and they aren't modified here. Only the components
        data_dict_output = deepcopy(data_dict)
        newComps = ['B_e', 'B_p', 'B_r'] if wInstr in ['RingCore'] else ['E_e', 'E_p', 'E_r']
        data_for_output_FA = np.array([[data_for_output[i][0], data_for_output[i][1], data_for_output[i][2]] for i in range(len(data_for_output))])

        # --- Insert Magnetic or Electric Components ---
        # get the attributes of the old components and replace them

        for i, key in enumerate(compNames):
            newAttrs = deepcopy(data_dict_output[key][1])
            newAttrs['LABLAXIS'] = newComps[i]
            newAttrs['DEPEND_0'] = 'Epoch'

            # remove the old key
            del data_dict_output[key]

            # append the new key
            data_dict_output = {**data_dict_output, **{newComps[i]: [data_for_output_FA[:, i], newAttrs]}}


        # --- Insert the geographic/geomagnetic coordinates ---

        # geograph altitude
        varAttrs = deepcopy(data_dict_attitude['Alt'][1])
        varAttrs['LABLAXIS'] = 'Alt'
        varAttrs['UNITS'] = 'km'
        data_dict_output = {**data_dict_output, **{'Alt': [np.array(data_dict_attitudeInterp['Alt'][0]), varAttrs]}}

        # geograph lat
        varAttrs = deepcopy(data_dict_attitude['Alt'][1])
        varAttrs['LABLAXIS'] = 'Latitude'
        varAttrs['UNITS'] = 'deg'
        data_dict_output = {**data_dict_output, **{'Lat': [np.array(data_dict_attitudeInterp['Lat'][0]), varAttrs]}}

        # geograph long
        varAttrs = deepcopy(data_dict_attitude['Alt'][1])
        varAttrs['LABLAXIS'] = 'Longitude'
        varAttrs['UNITS'] = 'deg'
        data_dict_output = {**data_dict_output, **{'Long': [np.array(data_dict_attitudeInterp['Long'][0]), varAttrs]}}

        # geomagnetic altitude
        varAttrs = deepcopy(data_dict_attitude['Alt'][1])
        varAttrs['LABLAXIS'] = 'geomagnetic Alt'
        varAttrs['UNITS'] = 'km'
        data_dict_output = {**data_dict_output, **{'Alt_geom': [np.array(geomagAlt), varAttrs]}}

        # geomagnetic latitude
        varAttrs = deepcopy(data_dict_attitude['Alt'][1])
        varAttrs['LABLAXIS'] = 'geomagnetic Lat'
        varAttrs['UNITS'] = 'deg'
        data_dict_output = {**data_dict_output, **{'Lat_geom': [np.array(geomagLat), varAttrs]}}

        # geomagnetic Longitude
        varAttrs = deepcopy(data_dict_attitude['Alt'][1])
        varAttrs['LABLAXIS'] = 'geomagnetic Long'
        varAttrs['UNITS'] = 'deg'
        data_dict_output = {**data_dict_output, **{'Long_geom': [np.array(geomagLong), varAttrs]}}

        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

        outputCDFdata(outputPath, data_dict_output,instrNam='E_Field')

        Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        ENU_to_Field_Aligned(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            ENU_to_Field_Aligned(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            ENU_to_Field_Aligned(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)