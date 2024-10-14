# --- PlotExtra.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Use the CHAOS model to general the Lat/Long/Alt coordinates for "n" B_geo lines and store them in a .cdf file

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False
wRocket = 4
wFile = 0
modifier = ''
inputPath_modifier = 'attitude' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\CHAOSBgeoProjections' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder



tTimes_high = [dt.datetime(2022,11,20,17,24,45,000000),
               dt.datetime(2022,11,20,17,25,00,000000),
               dt.datetime(2022,11,20,17,25,15,000000) ]

tTimes_low = [
               dt.datetime(2022,11,20,17,24,45,000000),
               dt.datetime(2022,11,20,17,25,00,000000),
               dt.datetime(2022,11,20,17,25,15,000000)]
tTimes = [tTimes_high,tTimes_low]
t_Params = [-1000 + 200*i for i in range(10)]
# ---------------------------
outputData = True
# ---------------------------

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.class_var_func import CHAOS,ENUtoECEF,Re,ECEF_to_Geodedic
from spacepy import coordinates as coord
from spacepy.time import Ticktock #used to determine the time I'm choosing the reference geomagentic field

def PlotExtra_conjugactconcept(wRocket, rocketFolderPath, justPrintFileNames):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wRocket-4]
    inputFiles, input_names, input_names_searchable = getInputFiles(rocketFolderPath=rocketFolderPath,wRocket=wRocket,inputPath_modifier=inputPath_modifier)
    fileoutName = f'ACESII_{rocketID}_CHAOS_Bgeo_LatLongAlt.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return


    # --- --- --- --- --- -
    # --- LOAD THE DATA ---
    # --- --- --- --- --- -

    # --- get the data from the file ---
    prgMsg(f'Loading data from {inputPath_modifier} Files')
    data_dict_attitude = loadDictFromFile(inputFiles[wFile])
    Done(start_time)

    exampleAttrs = {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': rocketAttrs.epoch_fillVal,
                    'FORMAT': 'E12.2', 'UNITS': None, 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data', 'SCALETYP': 'linear',
                    'LABLAXIS': None}

    targetTimes = tTimes[wRocket-4]

    # --- --- --- --- --- --- --
    # --- GET THE CHAOS MODE ---
    # --- --- --- --- --- --- --
    prgMsg('Calculating Bgeo Lines')

    pointsOfInterest_Geographic = [[] for tme in targetTimes]

    for tIDX, tme in enumerate(targetTimes):

        # STEP 1: Get B_geo in ENU coordinates
        idx = np.abs(data_dict_attitude['Epoch'][0] - tme).argmin()
        rktLat = data_dict_attitude['Lat'][0][idx]
        rktLong = data_dict_attitude['Long'][0][idx]
        rktLAlt = data_dict_attitude['Alt'][0][idx]

        B_model = CHAOS(lat=np.array([rktLat]),
                        long=np.array([rktLong]),
                        alt=np.array([rktLAlt/1000]),
                        times=np.array([tme]))  # CHAOS in ENU coordinates

        b = np.array(B_model[0]) / np.linalg.norm(B_model[0])
        print(b)



        # STEP 2: Get B_geo in ECEF coordinates
        b_ECEF = np.matmul(b, ENUtoECEF(rktLat, rktLong))
        print(b_ECEF)


        # STEP 3: Get Position of Rocket in ECEF coordinates
        coord.DEFAULTS.set_values(use_irbem=False, itol=5)
        rkt0_geodetic = np.array([rktLAlt/1000, rktLat, rktLong])
        ISOtime = tme.isoformat()
        cvals_GDZ = coord.Coords(rkt0_geodetic, 'GDZ', 'sph')
        cvals_GDZ.ticks = Ticktock(ISOtime, 'ISO')
        cvals_GEO = cvals_GDZ.convert('GEO', 'car')
        rkt0_ECEF = np.array([cvals_GEO.x[0],cvals_GEO.y[0],cvals_GEO.z[0]]) # In units of Earth Radii measured from Earth's Center

        # STEP 4: Get B_geo points of interest in ECEF coordinates
        # note: With b being a unit vector in ECEF space --> each movement in b_ECEF is 1 kilometer! But this is 1kilometer in the direction of hat{b}
        pointsOfInterest_ECEF = []
        for t in t_Params:
            pointsOfInterest_ECEF.append(t*b_ECEF/Re + rkt0_ECEF)


        # STEP 5: Convert points of interest to Geodedic coordinates and store them

        for point_set in pointsOfInterest_ECEF:
            convertedData = ECEF_to_Geodedic(point_set[0]*Re*1000,point_set[1]*Re*1000,point_set[2]*Re*1000)
            pointsOfInterest_Geographic[tIDX].append(
                [(2*rktLat-(180/np.pi)*convertedData['lat']), (2*rktLong - (180/np.pi)*convertedData['lon']),convertedData['height']]
            )
            # print([(180/np.pi)*convertedData['lat'],(180/np.pi)*convertedData['lon'],convertedData['height']/1000])


    Done(start_time)




    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        prgMsg('Creating output file')

        data_dict_output={}
        for tme in range(len(pointsOfInterest_Geographic)):
            data_dict_output = {**data_dict_output, **{f'Line {tme}': [np.array(pointsOfInterest_Geographic[tme]),deepcopy(exampleAttrs)]}}


        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wRocket-4]}\\{fileoutName}'

        outputCDFdata(outputPath, data_dict_output)

        Done(start_time)





# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

rocketFolderPath = ACES_data_folder


if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    PlotExtra_conjugactconcept(wRocket, rocketFolderPath, justPrintFileNames)

