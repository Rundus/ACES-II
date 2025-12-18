# --- rawASI_to_cdf.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Loads in the AllSky data, uses the calibration file to determine position
# finally loads in traj data to determine rocket trajectory


# assumes all light comes from these altitudes:
# 557.7nm: 150km
# 630nm: 250km


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
from src.ACESII.my_imports import *

start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# --- Select the Site ---
# 0 -> Skibotn
wSite = 0
site_nam = ['skibotn']

# --- Select the color ---
# 0 -> green 5570
# 1 -> red 6300
wColor = 1
wLengths = ['5577', '6300']
color_alt = [150, 250] # in kilometers

# additional toggles
architecture = 65535 # type of bit-architecture used for the allSky images. Nomially 16-bit
elevlimits = [20, 20] # value (in deg) of the cutoff elevation angle for the 5577A and 6300A line
K = 1 # Calibration value. ==1 If unlisted at http://tid.uio.no/plasma/aurora/tech.html

# writing out data
outputData = True



# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.io import readsav
import matplotlib.pyplot as plt
import math
import h5py
from itertools import product

def rawASI_to_cdf(wSite, wColor, rocketFolderPath):

    # --- load attributes for ACESII traj data ---
    folder_path = f'{rocketFolderPath}\\all_sky\\{site_nam[wSite]}\\'

    # --- get the input files ---
    cal_files = readsav(glob(folder_path+f'\\{wLengths[wColor]}\\*.dat*')[0])
    photo_files = glob(folder_path + f'\\{wLengths[wColor]}\\*.png*')

    # --- prepare the output ---
    data_dict_output = {}
    data_dict_globalAttrs = {}

    ############################################
    # --- COLLECT IMAGE FILES AND TIMESTAMPS ---
    ############################################

    # get the image time series and the data itself into single variables
    Epoch_AllSky = []
    imageData = []
    stl.prgMsg('Collecting Image Data')

    for imageStr in photo_files:

        # get the timestamp
        strTimeStamp = imageStr.replace(f'{folder_path}\\','').replace(f'{wLengths[wColor]}\\','').replace('skn4_','').replace(f'_{wLengths[wColor]}_cal.png','')
        year = int(strTimeStamp[0:4])
        month = int(strTimeStamp[4:6])
        day = int(strTimeStamp[6:8])
        hour = int(strTimeStamp[9:11])
        minute = int(strTimeStamp[11:13])
        second = int(strTimeStamp[13:15])
        dtTimeStamp = dt.datetime(year, month, day, hour, minute, second)
        Epoch_AllSky.append(dtTimeStamp)

        # get the grayscale data
        imageData.append(plt.imread(imageStr))

    # --- store data ---
    for key in cal_files.keys():
        if key in ['glats','glons','mlats','mlons','alts','gazms','mazms','elevs']:
            data_dict_output = {**data_dict_output,
                                **{
                                    f'{key}':[np.array(deepcopy(cal_files[f'{key}']))[::-1],{}]
                                }}
        else:
            data_dict_globalAttrs = {**data_dict_globalAttrs,
                                **{
                                    f'{key}': cal_files[f'{key}']
                                }}
    allImages = imageData
    stl.Done(start_time)

    # -- collect and process all aurora image data ---
    # remove nan values from data and replace with garbage AND convert all Images into Rayleighs
    # description: the image data is a 16-bit counts value normalized by 65535.
    # Invert this and use the calibration factor of 1 R/count given in the cal file
    stl.prgMsg('Applying Calibrations')

    # # --- correct the calibration data ---
    # for j in range(len(allGlats)):  # array rows
    #     for k in range(len(allGlats[j])):  # row values
    #         if math.isnan(allGlats[j][k]):
    #             allGlats[j][k] = 70
    #             for i in range(len(allImages)):  # correct this j,k point in all auroral images
    #                 allImages[i][j][k] = np.nan
    #
    #         if math.isnan(allGLongs[j][k]):
    #             allGLongs[j][k] = 20
    #             for i in range(len(allImages)):  # correct this j,k point in all auroral images
    #                 allImages[i][j][k] = np.nan
    #
    #         if allElevs[j][k] <= elevlimits or math.isnan(allElevs[j][k]):
    #             for i in range(len(allImages)):  # correct this j,k point in all auroral images
    #                 allImages[i][j][k] = np.nan

    # --- convert images to rayleighs ---
    ranges = [range(len(allImages)),range(len(allImages[0])), range(len(allImages[0][0]))]
    for i, j, k in product(*ranges):
        if not math.isnan(allImages[i][j][k]):
            allImages[i][j][k] = int(allImages[i][j][k]*architecture*K)
    stl.Done(start_time)

    ############################
    # --- WRITE OUT THE DATA ---
    ############################
    if outputData:

        data_dict_output = {
            **data_dict_output,
            **{
                'Epoch': [np.array(Epoch_AllSky), {'DEPEND_0':'Epoch'}],
                'images': [np.array(allImages), {'DEPEND_0':'Epoch', 'UNITS':'Rayleigh'}]
            }
        }

        # adjust some labels
        data_dict_output['glons'][1] = {'UNITS': 'Degrees', 'CATDESC': 'geographic longitude'}
        data_dict_output['glats'][1] = {'UNITS': 'Degrees', 'CATDESC': 'geographic latitude'}
        data_dict_output['alts'][1] = {'UNITS': 'km'}
        data_dict_output['elevs'][1] = {'UNITS': 'Degrees', 'CATDESC': 'pixel elevation'}

        # Create the .CDF file
        output_path = f'{folder_path}\\{wLengths[wColor]}\\ACESII_AllSky_skibotn_{wLengths[wColor]}'
        stl.outputDataDict(outputPath=output_path + '.cdf',
                          data_dict=data_dict_output,
                           globalAttrsMod=data_dict_globalAttrs)

        # create the .hdf5 file
        with h5py.File(rf'{output_path}.hdf5', 'w') as f:
            for key, val in data_dict_output.items():
                if key in ['Epoch']:
                    f[f'{key}_tt2000'] = np.array([pycdf.lib.datetime_to_tt2000(v) for v in  val[0]])
                else:
                    f[f'{key}'] = val[0]

                for subkey, subval in val[1].items():
                    if key in ['Epoch']:
                        f[f'{key}_tt2000'].attrs[f'{subkey}'] = subval
                    else:
                        f[f'{key}'].attrs[f'{subkey}'] = subval








# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

rocketFolderPath = DataPaths.ACES_data_folder
rawASI_to_cdf(wSite, wColor, rocketFolderPath)
