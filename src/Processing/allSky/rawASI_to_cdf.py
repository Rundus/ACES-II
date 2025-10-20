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
from src.my_imports import *

start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# --- Select the Site ---
# 0 -> Skibotn
wSite = 0
architecture = 65535 # type of bit-architecture used for the allSky images. Nomially 16-bit
elevlimits = [20, 20] # value (in deg) of the cutoff elevation angle for the 5577A and 6300A line

wLengths = ['5577', '6300']
site_nam = ['skibotn']

outputData = True


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.io import readsav
import matplotlib.pyplot as plt

def rawASI_to_cdf(wSite, rocketFolderPath):

    # --- load attributes for ACESII traj data ---
    folder_path = f'{rocketFolderPath}\\all_sky\\{site_nam[wSite]}\\'

    # --- get the input files ---
    cal_files = [readsav(glob(folder_path+'\\5577\\*.dat*')[0]), readsav(glob(folder_path + '\\6300\\*.dat*')[0])]
    photo_files = [glob(folder_path + '\\5577\\*.png*'), glob(folder_path + '\\6300\\*.png*')]

    ############################################
    # --- COLLECT IMAGE FILES AND TIMESTAMPS ---
    ############################################

    # get the image time series and the data itself into single variables
    Epoch_AllSky = [[], []]
    imageData = [[], []]
    stl.prgMsg('Collecting Image Data')

    for i in range(len(photo_files)):
        for imageStr in photo_files[i]:

            # get the timestamp
            strTimeStamp = imageStr.replace(f'{folder_path}\\','').replace(f'{wLengths[i]}\\','').replace('skn4_','').replace(f'_{wLengths[i]}_cal.png','')
            year = int(strTimeStamp[0:4])
            month = int(strTimeStamp[4:6])
            day = int(strTimeStamp[6:8])
            hour = int(strTimeStamp[9:11])
            minute = int(strTimeStamp[11:13])
            second = int(strTimeStamp[13:15])
            dtTimeStamp = dt.datetime(year, month, day, hour, minute, second)
            Epoch_AllSky[i].append(dtTimeStamp)

            # get the grayscale data
            imageData[i].append(plt.imread(imageStr))

    stl.Done(start_time)

    ######################
    # --- prepare data ---
    ######################
    allGlats = [np.array(deepcopy(cal_files[0]['glats']))[::-1], np.array(deepcopy(cal_files[1]['glats']))[::-1]]
    allGLongs = [np.array(deepcopy(cal_files[0]['glons']))[::-1], np.array(deepcopy(cal_files[1]['glons']))[::-1]]
    allElevs = [np.array(deepcopy(cal_files[0]['elevs'])), np.array(deepcopy(cal_files[1]['elevs']))]
    allImages = imageData



    # stl.prgMsg('Collecting and processing/cleaning Image data')
    #
    # # -- collect and process all aurora image data ---
    # # remove nan values from data and replace with garbage AND convert all Images into Rayleighs
    # # description: the image data is a 16-bit counts value normalized by 65535.
    # # Invert this and use the calibration factor of 1 R/count given in the cal file
    #
    # for i in range(2):  # wavelength
    #
    #     # --- correct the calibration data ---
    #     for j in range(len(allGlats[i])):  # array rows
    #         for k in range(len(allGlats[i][j])):  # row values
    #             if math.isnan(allGlats[i][j][k]):
    #                 allGlats[i][j][k] = 70
    #                 for a in range(len(allImages[i])):  # correct this j,k point in all auroral images
    #                     allImages[i][a][j][k] = np.nan
    #
    #             if math.isnan(allGLongs[i][j][k]):
    #                 allGLongs[i][j][k] = 20
    #                 for a in range(len(allImages[i])):  # correct this j,k point in all auroral images
    #                     allImages[i][a][j][k] = np.nan
    #
    #             if allElevs[i][j][k] <= elevlimits[i] or math.isnan(allElevs[i][j][k]):
    #                 for a in range(len(allImages[i])):  # correct this j,k point in all auroral images
    #                     allImages[i][a][j][k] = np.nan
    #
    #     # --- convert images to rayleighs ---
    #     for a in range(len(allImages[i])): #number of images
    #         for j in range(len(allImages[i][a])): # arrays in particular image
    #             for k in range(len(allImages[i][a][j])): # for each value in each array of image
    #                 if not math.isnan(allImages[i][a][j][k]):
    #                     allImages[i][a][j][k] = int(allImages[i][a][j][k]*architecture)
    #
    # stl.Done(start_time)
    #
    # # --- determine the timestamps for each photo and rocket data ---
    # Epoch_AllSky_tt2000 = [[pycdf.lib.datetime_to_tt2000(time) for time in Epoch_AllSky[0]],
    #                        [pycdf.lib.datetime_to_tt2000(time) for time in Epoch_AllSky[1]]]
    #
    # # for each image timestamp, find the index of the closest Epoch_esa in the High flyer's epoch
    # imageIndiciesToRocket = [[], []]
    # for i in range(len(Epoch_AllSky_tt2000)):
    #     for stamp in Epoch_AllSky_tt2000[i]:
    #         imageIndiciesToRocket[i].append(np.abs(data_dicts_traj[0]['Epoch_esa'][0] - stamp).argmin())
    #
    # imageIndiciesToRocket = np.array(imageIndiciesToRocket)
    #
    # # fill in imageIndiciesToRocket with indicies to match the size of High Flyer Epoch data
    # imgIndicies = [[], []]
    #
    # for i in range(2):
    #     for j in range(len(imageIndiciesToRocket[i])):  # loop through the various images
    #         if j == 0:
    #
    #             for k in range(imageIndiciesToRocket[i][j + 1]):
    #                 imgIndicies[i].append(j)
    #
    #         elif j == len(imageIndiciesToRocket[i]) - 1:  # if you're at the last image index
    #             for k in range(len(EpochRocket[0]) - imageIndiciesToRocket[i][j]):
    #                 imgIndicies[i].append(j)
    #
    #         else:
    #             for k in range((imageIndiciesToRocket[i][j + 1] - imageIndiciesToRocket[i][j])):
    #                 imgIndicies[i].append(j)
    #
    # imgIndicies = np.array(imgIndicies)
    # EpochMovie = np.array([pycdf.lib.tt2000_to_datetime(EpochRocket[0][i]) for i in range(len(EpochRocket[0]))])
    #
    # ###############################
    # # --- EXTEND LOW FLYER DATA ---
    # ###############################
    #
    # # --- extend Low Flyer Rocket data to be the same length as High flyer in the beginning and end ---
    # highFlyerSize, lowFlyerSize = len(geoAlt[0]), len(geoAlt[1])
    #
    # # --- Append start Values to  ---
    # no_of_points_start = int((EpochRocket[1][0] - EpochRocket[0][0]) / (rocketAttrs.MinorFrameTime))
    # newAlt = [geoAlt[1][0] for i in range(no_of_points_start)]
    # newLat = [geoLat[1][0] for i in range(no_of_points_start)]
    # newLong = [geoLong[1][0] for i in range(no_of_points_start)]
    #
    # for i in range(len(geoAlt[1])):
    #     newAlt.append(geoAlt[1][i])
    #     newLat.append(geoLat[1][i])
    #     newLong.append(geoLong[1][i])
    #
    # # --- Append the ending values ---
    # remainingIndicies = highFlyerSize - (lowFlyerSize + no_of_points_start)
    #
    # for i in range(remainingIndicies):
    #     newAlt.append(geoAlt[1][-1])
    #     newLat.append(geoLat[1][-1])
    #     newLong.append(geoLong[1][-1])
    #
    # geoAlt[1], geoLat[1], geoLong[1] = np.array(newAlt), np.array(newLat), np.array(newLong)








# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

rocketFolderPath = DataPaths.ACES_data_folder
rawASI_to_cdf(wSite,rocketFolderPath)
