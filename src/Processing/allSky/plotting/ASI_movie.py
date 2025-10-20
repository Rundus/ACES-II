# --- AllSkyTrajecMovie.py ---
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

# Just print the names of files
justPrintSiteNames = False

# --- Select the Site ---
# 0 -> Skibotn
wSite = 0

modifier = ''
inputPath_modifier_AllSky = 'all_sky' # e.g. 'L1' or 'L1'. It's the name of the broader input folder inside data\ACESII
inputPath_modifier_traj = 'trajectories'
inputPath_modifier_magGeo = 'mag'
outputPath_modifier = '\\Movies_and_Media\\movies' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder inside data\ACESII\ACESII_matlab

architecture = 65535 # type of bit-architecture used for the allSky images. Nomially 16-bit
elevlimits = [20, 20] # value (in deg) of the cutoff elevation angle for the 5577A and 6300A line

# --- AllSky MOVIE ---
createAllSkyMovie = True
fps = 20 # fps of the video
projectionAltitude = [150, 250] # in km. Format: [green, red]

# plot specific locations
plotSpecificLocations = False
specificLocations = [i*10 for i in range(int(12371/10))]

# plot frame skips
doFrameSkips = True
frame_skips = 5 # how many Epoch_esa frames to skip

# Kenton wanted just the low flyer for a plot, so here's a toggle to do that
useOnlyLowFlyer = False

# outputFileName = 'ACESII_AllSky'
outputFileName = 'ACESII_AllSky_newProject'



# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import math
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import pyIGRF
from matplotlib import pyplot as plt, animation
from scipy.io import readsav


def AllSkyTrajecMovie(wSite, justPrintSiteNames, rocketFolderPath):

    # --- load attributes for ACESII traj data ---
    rocketAttrs, b, c = ACESII.ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID

    allSkySiteFolder = f'{rocketFolderPath}\\{inputPath_modifier_AllSky}'
    sitePaths = glob(f'{allSkySiteFolder}\\*')
    sites = [path.replace(f'{allSkySiteFolder}\\','') for path in sitePaths]
    wSiteName = sites[wSite]
    pathToSite = sitePaths[wSite]
    WLengths = ['5577', '6300']

    if justPrintSiteNames:
        for i, file in enumerate(sites):
            print('[{:.0f}] {:80s}'.format(i, sites[i]))
        return


    #############################
    # --- get the input files ---
    #############################

    # --- get cal files and convert to python---
    calFiles = [readsav(glob(pathToSite+'\\5577\\*.dat*')[0]), readsav(glob(pathToSite + '\\6300\\*.dat*')[0])]

    # --- get All Sky photos ---
    photoFiles = [glob(pathToSite + '\\5577\\*.png*'), glob(pathToSite + '\\6300\\*.png*')]

    # --- traj Data ---
    trajFolderPath = f'{DataPaths.ACES_data_folder}trajectories\\'
    inputFilesTraj = [glob(trajFolderPath + rf'{ACESII.fliers[0]}\\*_ILat_ILong*')[0],
                     glob(trajFolderPath + rf'{ACESII.fliers[1]}\\\\*_ILat_ILong*')[0]]

    print('\n')
    print(stl.color.UNDERLINE + f'Plotting allSky Data for {wSiteName} data' + stl.color.END)

    # --- get the data from the tmCDF file ---
    stl.prgMsg(f'Loading ACESII traj data')
    data_dicts_traj = []
    for i in range(2):
        data_dict_traj = stl.loadDictFromFile(inputFilesTraj[i])
        data_dict_traj['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_traj['Epoch_esa'][0][i]) for i in (range(len(data_dict_traj['Epoch_esa'][0])))])
        data_dicts_traj.append(data_dict_traj)
    stl.Done(start_time)

    ############################################
    # --- COLLECT IMAGE FILES AND TIMESTAMPS ---
    ############################################

    # get the image time series and the data itself into single variables
    Epoch_AllSky = [[], []]
    imageData = [[], []]
    stl.prgMsg('Collecting Image Data')
    for i in range(len(photoFiles)):
        for imageStr in photoFiles[i]:

            # get the timestamp
            strTimeStamp = imageStr.replace(f'{pathToSite}\\','').replace(f'{WLengths[i]}\\','').replace('skn4_','').replace(f'_{WLengths[i]}_cal.png','')
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

    ###################################
    # --- prepare data for plotting ---
    ###################################
    allGlats = [np.array(deepcopy(calFiles[0]['glats']))[::-1], np.array(deepcopy(calFiles[1]['glats']))[::-1]]
    allGLongs = [np.array(deepcopy(calFiles[0]['glons']))[::-1], np.array(deepcopy(calFiles[1]['glons']))[::-1]]
    allElevs = [np.array(deepcopy(calFiles[0]['elevs'])), np.array(deepcopy(calFiles[1]['elevs']))]
    allImages = imageData
    EpochRocket = [data_dicts_traj[0]['Epoch_esa'][0], data_dicts_traj[1]['Epoch_esa'][0]]
    geoAlt = [data_dicts_traj[0]['geoAlt'][0], data_dicts_traj[1]['geoAlt'][0]]
    geoLat = [data_dicts_traj[0]['geoLat'][0], data_dicts_traj[1]['geoLat'][0]]
    geoLong = [data_dicts_traj[0]['geoLong'][0], data_dicts_traj[1]['geoLong'][0]]

    stl.prgMsg('Collecting and processing/cleaning Image data')

    # -- collect and process all aurora image data ---
    # remove nan values from data and replace with garbage AND convert all Images into Rayleighs
    # description: the image data is a 16-bit counts value normalized by 65535.
    # Invert this and use the calibration factor of 1 R/count given in the cal file

    for i in range(2):  # wavelength

        # --- correct the calibration data ---
        for j in range(len(allGlats[i])):  # array rows
            for k in range(len(allGlats[i][j])):  # row values
                if math.isnan(allGlats[i][j][k]):
                    allGlats[i][j][k] = 70
                    for a in range(len(allImages[i])):  # correct this j,k point in all auroral images
                        allImages[i][a][j][k] = np.nan

                if math.isnan(allGLongs[i][j][k]):
                    allGLongs[i][j][k] = 20
                    for a in range(len(allImages[i])):  # correct this j,k point in all auroral images
                        allImages[i][a][j][k] = np.nan

                if allElevs[i][j][k] <= elevlimits[i] or math.isnan(allElevs[i][j][k]):
                    for a in range(len(allImages[i])):  # correct this j,k point in all auroral images
                        allImages[i][a][j][k] = np.nan

        # --- convert images to rayleighs ---
        for a in range(len(allImages[i])): #number of images
            for j in range(len(allImages[i][a])): # arrays in particular image
                for k in range(len(allImages[i][a][j])): # for each value in each array of image
                    if not math.isnan(allImages[i][a][j][k]):
                        allImages[i][a][j][k] = int(allImages[i][a][j][k]*architecture)

    stl.Done(start_time)

    # --- determine the timestamps for each photo and rocket data ---
    Epoch_AllSky_tt2000 = [[pycdf.lib.datetime_to_tt2000(time) for time in Epoch_AllSky[0]],
                           [pycdf.lib.datetime_to_tt2000(time) for time in Epoch_AllSky[1]]]

    # for each image timestamp, find the index of the closest Epoch_esa in the High flyer's epoch
    imageIndiciesToRocket = [[], []]
    for i in range(len(Epoch_AllSky_tt2000)):
        for stamp in Epoch_AllSky_tt2000[i]:
            imageIndiciesToRocket[i].append(np.abs(data_dicts_traj[0]['Epoch_esa'][0] - stamp).argmin())

    imageIndiciesToRocket = np.array(imageIndiciesToRocket)

    # fill in imageIndiciesToRocket with indicies to match the size of High Flyer Epoch data
    imgIndicies = [[], []]

    for i in range(2):
        for j in range(len(imageIndiciesToRocket[i])):  # loop through the various images
            if j == 0:

                for k in range(imageIndiciesToRocket[i][j + 1]):
                    imgIndicies[i].append(j)

            elif j == len(imageIndiciesToRocket[i]) - 1:  # if you're at the last image index
                for k in range(len(EpochRocket[0]) - imageIndiciesToRocket[i][j]):
                    imgIndicies[i].append(j)

            else:
                for k in range((imageIndiciesToRocket[i][j + 1] - imageIndiciesToRocket[i][j])):
                    imgIndicies[i].append(j)

    imgIndicies = np.array(imgIndicies)
    EpochMovie = np.array([pycdf.lib.tt2000_to_datetime(EpochRocket[0][i]) for i in range(len(EpochRocket[0]))])

    ###############################
    # --- EXTEND LOW FLYER DATA ---
    ###############################

    # --- extend Low Flyer Rocket data to be the same length as High flyer in the beginning and end ---
    highFlyerSize, lowFlyerSize = len(geoAlt[0]), len(geoAlt[1])

    # --- Append start Values to  ---
    no_of_points_start = int((EpochRocket[1][0] - EpochRocket[0][0]) / (rocketAttrs.MinorFrameTime))
    newAlt = [geoAlt[1][0] for i in range(no_of_points_start)]
    newLat = [geoLat[1][0] for i in range(no_of_points_start)]
    newLong = [geoLong[1][0] for i in range(no_of_points_start)]

    for i in range(len(geoAlt[1])):
        newAlt.append(geoAlt[1][i])
        newLat.append(geoLat[1][i])
        newLong.append(geoLong[1][i])

    # --- Append the ending values ---
    remainingIndicies = highFlyerSize - (lowFlyerSize + no_of_points_start)

    for i in range(remainingIndicies):
        newAlt.append(geoAlt[1][-1])
        newLat.append(geoLat[1][-1])
        newLong.append(geoLong[1][-1])

    geoAlt[1], geoLat[1], geoLong[1] = np.array(newAlt), np.array(newLat), np.array(newLong)

    ############################################
    # --- calculate the ILat/ILong variables ---
    ############################################

    stl.prgMsg('Calculating ILat/ILong')

    # -- Output order forpyIGRF.igrf_value ---
    # [0] Declination (+ E | - W)
    # [1] Inclination (+ D | - U)
    # [2] Horizontal Intensity
    # [3] North Comp (+ N | - S)
    # [4] East Comp (+ E | - W)
    # [5] Vertical Comp (+ D | - U)
    # [6] Total Field

    # re-create the geoILat,geoIlong data
    data_dict_IonoProj = {
        'geoILat_RED':[[],[]],
        'geoILong_RED': [[],[]],
        'geoILat_GREEN': [[],[]],
        'geoILong_GREEN': [[],[]],
        'geoILat_AltvsLat_430km':[[],[]],
        'geoILat_AltvsLat_0km': [[], []]
    }
    colors = ['GREEN','RED']

    lat_to_meter = 111.319488  # 1 deg latitude to kilometers on Earth
    def long_to_meter(lat):
        return 111.319488 * math.cos(lat * math.pi / 180)

    date = 2022 + 323 / 365  # Corresponds to 11/20/2022


    for j in range(2): # WAVELENGTH
        colorWL = colors[j]

        for i in range(2): # ROCKET

            for tme in range(len(geoAlt[i])):
                B = pyIGRF.igrf_value(geoLat[i][tme], geoLong[i][tme], geoAlt[i][tme], date)

                data_dict_IonoProj[f'geoILat_{colorWL}'][i].append(
                    geoLat[i][tme]+((geoAlt[i][tme] - projectionAltitude[j])*np.abs(B[3]/B[5]))/lat_to_meter
                )

                data_dict_IonoProj[f'geoILong_{colorWL}'][i].append(
                    geoLong[i][tme] + ((geoAlt[i][tme] - projectionAltitude[j]) * np.abs(B[4] / B[5])) / long_to_meter(geoLong[i][tme])
                )

                if j == 0:
                    data_dict_IonoProj['geoILat_AltvsLat_430km'][i].append(
                        geoLat[i][tme] + ((geoAlt[i][tme] - 430) * np.abs(B[3] / B[5])) / lat_to_meter
                    )

                    data_dict_IonoProj['geoILat_AltvsLat_0km'][i].append(
                        geoLat[i][tme] + ((geoAlt[i][tme] - 0) * np.abs(B[3] / B[5])) / lat_to_meter
                    )


    stl.Done(start_time)

    # Format the IlatILong data to be plotted on the top plot: ALTvsILAT
    # output is a set of x,y points, Format: [ [[x, geolat] [0 km, geoAlt]], ], ]

    projectB = [[], []]

    for i in range(2):
        for j in range(len(geoAlt[i])):
            projectB[i].append(
                [
                    [data_dict_IonoProj['geoILat_AltvsLat_0km'][i][j], data_dict_IonoProj['geoILat_AltvsLat_430km'][i][j]], [0,430]
                ]
            )


    ###########################
    # --- PLOT ALLSKY MOVIE ---
    ###########################
    if createAllSkyMovie:

        ##########################################
        # --- ASSIGN IMAGE TIMESTAMP TO ROCKET ---
        ##########################################

        # --- determine the timestamps for each photo and rocket data ---
        # Note: The image timestamps are taken at the beginning of the image's collection period (30seconds),
        # so we will adjust the time tag to be the middle of the integration_tad_files period by adding 15seconds
        Epoch_AllSky_tt2000 = [[pycdf.lib.datetime_to_tt2000(time)+15E9 for time in Epoch_AllSky[0]],
                               [pycdf.lib.datetime_to_tt2000(time)+15E9 for time in Epoch_AllSky[1]]]

        # for each image timestamp, find the index of the closest Epoch_esa in the High flyer's epoch
        imageIndiciesToRocket = [[], []]
        for i in range(len(Epoch_AllSky_tt2000)):
            for stamp in Epoch_AllSky_tt2000[i]:
                imageIndiciesToRocket[i].append(np.abs(data_dicts_traj[0]['Epoch_esa'][0] - stamp).argmin())

        imageIndiciesToRocket = np.array(imageIndiciesToRocket)

        # fill in imageIndiciesToRocket with indicies to match the size of High Flyer Epoch data
        imgIndicies = [[], []]

        for i in range(2):
            for j in range(len(imageIndiciesToRocket[i])): # loop through the various images
                if j == 0:

                    for k in range(imageIndiciesToRocket[i][j+1]):
                        imgIndicies[i].append(j)

                elif j == len(imageIndiciesToRocket[i])-1: # if you're at the last image index
                    for k in range(len(EpochRocket[0]) - imageIndiciesToRocket[i][j]):
                        imgIndicies[i].append(j)

                else:
                    for k in range((imageIndiciesToRocket[i][j+1] - imageIndiciesToRocket[i][j]) ):
                        imgIndicies[i].append(j)

        imgIndicies = np.array(imgIndicies)
        EpochMovie = np.array([pycdf.lib.tt2000_to_datetime(EpochRocket[0][i]) for i in range(len(EpochRocket[0])) ])

        # Target time of interest
        targTime = dt.datetime(2022,11,20,17,25,000)

        # determine target time index for trajectory - High Flyer
        # HF_targTime_index = np.abs(EpochMovie - targTime).argmin()
        HF_targTime_index = 0

        # determine target time index for trajectory - Low Flyer
        # LF_targTime_index = np.abs(EpochMovie - targTime).argmin()
        LF_targTime_index = 0

        # determine target time index for all skies - 5570
        # image_tarTime_5570 = imgIndicies[0][HF_targTime_index]
        image_tarTime_5570 = 0

        # determine target time index for all skies - 6300
        # image_tarTime_6300 = imgIndicies[1][HF_targTime_index]
        image_tarTime_6300 = 0

        ##########################
        # --- INITIALIZE PLOTS ---
        ##########################
        stl.prgMsg('Initializing Allsky Plots')

        # --- prepare map information ---
        projPC = ccrs.PlateCarree()  # MUST be kept on the set_extent, crs =, command AND pcolormesh transform command
        lonW = 12
        lonE = 18
        latS = 68
        latN = 74.5
        cLat = (latN + latS) / 2
        cLon = (lonW + lonE) / 2
        res = '10m'
        # projType = ccrs.Mollweide(central_longitude=cLon, central_latitude=cLat)
        projType = ccrs.Stereographic(central_longitude=cLon, central_latitude=cLat)

        allSkyCmap = 'inferno'
        LowFlyerColor = 'royalblue'
        HighFlyerColor = 'red'
        LowFlyerProjectionColor = 'cyan'
        HighFlyerProjectionColor = 'darkorange'

        AltvsLatLineWidth = 4
        LatvsLongLineWidth = 4
        AllSkyLineWidth = 4

        # ---------------------
        # --- START PLOTTING ---
        # ---------------------

        # figure size
        figure_height = 20
        figure_width = 36
        plt.style.use('dark_background')
        fig = plt.figure(dpi=60)
        fig.set_figwidth(figure_width)
        fig.set_figheight(figure_height)

        # The overall plot shape
        nRows = 20
        nCols = 51 # MUST BE EVEN VALUE
        smallPlotSize = 4
        verticalSpacing = 3 #spacing between the top row and bottom row plots
        gs = gridspec.GridSpec(nrows=nRows, ncols=nCols, figure=fig)

        # define the axes
        axAlt = fig.add_subplot(gs[0:smallPlotSize, 0:51])
        ax5577 = fig.add_subplot(gs[smallPlotSize+verticalSpacing:nRows, 0:22+1], projection=projType)
        ax5577cbar = fig.add_subplot(gs[smallPlotSize+verticalSpacing:nRows, 23])
        ax6300 = fig.add_subplot(gs[smallPlotSize+verticalSpacing:nRows, 27:50], projection=projType)
        ax6300cbar = fig.add_subplot(gs[smallPlotSize+verticalSpacing:nRows, 50])

        # --- initialize the title ---
        supTitle = fig.suptitle(f'ACESII\n'
                                f'{EpochMovie[0]}\n', fontsize=40,color='white')

        # --- plot the norwegian map data ---
        # 5577A
        ax5577.set_title('GREEN 557.7 nm (150 km)', fontsize=30)
        gl5577 = ax5577.gridlines(draw_labels=True, linewidth=3, alpha=0.4, linestyle='--')
        gl5577.xlabel_style = {'size': 20, 'color': 'white', 'weight': 'bold'}
        gl5577.ylabel_style = {'size': 20, 'color': 'white', 'weight': 'bold'}
        ax5577.set_extent([lonW, lonE, latS, latN], crs=projPC)  # controls lat/long axes display
        ax5577.coastlines(resolution=res, color='white', alpha=0.8)  # adds coastlines with resolution
        ax5577.set_aspect(0.3)

        # 6300A
        ax6300.set_title('RED 630.0 nm (250 km)', fontsize=30)
        gl6300 = ax6300.gridlines(draw_labels=True, linewidth=3, color='white', alpha=0.4, linestyle='--')
        gl6300.xlabel_style = {'size': 20, 'color': 'white', 'weight': 'bold'}
        gl6300.ylabel_style = {'size': 20, 'color': 'white', 'weight': 'bold'}
        ax6300.set_extent([lonW, lonE, latS, latN], crs=projPC)  # controls lat/long axes display
        ax6300.coastlines(resolution=res, color='white', alpha=0.8)  # adds coastlines with resolution
        ax6300.set_aspect(0.3)

        # --- initialize the All Sky Image ---
        cbarVmin = 0
        cbarVmax = 16000
        # 5577A
        # cmap5577 = ax5577.pcolormesh(allGLongs[0], allGlats[0], allImages[0][0], cmap=allSkyCmap, transform=projPC,vmin=cbarVmin,vmax=cbarVmax)
        cmap5577 = ax5577.pcolormesh(allGLongs[0], allGlats[0], allImages[0][image_tarTime_5570], cmap='viridis', transform=projPC, vmin=cbarVmin, vmax=cbarVmax)
        plt.colorbar(mappable=cmap5577, cax=ax5577cbar, orientation='vertical',fraction=0.046, pad=0.04)
        ax5577cbar.set_ylabel('Brightness (Rayleigh)',fontsize=20)
        ax5577cbar.tick_params(labelsize=20)

        # 6300A
        cmap6300 = ax6300.pcolormesh(allGLongs[1], allGlats[1], allImages[1][image_tarTime_6300], cmap=allSkyCmap, transform=projPC,vmin=cbarVmin,vmax=cbarVmax)
        plt.colorbar(mappable=cmap6300, cax=ax6300cbar, orientation='vertical',fraction=0.046, pad=0.04)
        ax6300cbar.set_ylabel('Brightness (Rayleigh)',fontsize=20)
        ax6300cbar.tick_params(labelsize=20)

        # --- initialize total trajectory on the All Sky Image ---

        # --- GREEN ---
        # plot projected B on 5577 Allsky
        # ax5577.plot(data_dict_IonoProj['geoILong_GREEN'][1], data_dict_IonoProj['geoILat_GREEN'][1], color=LowFlyerProjectionColor, transform=projPC, linewidth=AllSkyLineWidth,label='Low Flyer Projected IGRF - ')
        # ax5577.plot(data_dict_IonoProj['geoILong_GREEN'][0], data_dict_IonoProj['geoILat_GREEN'][0], color=HighFlyerProjectionColor, transform=projPC, linewidth=AllSkyLineWidth,label='High Flyer Projected IGRF') if not useOnlyLowFlyer else [0,]

        # plot Ilat/Ilong marker on 5577 on allsky
        # AllSky5577_projected_marker_low = ax5577.scatter(data_dict_IonoProj['geoILong_GREEN'][1][0], data_dict_IonoProj['geoILat_GREEN'][1][0], color=LowFlyerProjectionColor, marker='x', linewidth=10, transform=projPC)
        # AllSky5577_projected_marker_high = ax5577.scatter(data_dict_IonoProj['geoILong_GREEN'][0][0], data_dict_IonoProj['geoILat_GREEN'][0][0], color=HighFlyerProjectionColor, marker='x', linewidth=10, transform=projPC) if not useOnlyLowFlyer else [0,]

        # plot on 5577 Allsky
        ax5577.plot(geoLong[1], geoLat[1], color=LowFlyerColor,linewidth=4, transform=projPC, label='Low Flyer Trajectory')
        ax5577.plot(geoLong[0], geoLat[0], color=HighFlyerColor,linewidth=4, transform=projPC,label='High Flyer Trajectory') if not useOnlyLowFlyer else [0,]

        # plot marker on 5577 Allsky
        # AllSky5577_marker_low = ax5577.scatter(geoLong[1][0], geoLat[1][0], color=LowFlyerColor, marker='x', linewidth=10, transform=projPC)
        # AllSky5577_marker_high = ax5577.scatter(geoLong[0][0], geoLat[0][0], color=HighFlyerColor, marker='x', linewidth=10, transform=projPC) if not useOnlyLowFlyer else [0,]
        AllSky5577_marker_low = ax5577.scatter(geoLong[1][LF_targTime_index], geoLat[1][LF_targTime_index], color=LowFlyerColor, marker='x', linewidth=10, transform=projPC)
        AllSky5577_marker_high = ax5577.scatter(geoLong[0][HF_targTime_index], geoLat[0][HF_targTime_index], color=HighFlyerColor, marker='x', linewidth=10, transform=projPC) if not useOnlyLowFlyer else [0, ]

        # --- RED ---
        # plot projected B on 6300 Allsky
        # ax6300.plot(data_dict_IonoProj['geoILong_RED'][1], data_dict_IonoProj['geoILat_RED'][1], color=LowFlyerProjectionColor, transform=projPC, linewidth=AllSkyLineWidth)
        # ax6300.plot(data_dict_IonoProj['geoILong_RED'][0], data_dict_IonoProj['geoILat_RED'][0], color=HighFlyerProjectionColor, transform=projPC, linewidth=AllSkyLineWidth) if not useOnlyLowFlyer else [0,]

        # plot Ilat/Ilong marker on 6300 on allsky
        # AllSky6300_projected_marker_low = ax6300.scatter(data_dict_IonoProj['geoILong_RED'][1][0], data_dict_IonoProj['geoILat_RED'][1][0], color=LowFlyerProjectionColor, marker='x', linewidth=10, transform=projPC)
        # AllSky6300_projected_marker_high = ax6300.scatter(data_dict_IonoProj['geoILong_RED'][0][0], data_dict_IonoProj['geoILat_RED'][0][0], color=HighFlyerProjectionColor, marker='x', linewidth=10, transform=projPC) if not useOnlyLowFlyer else [0,]

        # plot on 6300 Allsky
        ax6300.plot(geoLong[1], geoLat[1], color=LowFlyerColor,linewidth=4, transform=projPC)
        ax6300.plot(geoLong[0], geoLat[0], color=HighFlyerColor,linewidth=4, transform=projPC) if not useOnlyLowFlyer else [0,]

        # plot marker on 6300 Allsky
        AllSky6300_marker_low = ax6300.scatter(geoLong[1][LF_targTime_index], geoLat[1][LF_targTime_index], color=LowFlyerColor, linewidth=10, marker='x', transform=projPC)
        AllSky6300_marker_high = ax6300.scatter(geoLong[0][HF_targTime_index], geoLat[0][HF_targTime_index], color=HighFlyerColor, linewidth=10, marker='x', transform=projPC) if not useOnlyLowFlyer else [0,]
        # AllSky6300_marker_low = ax6300.scatter(geoLong[1][0], geoLat[1][0], color=LowFlyerColor, linewidth=10,
        #                                        marker='x', transform=projPC)
        # AllSky6300_marker_high = ax6300.scatter(geoLong[0][0], geoLat[0][0], color=HighFlyerColor, linewidth=10,
        #                                         marker='x', transform=projPC) if not useOnlyLowFlyer else [0, ]

        # -----------------------------------------------------------------------------------
        # stylization parameters of text labels in alt vs lat and lat vs long plots
        textOffset_alt = 0.15
        textOffset_latlong = [[-0.15, -0.15], [0.15, 0.25]]
        text_alignment = ['left', 'right']
        textUTC_style_lat = [dict(size=20, color=HighFlyerColor), dict(size=20, color=LowFlyerColor)]
        textUTC_style_project_lat = [dict(size=20, color=HighFlyerProjectionColor), dict(size=20, color=LowFlyerProjectionColor)]
        textUTC_style_alt = [dict(size=20, color=HighFlyerColor), dict(size=20, color=LowFlyerColor)]
        altrounding = 1  # how many decimals to round to
        # -----------------------------------------------------------------------------------

        # --- initialize Altitude vs Lat Plot ---
        axAlt.set_ylabel('Altitude [km]', fontsize=20)
        axAlt.set_xlabel('Geo Lat', fontsize=20)
        axAlt.tick_params(axis='x', which='both', labelsize=25)
        axAlt.tick_params(axis='y', which='both', labelsize=25)
        axAlt.axhline(projectionAltitude[0], color='limegreen', linestyle='--', alpha=0.5, linewidth=AltvsLatLineWidth)
        axAlt.axhline(projectionAltitude[1], color='limegreen', linestyle='--', alpha=0.5, linewidth=AltvsLatLineWidth)
        axAlt.text(69.05, 170, '150 km', color='limegreen', fontsize=20)
        axAlt.text(69.05, 270, '250 km', color='limegreen', fontsize=20)

        # trajectory
        axAlt.plot(geoLat[0], geoAlt[0], color=HighFlyerColor, linewidth=AltvsLatLineWidth) if not useOnlyLowFlyer else [0,]
        axAlt.plot(geoLat[1], geoAlt[1], color=LowFlyerColor, linewidth=AltvsLatLineWidth)

        # B-projection
        legendAltLatLabel = f'Projected IGRF'
        BprojectionHigh, = axAlt.plot(projectB[0][HF_targTime_index][0], projectB[0][HF_targTime_index][1], color='white', alpha=0.7, linestyle='--') if not useOnlyLowFlyer else [0,]  # high flyer
        BprojectionLow, = axAlt.plot(projectB[1][LF_targTime_index][0], projectB[1][LF_targTime_index][1], color='white', alpha=0.7, linestyle='--', label=legendAltLatLabel)  # low flyer

        # legend for AltvsLat projection
        axAlt.legend(loc="upper left", fontsize=15)

        # marker
        Altlat_marker_high = axAlt.scatter(geoLat[0][HF_targTime_index], geoAlt[0][HF_targTime_index], color=HighFlyerColor, marker='x', linewidth=15) if not useOnlyLowFlyer else [0,]
        Altlat_marker_low = axAlt.scatter(geoLat[1][LF_targTime_index], geoAlt[1][LF_targTime_index], color=LowFlyerColor, marker='x', linewidth=15)

        # text labels
        Altlat_text_high = axAlt.text(geoLat[0][HF_targTime_index] - textOffset_alt, geoAlt[0][HF_targTime_index],
                                   f'{round(geoAlt[0][HF_targTime_index], altrounding)} km', ha=text_alignment[1],
                                   **textUTC_style_alt[0]) if not useOnlyLowFlyer else [0,]
        Altlat_text_low = axAlt.text(geoLat[1][LF_targTime_index] - textOffset_alt, geoAlt[1][LF_targTime_index],
                                  f'{round(geoAlt[1][LF_targTime_index], altrounding)} km', ha=text_alignment[1],
                                  **textUTC_style_alt[1])

        # --- intialize Trajectory ILat vs Ilong projection ---

        # plot the lat/long coordinates that update on the 5577 graph
        # IlatILong_text_low_GREEN = ax5577.text(12, 73.25, f'({round(data_dict_IonoProj["geoILong_GREEN"][1][0], altrounding)}' + '$^{\circ}$' + f', {round(data_dict_IonoProj["geoILat_GREEN"][1][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_project_lat[1], transform=projPC)
        # IlatILong_text_high_GREEN = ax5577.text(12.05, 73, f'({round(data_dict_IonoProj["geoILong_GREEN"][0][0], altrounding)}' + '$^{\circ}$' + f', {round(data_dict_IonoProj["geoILat_GREEN"][0][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_project_lat[0],transform=projPC) if not useOnlyLowFlyer else [0, ]
        # latLong_text_low_GREEN = ax5577.text(12.1, 72.75, f'({round(geoLong[1][0], altrounding)}' + '$^{\circ}$' + f', {round(geoLat[1][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_lat[1], transform=projPC)
        # latLong_text_high_GREEN = ax5577.text(12.15, 72.5, f'({round(geoLong[0][0], altrounding)}' + '$^{\circ}$' + f', {round(geoLat[0][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_lat[0],transform=projPC) if not useOnlyLowFlyer else [0, ]

        ax5577.legend(loc='upper left', prop={'size': 14})

        # plot the lat/long coordinates that update on the 6300 graph
        # IlatILong_text_low_RED = ax6300.text(12.5, 70.25, f'({round(data_dict_IonoProj["geoILong_RED"][1][0], altrounding)}' + '$^{\circ}$' + f', {round(data_dict_IonoProj["geoILat_RED"][1][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_project_lat[1], transform=projPC)
        # IlatILong_text_high_RED = ax6300.text(12.55, 70, f'({round(data_dict_IonoProj["geoILong_RED"][0][0], altrounding)}' + '$^{\circ}$' + f', {round(data_dict_IonoProj["geoILat_RED"][0][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_project_lat[0], transform=projPC) if not useOnlyLowFlyer else [0, ]
        # latLong_text_low_RED = ax6300.text(12.6, 69.75, f'({round(geoLong[1][0], altrounding)}' + '$^{\circ}$' + f', {round(geoLat[1][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_lat[1], transform=projPC)
        # latLong_text_high_RED = ax6300.text(12.65, 69.5, f'({round(geoLong[0][0], altrounding)}' + '$^{\circ}$' + f', {round(geoLat[0][0], altrounding)}' + '$^{\circ}$)', ha='center', **textUTC_style_lat[0], transform=projPC) if not useOnlyLowFlyer else [0, ]

        # ax6300.legend(loc='upper left', prop={'size': 14})

        stl.Done(start_time)

        # fig.savefig(r'C:\Users\cfelt\OneDrive\Desktop\TRACERS_at_SSL_trajectory.png')
        # plt.show()

        ###########################
        # --- ANIMATE THE MOVIE ---
        # ###########################
        def animatePlot(i):
            print('', end='\r' + stl.color.RED + f'{round(100 * i / len(EpochMovie), 1)} %' + stl.color.END)

            # update the title
            supTitle.set_text(f'ACESII\n'
                              f'{EpochMovie[i]}\n')

            # --- ALL SKY ---
            # update all sky images
            wImage5577 = imgIndicies[0][i]
            wImage6300 = imgIndicies[1][i]

            cmap5577.set_array(allImages[0][imgIndicies[0][i]].ravel())
            cmap6300.set_array(allImages[1][imgIndicies[1][i]].ravel())

            # update the allsky5577 rocket marker positions
            AllSky5577_marker_high.set_offsets([geoLong[0][i], geoLat[0][i]]) if not useOnlyLowFlyer else [0, ]
            AllSky5577_marker_low.set_offsets([geoLong[1][i], geoLat[1][i]])
            # AllSky5577_projected_marker_high.set_offsets([data_dict_IonoProj['geoILong_GREEN'][0][i], data_dict_IonoProj['geoILat_GREEN'][0][i]]) if not useOnlyLowFlyer else [0, ]
            # AllSky5577_projected_marker_low.set_offsets([data_dict_IonoProj['geoILong_GREEN'][1][i], data_dict_IonoProj['geoILat_GREEN'][1][i]])

            # update the allsky6300 rocket marker positions
            AllSky6300_marker_high.set_offsets([geoLong[0][i], geoLat[0][i]]) if not useOnlyLowFlyer else [0, ]
            AllSky6300_marker_low.set_offsets([geoLong[1][i], geoLat[1][i]])
            # AllSky6300_projected_marker_high.set_offsets([data_dict_IonoProj['geoILong_RED'][0][i], data_dict_IonoProj['geoILat_RED'][0][i]]) if not useOnlyLowFlyer else [0, ]
            # AllSky6300_projected_marker_low.set_offsets([data_dict_IonoProj['geoILong_RED'][1][i], data_dict_IonoProj['geoILat_RED'][1][i]])

            # --- ALT VS LAT ---
            # update Alt vs Lat marker
            Altlat_marker_high.set_offsets([geoLat[0][i], geoAlt[0][i]]) if not useOnlyLowFlyer else [0, ]
            Altlat_marker_low.set_offsets([geoLat[1][i], geoAlt[1][i]])

            # update B-project lines on Alt vs lat plot
            BprojectionHigh.set_xdata(projectB[0][i][0]) if not useOnlyLowFlyer else [0, ]
            BprojectionHigh.set_ydata(projectB[0][i][1]) if not useOnlyLowFlyer else [0, ]
            BprojectionLow.set_xdata(projectB[1][i][0])
            BprojectionLow.set_ydata(projectB[1][i][1])

            # update Alt vs lat text
            Altlat_text_high.set_x(geoLat[0][i]) if not useOnlyLowFlyer else [0, ]
            Altlat_text_high.set_y(geoAlt[0][i]) if not useOnlyLowFlyer else [0, ]
            Altlat_text_high.set_text(f'{round(geoAlt[0][i], altrounding)} km') if not useOnlyLowFlyer else [0, ]

            Altlat_text_low.set_x(geoLat[1][i])
            Altlat_text_low.set_y(geoAlt[1][i])
            Altlat_text_low.set_text(f'{round(geoAlt[1][i], altrounding)} km')

            # --- GREEN ---
            # update lat vs long text
            # latLong_text_high_GREEN.set_text(f'({round(geoLong[0][i], altrounding)}' + '$^{\circ}$' +
            #                            f', {round(geoLat[0][i], altrounding)}' + '$^{\circ}$)') if not useOnlyLowFlyer else [0, ]
            # latLong_text_low_GREEN.set_text(f'({round(geoLong[1][i], altrounding)}' + '$^{\circ}$' +
            #                            f', {round(geoLat[1][i], altrounding)}' + '$^{\circ}$)')

            # update Ilat vs Ilong text
            # IlatILong_text_high_GREEN.set_text(f'({round(data_dict_IonoProj["geoILong_GREEN"][0][i], altrounding)}' + '$^{\circ}$' +
            #                            f', {round(data_dict_IonoProj["geoILat_GREEN"][0][i], altrounding)}' + '$^{\circ}$)') if not useOnlyLowFlyer else [0,]
            #
            # IlatILong_text_low_GREEN.set_text(f'({round(data_dict_IonoProj["geoILong_GREEN"][1][i], altrounding)}' + '$^{\circ}$' +
            #                           f', {round(data_dict_IonoProj["geoILat_GREEN"][1][i], altrounding)}' + '$^{\circ}$)')

            # --- RED ---
            # update lat vs long text
            # latLong_text_high_RED.set_text(f'({round(geoLong[0][i], altrounding)}' + '$^{\circ}$' +
            #                                  f', {round(geoLat[0][i], altrounding)}' + '$^{\circ}$)') if not useOnlyLowFlyer else [0, ]
            # latLong_text_low_RED.set_text(f'({round(geoLong[1][i], altrounding)}' + '$^{\circ}$' +
            #                                 f', {round(geoLat[1][i], altrounding)}' + '$^{\circ}$)')

            # update Ilat vs Ilong text
            # IlatILong_text_high_RED.set_text(f'({round(data_dict_IonoProj["geoILong_RED"][0][i], altrounding)}' + '$^{\circ}$' +
            #                                    f', {round(data_dict_IonoProj["geoILat_RED"][0][i], altrounding)}' + '$^{\circ}$)') if not useOnlyLowFlyer else [0, ]
            #
            # IlatILong_text_low_RED.set_text(f'({round(data_dict_IonoProj["geoILong_RED"][1][i], altrounding)}' + '$^{\circ}$' +
            #                                   f', {round(data_dict_IonoProj["geoILat_RED"][0][i], altrounding)}' + '$^{\circ}$)')

        stl.prgMsg('Creating AllSky Movie')

        if plotSpecificLocations:
            locations = [i for i in specificLocations]
        elif doFrameSkips:
            locations = [i for i in range(0, len(EpochMovie), frame_skips)]  # NEEDS TO BE THE HIGH FLYER LENGTH
        else:
            locations = [i for i in range(len(EpochMovie))]


        anim = animation.FuncAnimation(fig=fig, func=animatePlot, interval=1000 / fps, frames=locations)

        # need a .mp4 writer, following code points to it
        writervideo = animation.FFMpegWriter(fps=fps)
        import matplotlib as mpl
        mpl.rcParams['animation.ffmpeg_path'] = r'C:\Users\cfelt\PycharmProjects\UIOWA_CDF_operator\ACESII_code\supportCode\ffmpeg\bin\ffmpeg.exe'
        anim.save(rf'C:\Data\ACESII\trajectories\trajectory_plots\movies\{outputFileName}.mp4',writer=writervideo)


        stl.Done(start_time)








# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

rocketFolderPath = DataPaths.ACES_data_folder

if justPrintSiteNames:
    AllSkyTrajecMovie(wSite, justPrintSiteNames,rocketFolderPath)
else:
    AllSkyTrajecMovie(wSite, justPrintSiteNames,rocketFolderPath)
