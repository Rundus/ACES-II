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
inputPath_modifier_AllSky = 'all_sky\skibotn' # e.g. 'L1' or 'L1'. It's the name of the broader input folder inside data\ACESII
inputPath_modifier_attitude = 'attitude'


# reduce the data to speed up the code in order to help with debugging
fracReduction = 1

# --- Plot/Anime AllSky MOVIE ---
Plot_AllSkyMovie = True # MUST BE TRUE TO ANIMATE MOVIE
projectionAltitude = [150, 250] # in km. Format: [green, red]

Animate_AllSkyMovie = True
fps = 20 # fps of the video

plotSpecificLocations = False # plot specific locations
specificLocations = [0 + 500*i for i in range(29)]

doFrameSkips = True # plot frame skips
frame_skips = 20

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import spaceToolsLib as stl
from matplotlib import pyplot as plt, animation
# --- --- --- ---  --- --
# --- PLOT PARAMETERS ---
# --- --- --- ---  --- --

# Figure toggles
dpi = 60
figure_height = 16
figure_width = 30

# formatting toggles
allSkyCmap = stl.matlab_parula_cmap()
HighFlyerColor = 'tab:red'
LowFlyerColor = 'tab:orange'
LowFlyerProjectionColor = 'cyan'
HighFlyerProjectionColor = 'darkorange'
plt.rcParams["font.family"] = "Arial"

# Projection toggles
central_longitude = 15
central_latitude = 70
lonW = 5
lonE = 25
latS = 67
latN = 75
cLat = (latN + latS) / 2
cLon = (lonW + lonE) / 2
res = '50m'

# Alt_Lat_Plot
AltvsLat_LineWidth = 4
AltvsLat_scatterLineWidth = 15
AltvsLat_textOffset = 0.15
AltvsLat_text_alignment = ['left', 'right']
AltvsLat_textUTC_style = [dict(size=20, color=HighFlyerColor), dict(size=20, color=LowFlyerColor)]
AltvsLat_textUTC_style_project = [dict(size=20, color=HighFlyerProjectionColor), dict(size=20, color=LowFlyerProjectionColor)]
AltvsLat_textUTC_style = [dict(size=20, color=HighFlyerColor), dict(size=20, color=LowFlyerColor)]
AltvsLat_rounding = 1  # how many decimals to round to
AltvsLat_LineWidth = 4
AltVsLat_xTickSize = 30
AltVsLat_yTickSize = 20
AltvsLat_LabelFontsize = 30
AltvsLat_textFontSize = 20
AltvsLat_legendFontsize = 25
AltvsLat_Bproject_low = 0# in kilometers
AltvsLat_Bproject_High = 600 # in kilometers

# AllSky Plot toggles
AllSky_LineWidth = 4
AllSky_gridLineWidth = 5
AllSky_LabelSize = 26
AllSky_TitleSize = 45
AllSky_markerSize = 10
AllSky_legendSize = 5
AllSky_TitlePadding = 0

# colorbar toggles
cbar_Vmin = 0
cbar_Vmax = 12
cbar_labelFontSize = 30
cbar_tickFontSize = 20



#######################
# --- BEGIN PROGRAM ---
#######################
def PlotExtra_AllSkyMovie(justPrintSiteNames, rocketFolderPath):

    # --- load attributes for ACESII traj data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID
    allSkySiteFolder = f'{rocketFolderPath}\\{inputPath_modifier_AllSky}'
    sites = [path.replace(f'{allSkySiteFolder}\\','') for path in glob(f'{allSkySiteFolder}\\*')]

    if justPrintSiteNames:
        for i, file in enumerate(sites):
            print('[{:.0f}] {:80s}'.format(i, sites[i]))
        return

    #############################
    # --- get the input files ---
    #############################

    # --- get All Sky data ---
    inputFiles_allsky = [glob(f'{ACES_data_folder}\\{inputPath_modifier_AllSky}\\5577\\*.cdf*')[0], glob(f'{ACES_data_folder}\\{inputPath_modifier_AllSky}\\\\6300\\*.cdf*')[0]]
    data_dicts_allSky = [loadDictFromFile(inputFiles_allsky[0]), loadDictFromFile(inputFiles_allsky[1])]
    Epoch_AllSky = [data_dicts_allSky[0]['Epoch'][0], data_dicts_allSky[1]['Epoch'][0]] # get the image time series and the data itself into single variables
    allImages = [data_dicts_allSky[0]['AllSkyImages'][0], data_dicts_allSky[1]['AllSkyImages'][0]]

    # --- traj Data ---
    inputFiles_Attitude = [glob(f'{rocketFolderPath}\\{inputPath_modifier_attitude}\\{fliers[0]}\\*.cdf')[0], glob(f'{rocketFolderPath}\\{inputPath_modifier_attitude}\\{fliers[1]}\\*.cdf')[0]]

    print('\n')
    print(color.UNDERLINE + f'Creating AllSky Movie' + color.END)

    # --- get the data from the tmCDF file ---
    prgMsg(f'Loading attitude data')
    data_dicts_attitude = [loadDictFromFile(inputFiles_Attitude[0]),loadDictFromFile(inputFiles_Attitude[1])]
    Done(start_time)

    ###################################
    # --- prepare data for plotting ---
    ###################################
    prgMsg('Preparing Data for Plotting')
    allGlats = [data_dicts_allSky[0]['GLats'][0], data_dicts_allSky[1]['GLats'][0]]
    allGLongs = [data_dicts_allSky[0]['GLongs'][0], data_dicts_allSky[1]['GLongs'][0]]
    EpochRocket = [data_dicts_attitude[0]['Epoch'][0], data_dicts_attitude[1]['Epoch'][0]]

    Alt = [data_dicts_attitude[0]['Alt'][0]/1000, data_dicts_attitude[1]['Alt'][0]/1000]
    Lat = [data_dicts_attitude[0]['Lat'][0], data_dicts_attitude[1]['Lat'][0]]
    Long = [data_dicts_attitude[0]['Long'][0], data_dicts_attitude[1]['Long'][0]]

    # reduce the data to speed up the code in order to help with debugging
    for i in range(2):
        EpochRocket[i] = EpochRocket[i][0: int(len(EpochRocket[i])*fracReduction)]
        Alt[i] = Alt[i][0:int(len(Alt[i])*fracReduction)]
        Lat[i] = Lat[i][0:int(len(Lat[i])*fracReduction)]
        Long[i] = Long[i][0:int(len(Long[i])*fracReduction)]

    # --- determine the timestamps for each photo and rocket data ---
    # Note: The image timestamps are taken at the beginning of the image's collection period (30seconds),
    # so we will adjust the time tag to be the middle of the integration_tad_files period by adding 15 seconds
    Epoch_AllSky_centered = [
        [pycdf.lib.tt2000_to_datetime(int(pycdf.lib.datetime_to_tt2000(time) + 15E9)) for time in Epoch_AllSky[0]],
        [pycdf.lib.tt2000_to_datetime(int(pycdf.lib.datetime_to_tt2000(time) + 15E9)) for time in Epoch_AllSky[1]]]

    # for each image timestamp, find the index of the closest Epoch in the High flyer's epoch
    imageIndicies_ToHFRocketEpoch = np.array([[np.abs(EpochRocket[0] - stamp).argmin() for stamp in Epoch_AllSky_centered[i]] for i in range(len(Epoch_AllSky_centered))])
    imageIndicies = [[],[]]

    for i in range(len(imageIndicies)):

        for j,idxEpoch_HF in enumerate(imageIndicies_ToHFRocketEpoch[i]):

            if j == 0:
                for k in range(imageIndicies_ToHFRocketEpoch[i][j+1]):
                    imageIndicies[i].append(j)
            elif j == len(imageIndicies_ToHFRocketEpoch[i])-1:
                for k in range(len(EpochRocket[0]) - imageIndicies_ToHFRocketEpoch[i][j]):
                    imageIndicies[i].append(j)
            else:
                for k in range(imageIndicies_ToHFRocketEpoch[i][j+1] - imageIndicies_ToHFRocketEpoch[i][j]):
                    imageIndicies[i].append(j)

    # from collections import Counter
    # # print(imageIndicies_ToHFRocketEpoch[0])
    # # print(Counter(imageIndicies[0]).keys())
    # # print(Counter(imageIndicies[0]).values())
    # # print(len(imageIndicies[0]))
    # # print(imageIndicies_ToHFRocketEpoch[1])
    # # print(Counter(imageIndicies[1]).keys())
    # # print(Counter(imageIndicies[1]).values())
    # # print(len(imageIndicies[1]))

    Done(start_time)

    ###############################
    # --- EXTEND LOW FLYER DATA ---
    ###############################
    prgMsg('Extending low flyer data')

    # --- extend Low Flyer Rocket data to be the same length as High flyer in the beginning and end ---

    # --- Append start Values to  ---
    no_of_points_start = np.abs(EpochRocket[0] - EpochRocket[1][0]).argmin()
    newAlt = [Alt[1][0] for i in range(no_of_points_start)]
    newLat = [Lat[1][0] for i in range(no_of_points_start)]
    newLong = [Long[1][0] for i in range(no_of_points_start)]
    Alt[1] = newAlt + list(Alt[1])
    Lat[1] = newLat + list(Lat[1])
    Long[1] = newLong + list(Long[1])

    # --- Append the ending values ---
    remainingIndicies = len(EpochRocket[0]) - (len(EpochRocket[1]) + no_of_points_start)
    newAlt = [Alt[1][-1] for i in range(remainingIndicies)]
    newLat = [Lat[1][-1] for i in range(remainingIndicies)]
    newLong = [Long[1][-1] for i in range(remainingIndicies)]
    Alt[1] = list(Alt[1]) + newAlt
    Lat[1] = list(Lat[1]) + newLat
    Long[1] = list(Long[1]) + newLong
    Alt[1], Lat[1], Long[1] = np.array(Alt[1]), np.array(Lat[1]), np.array(Long[1])
    Done(start_time)


    ###########################################
    # --- DETERMINE THE MAGNETIC PROJECTION ---
    ###########################################


    # Get the ENU direction of the magnetic field for both flyers throughout the flight
    prgMsg('Calculating CHAOS model')
    # Get the geomagnetic field direction
    Bgeo = [stl.CHAOS(lat=Lat[0], long=Long[0], alt=Alt[0], times=EpochRocket[0]),
            stl.CHAOS(lat=Lat[1], long=Long[1], alt=Alt[1], times=EpochRocket[0])]

    # Normalize the arrows
    Bgeo_norm = [[vec/np.linalg.norm(vec) for vec in bgeo] for bgeo in Bgeo]
    Done(start_time)

    prgMsg('Determining Bgeo Projection')
    projectB_lats = []
    for wRkt in range(2):
        ILat_highAlt, Ilong_highAlt = stl.ILatILong_Projection(Alt=Alt[wRkt],Lat=Lat[wRkt],Long=Long[wRkt], Zproject=AltvsLat_Bproject_High,b=Bgeo_norm[wRkt])
        ILat_lowAlt, Ilong_lowAlt = stl.ILatILong_Projection(Alt=Alt[wRkt], Lat=Lat[wRkt], Long=Long[wRkt], Zproject=AltvsLat_Bproject_low, b=Bgeo_norm[wRkt])
        projectB_lats.append(np.transpose(np.array([ILat_highAlt, ILat_lowAlt])))

    Done(start_time)

    #--------------------------
    ###########################
    # --- PLOT ALLSKY MOVIE ---
    ###########################
    # --------------------------

    if Plot_AllSkyMovie:

        ##########################
        # --- INITIALIZE PLOTS ---
        ##########################
        prgMsg('Initializing Allsky Plots')

        # ----------------------
        # --- START PLOTTING ---
        # ----------------------

        projProjection = ccrs.Orthographic(central_longitude=central_longitude, central_latitude=central_latitude)
        projTransform = ccrs.PlateCarree()

        # figure size
        plt.style.use('dark_background')
        fig = plt.figure(dpi=dpi)
        fig.set_figwidth(figure_width)
        fig.set_figheight(figure_height)

        # The overall plot shape
        gs0 = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1, 4],hspace=0.2)

        # define the axes
        axAlt = fig.add_subplot(gs0[0])
        gs00 = gs0[1].subgridspec(1,2)
        ax5577 = fig.add_subplot(gs00[0], projection=projProjection)
        ax6300 = fig.add_subplot(gs00[1], projection=projProjection)

        # --- --- --- --- --- --- --- --- --- ---
        # --- initialize Altitude vs Lat Plot ---
        # --- --- --- --- --- --- --- --- --- ---

        # trajectory
        axAlt.plot(Lat[0], Alt[0], color=HighFlyerColor, linewidth=AltvsLat_LineWidth)
        axAlt.plot(Lat[1], Alt[1], color=LowFlyerColor, linewidth=AltvsLat_LineWidth)

        # marker
        Altlat_marker_high = axAlt.scatter(Lat[0][0], Alt[0][0], color=HighFlyerColor, marker='x', linewidth=AltvsLat_scatterLineWidth)
        Altlat_marker_low = axAlt.scatter(Lat[1][0], Alt[1][0], color=LowFlyerColor, marker='x', linewidth=AltvsLat_scatterLineWidth)

        # text
        Altlat_text_high = axAlt.text(Lat[0][0] - AltvsLat_textOffset, Alt[0][0], f'{round(Alt[0][0], AltvsLat_rounding)} km', ha=AltvsLat_text_alignment[1], **AltvsLat_textUTC_style[0])
        Altlat_text_low = axAlt.text(Lat[1][0] - AltvsLat_textOffset, Alt[1][0],f'{round(Alt[1][0], AltvsLat_rounding)} km', ha=AltvsLat_text_alignment[1], **AltvsLat_textUTC_style[1])

        # Labels/Limits
        axAlt.set_ylabel('Altitude [km]', fontsize=AltvsLat_LabelFontsize)
        axAlt.set_xlabel('Geographic Lat', fontsize=AltvsLat_LabelFontsize)
        axAlt.tick_params(axis='x', which='both', labelsize=AltVsLat_xTickSize)
        axAlt.tick_params(axis='y', which='both', labelsize=AltVsLat_yTickSize)
        axAlt.axhline(projectionAltitude[0], color='limegreen', linestyle='--', alpha=0.5, linewidth=AltvsLat_LineWidth)
        axAlt.axhline(projectionAltitude[1], color='limegreen', linestyle='--', alpha=0.5, linewidth=AltvsLat_LineWidth)
        axAlt.text(69.05, 170, '150 km', color='limegreen', fontsize=AltvsLat_textFontSize)
        axAlt.text(69.05, 270, '250 km', color='limegreen', fontsize=AltvsLat_textFontSize)
        axAlt.set_ylim(0, 450)
        axAlt.set_xlim(68, 75)

        # --- --- --- --- --- --- --- --
        # --- INITIALIZE B_GEO FIELD ---
        # --- --- --- --- --- --- --- --
        BprojectionHigh, = axAlt.plot(projectB_lats[0][0], [AltvsLat_Bproject_High, AltvsLat_Bproject_low], color='white', alpha=0.7, linestyle='--')  # high flyer
        BprojectionLow, = axAlt.plot(projectB_lats[1][0], [AltvsLat_Bproject_High, AltvsLat_Bproject_low], color='white', alpha=0.7, linestyle='--', label='$B_{geo}$')  # low flyer
        axAlt.legend(loc='upper left',fontsize=AltvsLat_legendFontsize)

        ######################################
        # --- --- --- --- --- --- --- --- ----
        # --- initialize the All Sky Image ---
        # --- --- --- --- --- --- --- --- ----
        ######################################
        supTitle = fig.suptitle(f'ACESII\n {EpochRocket[0]}\n', fontsize=40,color='white')

        # --- plot the norwegian map data ---
        # 5577A
        ax5577.set_facecolor("gray")
        gl5577 = ax5577.gridlines(draw_labels=True, linewidth=AllSky_gridLineWidth, alpha=0.4, linestyle='--')
        gl5577.xlabel_style = {'size': AllSky_LabelSize, 'color': 'white', 'weight': 'bold'}
        gl5577.ylabel_style = {'size': AllSky_LabelSize, 'color': 'white', 'weight': 'bold'}
        # gl5577.left_labels = False
        gl5577.top_labels = False
        ax5577.set_title('Skibotn 5577$\AA$ - 150 km', fontsize=AllSky_TitleSize,pad=AllSky_TitlePadding)
        ax5577.set_extent([lonW, lonE, latS, latN])  # controls lat/long axes display
        ax5577.coastlines(resolution=res, color='white', alpha=0.8)  # adds coastlines with resolution

        # 6300A
        ax6300.set_facecolor("gray")
        ax6300.set_title('Skibotn 6300$\AA$ - 250 km', fontsize=AllSky_TitleSize,pad=AllSky_TitlePadding)
        gl6300 = ax6300.gridlines(draw_labels=True, linewidth=AllSky_gridLineWidth, color='white', alpha=0.4, linestyle='--')
        gl6300.xlabel_style = {'size': AllSky_LabelSize, 'color': 'white', 'weight': 'bold'}
        gl6300.ylabel_style = {'size': AllSky_LabelSize, 'color': 'white', 'weight': 'bold'}
        gl6300.right_labels = False
        gl6300.top_labels = False
        ax6300.set_extent([lonW, lonE, latS, latN])  # controls lat/long axes display
        ax6300.coastlines(resolution=res, color='white', alpha=0.8)  # adds coastlines with resolution

        # 5577A
        allSky5577 = ax5577.pcolormesh(allGLongs[0], allGlats[0], allImages[0][imageIndicies[0][0]], cmap=allSkyCmap, transform=projTransform, vmin=cbar_Vmin, vmax=cbar_Vmax)
        #
        # # 6300A
        allSky6300 = ax6300.pcolormesh(allGLongs[1], allGlats[1], allImages[1][imageIndicies[0][0]], cmap=allSkyCmap, transform=projTransform, vmin=cbar_Vmin, vmax=cbar_Vmax)


        # --- --- --- --- --- --- --- --- --- --- --- --- --- ----
        # --- initialize total trajectory on the All Sky Image ---
        # --- --- --- --- --- --- --- --- --- --- --- --- --- ----
        # --- GREEN ---
        # plot on 5577 Allsky
        ax5577.plot(Long[1], Lat[1], color=LowFlyerColor, linewidth=AllSky_LineWidth, transform=projTransform)
        ax5577.plot(Long[0], Lat[0], color=HighFlyerColor, linewidth=AllSky_LineWidth, transform=projTransform)
        # ax5577.legend(loc='upper left', prop={'size': AllSky_legendSize})

        # plot marker on 5577 Allsky
        AllSky5577_marker_low = ax5577.scatter(Long[1][0], Lat[1][0], color=LowFlyerColor, marker='x', linewidth=AllSky_markerSize, transform=projTransform)
        AllSky5577_marker_high = ax5577.scatter(Long[0][0], Lat[0][0], color=HighFlyerColor, marker='x', linewidth=AllSky_markerSize, transform=projTransform)

        # --- RED ---
        # plot on 6300 Allsky
        ax6300.plot(Long[1], Lat[1], color=LowFlyerColor,linewidth=4, transform=projTransform)
        ax6300.plot(Long[0], Lat[0], color=HighFlyerColor,linewidth=4, transform=projTransform)

        # plot marker on 6300 Allsky
        AllSky6300_marker_low = ax6300.scatter(Long[1][0], Lat[1][0], color=LowFlyerColor, linewidth=10, marker='x', transform=projTransform)
        AllSky6300_marker_high = ax6300.scatter(Long[0][0], Lat[0][0], color=HighFlyerColor, linewidth=10, marker='x', transform=projTransform)

        # --- Add the cbar ---
        cax = fig.add_axes([0.925, 0.03, 0.02, 0.632])
        cbar = plt.colorbar(mappable=allSky6300, cax=cax, orientation='vertical', fraction=0.046, pad=0.04)
        cbar.set_label('Brightness (kR)', fontsize=cbar_labelFontSize)
        cbar.ax.tick_params(labelsize=cbar_tickFontSize)

        # SAVE FIGURE
        plt.subplots_adjust(left=0.04, bottom=0.03, right=0.92, top=0.9)
        fig.savefig(r'C:\Users\cfelt\Desktop\rockets\ACES-II\Media\movies\AllSky_Movie_initializationFrame.png',dpi=dpi)
        Done(start_time)

        if Animate_AllSkyMovie:
            prgMsg('Creating Allsky Animation')


            ############################
            # --- ANIMATION FUNCTION ---
            ############################
            counter5577 = 0
            counter6300 =0
            def animatePlot(i):

                print('', end='\r' + color.RED + f'{round(100 * i / len(EpochRocket[0]), 1)} %' + color.END)

                # update the title
                supTitle.set_text(f'ACESII\n {EpochRocket[0][i]} UTC\n')

                # --- ALL SKY ---
                # update all sky images
                allSky5577.set_array(allImages[0][imageIndicies[0][i]].ravel())
                allSky6300.set_array(allImages[1][imageIndicies[1][i]].ravel())
                # global counter5577
                # if imageIndicies[0][i] != counter5577:
                #     allSky5577.set_array(allImages[0][imageIndicies[0][i]].ravel())
                #     counter5577 = imageIndicies[0][i]
                # global counter6300
                # if imageIndicies[0][i] != counter6300:
                #
                #     allSky6300.set_array(allImages[1][imageIndicies[1][i]].ravel())
                #     counter6300 = imageIndicies[0][i]

                # update the allsky5577 rocket marker positions
                AllSky5577_marker_high.set_offsets([Long[0][i], Lat[0][i]])
                AllSky5577_marker_low.set_offsets([Long[1][i], Lat[1][i]])

                # update the allsky6300 rocket marker positions
                AllSky6300_marker_high.set_offsets([Long[0][i], Lat[0][i]])
                AllSky6300_marker_low.set_offsets([Long[1][i], Lat[1][i]])

                # --- ALT VS LAT ---
                # update Alt vs Lat marker
                Altlat_marker_high.set_offsets([Lat[0][i], Alt[0][i]])
                Altlat_marker_low.set_offsets([Lat[1][i], Alt[1][i]])

                # update B-project lines on Alt vs lat plot
                BprojectionHigh.set_xdata(projectB_lats[0][i])
                BprojectionLow.set_xdata(projectB_lats[1][i])

                # update Alt vs lat text
                Altlat_text_high.set_x(Lat[0][i])
                Altlat_text_high.set_y(Alt[0][i])
                Altlat_text_high.set_text(f'{round(Alt[0][i], AltvsLat_rounding)} km')

                Altlat_text_low.set_x(Lat[1][i])
                Altlat_text_low.set_y(Alt[1][i])
                Altlat_text_low.set_text(f'{round(Alt[1][i], AltvsLat_rounding)} km')

            ###########################
            # --- ANIMATE THE MOVIE ---
            ###########################
            prgMsg('Creating AllSky Movie')

            if plotSpecificLocations:
                locations = [i for i in specificLocations]
            elif doFrameSkips:
                locations = [i for i in range(0, len(EpochRocket[0]), frame_skips)]  # NEEDS TO BE THE HIGH FLYER LENGTH
            else:
                locations = [i for i in range(len(EpochRocket[0]))]

            anim = animation.FuncAnimation(fig=fig, func=animatePlot, interval=1000 / fps, frames=locations)

            # need a .mp4 writer, following code points to it
            writervideo = animation.FFMpegWriter(fps=fps)

            plt.rcParams['animation.ffmpeg_path'] = r'C:\Users\cfelt\myExecutables\ffmpeg-master-latest-win64-gpl-shared\bin\ffmpeg.exe'
            anim.save(rf'C:\Users\cfelt\Desktop\rockets\ACES-II\Media\movies\ACESII_AllSky_Trajectory.mp4',writer=writervideo)

            Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder

if justPrintSiteNames:
    PlotExtra_AllSkyMovie(justPrintSiteNames,rocketFolderPath)
else:
    PlotExtra_AllSkyMovie(justPrintSiteNames,rocketFolderPath)
