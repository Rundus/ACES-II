# --- Plot1_AllSky.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: All Sky Imager data plot



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
from src.my_imports import *
# 'C:\Users\cfelt\AppData\Local\Microsoft\Windows\Fonts'
# plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import cartopy.crs as ccrs
import spaceToolsLib as spl

print(color.UNDERLINE + f'Plot1_AllSky' + color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# -------------GENERAL PLOT TOGGLES-------------------
# cmapColor = 'viridis'
cmapColor = spl.matlab_parula_cmap()
faceColorChoice = (156 / 255, 156 / 255, 156 / 255, 0.5)  # in normalize RGBA
timeTargetsUTC = [dt.datetime(2022,11,20,17,23,20,100000),
                      dt.datetime(2022,11,20,17,24,00,100000),
                      dt.datetime(2022,11,20,17,24,40,100000),
                      dt.datetime(2022,11,20,17,25,20,100000),
                      dt.datetime(2022,11,20,17,26,00,100000),
                      dt.datetime(2022,11,20,17,26,40,100000),
                      dt.datetime(2022,11,20,17,27,20,100000)] # find the UTC dates times of the specifically sampled labels
# --------------ALTvsLAT------------------
altLatPlot = False
AltLat_Height = 20
AltLat_Width = 35
trajColors = ['tab:red', 'tab:orange']
altLat_labelFontSize = 65
altLat_textSize = 55
altLat_TickLabelSize = 65
altLat_TickLength = 30
altLat_TickWidth = 4
altLat_scatterSize = 800
altLat_lineThickness = 6.5
AltLat_LegendSize = 55
altLat_LabelPadding = 45
# ---------------BigAllSky-----------------
BigAllSkyPlot = False
# lonW = 10
# lonE = 23.5
lonW = 5
lonE = 25
latS = 67
latN = 75
res = '50m'
wImage = 10
cbarVmin,cbarVmax = 0, 12 # in kRayleigh
BigAllSky_textSize = 35
BigAllSky_tickLabelSize = 60
BigAllSky_scatterSize = 550
BigAllSky_lineThickness = 8
BigAllSky_GridSize = 5
BigAllSky_TitleSize = 70
BigAllSky_costLineSize = 3
# --------------------------------
ILatDiffPlot = True
ILatDiff_Height = 15
ILatDiff_Width = 35
ILatDiff_PlotLineWidth = 12
ILatDiff_LabelSize = 75
ILatDiff_TickLabelSize = 60
ILatDiff_TickLength = 25
ILatDiff_TickWidth = 4
ILatDiff_LabelPadding = -15
# --------------------------------
makeColorbarPlot = False


# --- --- --- --- --- --- -
# --- LOAD ALL THE DATA ---
# --- --- --- --- --- --- -

# trajectory
attitudeFolderPath = f'{ACES_data_folder}\\attitude'
inputFilesTraj = [glob(f'{attitudeFolderPath}\\{fliers[0]}\\*.cdf*')[0],
                  glob(f'{attitudeFolderPath}\\{fliers[1]}\\*.cdf*')[0]]

# --- GET TRAJECTORY DATA ---
prgMsg(f'Loading ACESII traj data')
data_dicts_attitude = [loadDictFromFile(inputFilesTraj[0]),loadDictFromFile(inputFilesTraj[1])]

# define some variables
EpochRocket = [data_dicts_attitude[0]['Epoch'][0], data_dicts_attitude[1]['Epoch'][0]]
geoAlt = [data_dicts_attitude[0]['Alt'][0], data_dicts_attitude[1]['Alt'][0]]
geoLat = [data_dicts_attitude[0]['Lat'][0], data_dicts_attitude[1]['Lat'][0]]
geoLong = [data_dicts_attitude[0]['Long'][0], data_dicts_attitude[1]['Long'][0]]
geoMagLat = [data_dicts_attitude[0]['Lat_geom'][0], data_dicts_attitude[1]['Lat_geom'][0]]
Done(start_time)

# Load AllSky data
prgMsg('Loading Allsky Data')
data_dict_allSky5577 = loadDictFromFile(glob(r'C:\Data\ACESII\all_sky\skibotn\5577\\*.cdf')[0])
data_dict_allSky6300 = loadDictFromFile(glob(r'C:\Data\ACESII\all_sky\skibotn\6300\\*.cdf')[0])
Done(start_time)


############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################


# --- ALTITUDE VS LATITUDE PLOT ---
if altLatPlot:
    prgMsg('Plotting AltLat')

    # --- --- --- --- ---
    # --- AltLat plot ---
    # --- --- --- --- ---
    # axAltLat = fig.add_subplot(gs_altLat_BigAllSky[0])
    fig, axAltLat = plt.subplots()
    figure_height = AltLat_Height
    figure_width = AltLat_Width
    fig.set_figwidth(figure_width)
    fig.set_figheight(figure_height)

    axAltLat.set_ylabel('Altitude [km]', fontsize=altLat_labelFontSize,weight='bold',labelpad=altLat_LabelPadding)
    axAltLat.set_xlabel('Geographic Lat [deg]',fontsize=altLat_labelFontSize,weight='bold',labelpad=altLat_LabelPadding-40)
    axAltLat.set_xlim(69, 75)
    axAltLat.set_ylim(0, 500)

    # plot the pseudo geomagnetic field line
    slope = -1 * (111 / np.sin(np.radians(90 - 78.13)))  # corresponds to line with -78.13deg inclination
    for i in range(31):
        axAltLat.axline(xy1=(69 + i * 0.25, 0), slope=slope, color='tab:blue', linewidth=altLat_lineThickness, linestyle='-.', alpha=0.3)
    axAltLat.legend(['B$_{Geo}$'], loc='upper right',fontsize=AltLat_LegendSize)

    # set the facecolor of the axAltLat plot
    axAltLat.set_facecolor(faceColorChoice)

    # plot the UTC labels
    axGeomLat = axAltLat.twiny()
    axGeomLat.plot(geoMagLat[0], geoAlt[0]/1000, color=trajColors[0], alpha=0)  # High
    axGeomLat.plot(geoMagLat[1], geoAlt[1]/1000, color=trajColors[1], alpha=0)  # Low
    axGeomLat.set_xlabel('Geomagnetic Latitude [deg]',fontsize=altLat_labelFontSize,weight='bold',labelpad=altLat_LabelPadding)
    AltLat_vertical_Alignments = ['bottom' for tme in timeTargetsUTC]
    AltLat_horizontal_Alignments = ['right', 'right', 'right', 'center', 'left', 'left', 'left']
    # vertical_text_label_adjustments = [-0.04, -0.03, 0.02, 0.06, 0.02, -0.04, -0.04]
    vertical_text_label_adjustments = [-0.09, -0.06, -0.005, 0.04, -0.01, -0.06, -0.09]
    horizontal_text_label_adjustments = [-0.002, -0.0015, -0.001, 0.0, 0.001, 0.0015, 0.002]

    # plot the scatterpoint of each of the timeTargetUTC_labels. Plot the text itself only for the High Flyer and
    # create a connecting line between the scatterpoints between the flyers
    for i in range(2):  # for each rocket
        for j, ttme in enumerate(timeTargetsUTC):
            Index = np.abs(EpochRocket[i] - ttme).argmin()
            xPos = geoLat[i][Index]
            yPos = geoAlt[i][Index]/1000

            if i == 0:
                # plot the text itself
                label = EpochRocket[i][Index]
                deltaY = vertical_text_label_adjustments[j] * yPos
                deltaX = horizontal_text_label_adjustments[j] * xPos
                axAltLat.text(x=xPos + deltaX, y=yPos + deltaY, s=label.strftime("%H:%M:%S"), color='black',weight='bold',
                              va=AltLat_vertical_Alignments[j], ha=AltLat_horizontal_Alignments[j], size=altLat_textSize)

                # plot the connecting line
                Index_LF = np.abs(EpochRocket[1] - ttme).argmin()
                xPos_LF = geoLat[1][Index_LF]
                yPos_LF = geoAlt[1][Index_LF]/1000
                axAltLat.plot([xPos, xPos_LF], [yPos, yPos_LF], color='green', linestyle='--', alpha=0.5,linewidth=altLat_lineThickness)
                # axAltLat.axline(xy1=(xPos, yPos), xy2=(xPos_LF, yPos_LF), color='green', linestyle='--', alpha=0.5)

            # plot a dot at the text label
            axAltLat.scatter(x=xPos, y=yPos, s=altLat_scatterSize, marker="o", color=trajColors[i])

    # adjust the tick label size
    axAltLat.tick_params(axis='both', labelsize=altLat_TickLabelSize, length=altLat_TickLength, width=altLat_TickWidth)
    axAltLat.tick_params(axis='both', which='minor', length=int(altLat_TickLength*0.65), width=altLat_TickWidth)
    axGeomLat.tick_params(axis='both', labelsize=altLat_TickLabelSize, length=altLat_TickLength, width=altLat_TickWidth)
    axGeomLat.tick_params(axis='both', which='minor', length=int(altLat_TickLength*0.65) , width=altLat_TickWidth)
    axAltLat.minorticks_on()
    axGeomLat.minorticks_on()

    # plot the trajectory over everything
    axAltLat.plot(geoLat[0], geoAlt[0]/1000, color=trajColors[0], label='High Flyer',linewidth=altLat_lineThickness)  # High
    axAltLat.plot(geoLat[1], geoAlt[1]/1000, color=trajColors[1], label='Low Flyer',linewidth=altLat_lineThickness)  # Low

    plt.tight_layout()
    plt.savefig(r'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot1\\AltLat.png')
    Done(start_time)
    # plt.show()

# --- BIG ALLSKYIMAGER PLOT ---
if BigAllSkyPlot:
    prgMsg('Plotting BigAllSky')

    BigAllSky_wScatterPoints = np.array([int((1E-9) * (pycdf.lib.datetime_to_tt2000(tme) - pycdf.lib.datetime_to_tt2000( dt.datetime(2022, 11, 20, 17, 20, 00, 000)))) for tme in timeTargetsUTC])

    # --- --- --- --- --- --
    # --- BigAllSky plot ---
    # --- --- --- --- --- --
    projProjection = ccrs.Orthographic(central_longitude=15, central_latitude=70)
    projTransform = ccrs.PlateCarree()

    for i in range(2):
        # --- PLOT MAP OF NORWAY ---
        fig, axBigAllSky = plt.subplots(1,subplot_kw=dict(projection=projProjection))
        figure_height = 20
        figure_width = 20
        fig.set_figwidth(figure_width)
        fig.set_figheight(figure_height)

        # get the elevation map data and plot it
        # fname = glob(r"C:\Users\cfelt\PycharmProjects\UIOWA_CDF_operator\ACESII_code\Papers\ACESII_Alfvenic_Observations\ElevationData\*.nc*")[0]
        # dataset = netcdf_dataset(fname) # Load data into separate arrays,if variables are not known, print(dataset.variables) to check them
        # elev = dataset.variables['elevation'][:]
        # lats = dataset.variables['lat'][:]
        # lons = dataset.variables['lon'][:]
        # vmin, vmax = -8000, 3000
        # v = np.linspace(vmin, vmax, 100, endpoint=True) # Set how many contour lines to display
        # axBigAllSky.contourf(lons, lats, elev, v, cmap="gray", vmin=vmin, vmax=vmax, transform=projTransform)

        # gridlines
        gl = axBigAllSky.gridlines(draw_labels=True, linewidth=BigAllSky_GridSize,
                                   alpha=0.4,
                                   linestyle='--',
                                   color='black')
        gl.xlabel_style = {'size': BigAllSky_tickLabelSize, 'color': 'black', 'weight': 'bold'}
        gl.ylabel_style = {'size': BigAllSky_tickLabelSize, 'color': 'black', 'weight': 'bold'}
        gl.top_labels = False

        # extent of map
        axBigAllSky.set_extent([lonW, lonE, latS, latN])  # controls lat/long axes display

        # coastlines
        axBigAllSky.coastlines(resolution=res, color='black',  alpha=1,linewidth=BigAllSky_costLineSize)  # adds coastlines with resolution

        if i == 0:
            #--- Plot the Big AllSky image ---
            cmapBigAllSky = axBigAllSky.pcolormesh(data_dict_allSky5577['GLongs'][0], data_dict_allSky5577['GLats'][0], data_dict_allSky5577['AllSkyImages'][0][wImage],
                                                   cmap=cmapColor,
                                                   transform=projTransform,
                                                   vmin=cbarVmin,
                                                   vmax=cbarVmax,
                                                   alpha=1)
            BigAllSky_outputPath = r'C:\Users\cfelt\Desktop\rockets\ACES-II\Papers\ACESII_Alfven_Observations\Plot1\\BigAllSky_5570.png'
            fig.suptitle('Skibotn 5577$\AA$ - 150 km\n' + data_dict_allSky5577['Epoch'][0][wImage].strftime("%Y-%B-%d %H:%M:%S") + ' UTC',fontsize=BigAllSky_TitleSize,weight='bold')
        elif i == 1:
            cmapBigAllSky = axBigAllSky.pcolormesh(data_dict_allSky6300['GLongs'][0], data_dict_allSky6300['GLats'][0],
                                                   data_dict_allSky6300['AllSkyImages'][0][wImage],
                                                   cmap=cmapColor,
                                                   transform=projTransform,
                                                   vmin=cbarVmin,
                                                   vmax=cbarVmax,
                                                   alpha=1)
            BigAllSky_outputPath = r'C:\Users\cfelt\Desktop\rockets\ACES-II\Papers\ACESII_Alfven_Observations\Plot1\\BigAllSky_6300.png'
            fig.suptitle('Skibotn 6300$\AA$ - 250 km\n' + data_dict_allSky6300['Epoch'][0][wImage].strftime("%Y-%B-%d %H:%M:%S")+ ' UTC',fontsize=BigAllSky_TitleSize,weight='bold')
        axBigAllSky.set_facecolor(faceColorChoice)

        # --- plot the rocket trajectory data on the large AllSky plot ---
        axBigAllSky.plot(geoLong[0], geoLat[0], color=trajColors[0], transform=projTransform,linewidth=BigAllSky_lineThickness) # High
        axBigAllSky.plot(geoLong[1], geoLat[1], color=trajColors[1], transform=projTransform,linewidth=BigAllSky_lineThickness) # Low

        # plot specific UTC times on the trajectory lines
        timeTargetsUTC_labels = [pycdf.lib.tt2000_to_datetime(int(tme*1E9 + pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 20, 00, 000)))) for tme in BigAllSky_wScatterPoints]

        alignment = ['left', 'right']
        for i in range(2): # for each rocket
            for tt,ttme in enumerate(timeTargetsUTC_labels):

                Index = np.abs(EpochRocket[i] - ttme).argmin()
                label = EpochRocket[i][Index]

                xPos = geoLong[i][Index]
                yPos = geoLat[i][Index]

                # plot a dot at the text label
                axBigAllSky.scatter(x=xPos,y=yPos,s=BigAllSky_scatterSize,marker="o",color=trajColors[i],transform=projTransform)

                # plot the text itself
                deltaX = 0.0075*xPos if i == 0 else -1*0.0075*xPos
                # if tt==0 or tt==6:
                #     axBigAllSky.text(x=xPos + deltaX, y=yPos, s=label.strftime("%H:%M:%S"), color='white', weight='bold',
                #                      ha=alignment[i], transform=projTransform, size=BigAllSky_textSize+4)

        plt.tight_layout()
        plt.savefig(BigAllSky_outputPath)
        # plt.show()
        Done(start_time)

if ILatDiffPlot:
    prgMsg('Plotting ILatDiff')

    ttIndicies = [np.abs(data_dicts_attitude[0]['Epoch'][0] - timeTargetsUTC[0]).argmin(),np.abs(data_dicts_attitude[0]['Epoch'][0] - timeTargetsUTC[-1]).argmin()]
    ILat = data_dicts_attitude[0]['ILat'][0][ttIndicies[0]:ttIndicies[1]]
    ILong = data_dicts_attitude[0]['ILong'][0][ttIndicies[0]:ttIndicies[1]]
    Epoch = data_dicts_attitude[0]['Epoch'][0][ttIndicies[0]:ttIndicies[1]]
    spatialDiff_ILat = data_dicts_attitude[0]['footPrint_ILat_km_Diff'][0][ttIndicies[0]:ttIndicies[1]]
    timeDiff_Ilat = data_dicts_attitude[0]['footPrint_lattitude_time_Difference'][0][ttIndicies[0]:ttIndicies[1]]
    spatialDiff_ILong = data_dicts_attitude[0]['footPrint_ILong_km_Diff'][0][ttIndicies[0]:ttIndicies[1]]
    timeDiff_Ilong = data_dicts_attitude[0]['footPrint_longitude_time_Difference'][0][ttIndicies[0]:ttIndicies[1]]

    ###########################
    # --- --- --- --- --- --- -
    # --- ILatDiffPlot plot ---
    # --- --- --- --- --- --- -
    ###########################
    fig, (axILatDiff_time,axILatDiff_space) = plt.subplots(2, sharex=True)
    figure_height = ILatDiff_Height
    figure_width = ILatDiff_Width
    fig.set_figwidth(figure_width)
    fig.set_figheight(figure_height)

    # --- SPATIAL---
    axILatDiff_space.plot(Epoch, spatialDiff_ILat, color='black', linewidth=ILatDiff_PlotLineWidth)
    axILatDiff_space.set_ylabel(r'ILat $\Delta \varphi$'+'\n[km]', fontsize=ILatDiff_LabelSize, labelpad=-75, weight='bold')
    axILatDiff_space.tick_params(axis='both', labelsize=ILatDiff_TickLabelSize, length=ILatDiff_TickLength, width=ILatDiff_TickWidth)
    axILatDiff_space.tick_params(axis='both', which='minor', length=int(ILatDiff_TickLength * 0.65), width=ILatDiff_TickWidth)
    axILatDiff_space.margins(0)
    axILatDiff_space.grid(which='both', linewidth=2.5, color='gray', alpha=0.6)
    axILatDiff_space.set_ylim(-110, 110)
    axILatDiff_space.axhline(y=0, color='red', linestyle='--', linewidth=ILatDiff_PlotLineWidth-5, alpha=0.4)
    axILatDiff_space.minorticks_on()
    axILatDiff_space.set_xlabel('Time [UTC]',fontsize=ILatDiff_LabelSize,labelpad=15, weight='bold')

    # extraAxis = axILatDiff_space.twiny()
    # extraAxis.plot(Epoch, spatialDiff_ILat, alpha=0)
    # extraAxis.set_xticks(timeTargetsUTC)
    # tickIndicies = [np.abs(Epoch - val).argmin() for val in timeTargetsUTC]
    # extraAxis.set_xticklabels([round(ILat[idx], 2) for idx in tickIndicies])
    # extraAxis.set_xlabel('Low Flyer ILat [deg]', fontsize=ILatDiff_LabelSize, labelpad=ILatDiff_LabelPadding+50)
    # extraAxis.tick_params(axis='both', labelsize=ILatDiff_TickLabelSize, length=ILatDiff_TickLength, width=ILatDiff_TickWidth)
    # extraAxis.tick_params(axis='both', which='minor', length=int(ILatDiff_TickLength * 0.65), width=ILatDiff_TickWidth)
    # extraAxis.margins(x=0)
    # extraAxis.set_xlim(timeTargetsUTC[0],timeTargetsUTC[-1])
    # set the ticks

    # --- TEMPORAL ---
    axILatDiff_time.plot(Epoch, timeDiff_Ilat, color='black', linewidth=ILatDiff_PlotLineWidth)
    axILatDiff_time.grid(which='both', linewidth=2.5, color='gray', alpha=0.6)
    axILatDiff_time.set_ylabel('ILat $\Delta t$ \n [s]', fontsize=ILatDiff_LabelSize, labelpad=-45, weight='bold')
    axILatDiff_time.set_ylim(-65, 65)
    axILatDiff_time.margins(0)

    # TICKS
    axILatDiff_time.set_xticks(timeTargetsUTC)
    tlabels = [label.strftime("%H:%M:%S") for label in timeTargetsUTC]
    axILatDiff_time.set_xticklabels(tlabels)
    axILatDiff_time.tick_params(axis='both', labelsize=ILatDiff_TickLabelSize, length=ILatDiff_TickLength, width=ILatDiff_TickWidth)
    axILatDiff_time.tick_params(axis='both', which='minor', length=int(ILatDiff_TickLength * 0.65), width=ILatDiff_TickWidth)
    axILatDiff_time.axhline(y=0, color='red', linestyle='--', linewidth=ILatDiff_PlotLineWidth - 5, alpha=0.4)
    axILatDiff_time.minorticks_on()
    # fig.align_ylabels([axILatDiff_space, axILatDiff_time])
    plt.tight_layout()
    plt.savefig(r'C:\Users\cfelt\Desktop\rockets\ACES-II\Papers\ACESII_Alfven_Observations\PLOTS\Plot1\ILatDiff.png')
    Done(start_time)

    ###########################
    # --- --- --- --- --- --- -
    # --- ILongDiffPlot plot ---
    # --- --- --- --- --- --- -
    ###########################
    # fig, (axILongDiff_space, axIongDiff_time) = plt.subplots(2, sharex=True)
    # figure_height = ILatDiff_Height
    # figure_width = ILatDiff_Width
    # fig.set_figwidth(figure_width)
    # fig.set_figheight(figure_height)
    #
    # # --- --- --- --- ---
    # # --- ILONG PLOTS ---
    # # --- --- --- --- ---
    #
    # # --- SPATIAL ---
    # axILongDiff_space.plot(Epoch, spatialDiff_ILong, color='tab:green', linewidth=ILatDiff_PlotLineWidth)
    # axILongDiff_space.set_ylabel(r'ILong $\Delta \lambda$' + '\n[km]', fontsize=ILatDiff_LabelSize, labelpad=ILatDiff_LabelPadding, weight='bold')
    # axILongDiff_space.tick_params(axis='both', labelsize=ILatDiff_TickLabelSize, length=ILatDiff_TickLength, width=ILatDiff_TickWidth)
    # axILongDiff_space.tick_params(axis='both', which='minor', length=int(ILatDiff_TickLength * 0.65), width=ILatDiff_TickWidth)
    # axILongDiff_space.margins(0)
    # axILongDiff_space.grid(which='both', linewidth=2.5, color='gray', alpha=0.6)
    # axILongDiff_space.set_ylim(-60, 60)
    # axILongDiff_space.axhline(y=0, color='red', linestyle='--', linewidth=ILatDiff_PlotLineWidth - 5, alpha=0.4)
    # axILongDiff_space.minorticks_on()
    #
    # extraAxis = axILongDiff_space.twiny()
    # extraAxis.plot(Epoch, spatialDiff_ILong, alpha=0)
    # extraAxis.set_xticks(timeTargetsUTC)
    # tickIndicies = [np.abs(Epoch - val).argmin() for val in timeTargetsUTC]
    # extraAxis.set_xticklabels([round(ILong[idx], 2) for idx in tickIndicies])
    # extraAxis.set_xlabel('Low Flyer ILong [deg]', fontsize=ILatDiff_LabelSize, labelpad=ILatDiff_LabelPadding + 50)
    # extraAxis.tick_params(axis='both', labelsize=ILatDiff_TickLabelSize, length=ILatDiff_TickLength, width=ILatDiff_TickWidth)
    # extraAxis.tick_params(axis='both', which='minor', length=int(ILatDiff_TickLength * 0.65), width=ILatDiff_TickWidth)
    # extraAxis.margins(x=0)
    # extraAxis.set_xlim(timeTargetsUTC[0], timeTargetsUTC[-1])
    #
    # # --- TEMPORAL ---
    # axIongDiff_time.plot(Epoch, timeDiff_Ilong, color='black', linewidth=ILatDiff_PlotLineWidth)
    # axIongDiff_time.grid(which='both', linewidth=2.5, color='gray', alpha=0.6)
    # axIongDiff_time.set_ylabel('ILong $\Delta t$ \n [s]', fontsize=ILatDiff_LabelSize, labelpad=ILatDiff_LabelPadding, weight='bold')
    # axIongDiff_time.set_ylim(40, 250)
    # axIongDiff_time.margins(0)
    # axIongDiff_time.set_xlabel('time (UTC)', fontsize=ILatDiff_LabelSize, labelpad=ILatDiff_LabelPadding + 50)
    #
    # # TICKS
    # axIongDiff_time.set_xticks(timeTargetsUTC)
    # tlabels = [label.strftime("%H:%M:%S") for label in timeTargetsUTC]
    # axIongDiff_time.set_xticklabels(tlabels)
    # axIongDiff_time.tick_params(axis='both', labelsize=ILatDiff_TickLabelSize, length=ILatDiff_TickLength, width=ILatDiff_TickWidth)
    # axIongDiff_time.tick_params(axis='both', which='minor', length=int(ILatDiff_TickLength * 0.65), width=ILatDiff_TickWidth)
    # axIongDiff_time.axhline(y=0, color='red', linestyle='--', linewidth=ILatDiff_PlotLineWidth - 5, alpha=0.4)
    # axIongDiff_time.minorticks_on()
    # # axILatDiff_time.margins(0)
    #
    #
    # fig.align_ylabels([axILongDiff_space,axIongDiff_time])
    # plt.tight_layout()
    # plt.savefig(r'C:\Users\cfelt\Desktop\rockets\ACES-II\Papers\ACESII_Alfven_Observations\PLOTS\Plot1\ILongDiff.png')
    # Done(start_time)


# --- --- --- ----
# --- COLORBAR ---
# --- --- --- ----
if makeColorbarPlot:

    # --- PLOT MAP OF NORWAY ---
    projTransform = ccrs.PlateCarree()
    fig, axColorbarPlot = plt.subplots(1, subplot_kw=dict(projection=projTransform))
    figure_height = 8
    figure_width = 40
    fig.set_figwidth(figure_width)
    fig.set_figheight(figure_height)

    # gridlines
    gl = axColorbarPlot.gridlines(draw_labels=True, linewidth=BigAllSky_GridSize,
                               alpha=0.0,
                               linestyle='--',
                               color='black')
    gl.xlabel_style = {'size': BigAllSky_tickLabelSize, 'color': 'black', 'weight': 'bold'}
    gl.ylabel_style = {'size': BigAllSky_tickLabelSize, 'color': 'black', 'weight': 'bold'}
    gl.top_labels = False
    gl.bottom_labels = False

    # extent of map
    axColorbarPlot.set_extent([lonW, lonE, latS, latN])  # controls lat/long axes display

    # --- Plot the Big AllSky image ---
    cmapBigAllSky = axColorbarPlot.pcolormesh(data_dict_allSky5577['GLongs'][0],
                                              data_dict_allSky5577['GLats'][0],
                                           data_dict_allSky5577['AllSkyImages'][0][wImage],
                                           cmap=cmapColor,
                                           transform=projTransform,
                                           vmin=cbarVmin,
                                           vmax=cbarVmax,
                                           alpha=1)

    cbar = plt.colorbar(mappable=cmapBigAllSky, orientation='horizontal',fraction=0.8, pad=2)
    cbar.set_label('Intensity [kR]', fontsize=65, weight='bold')
    cbar.ax.tick_params(labelsize=60)

    # plt.tight_layout()
    plt.savefig(r'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot1\colorbar.png')
    # plt.show()
    Done(start_time)
