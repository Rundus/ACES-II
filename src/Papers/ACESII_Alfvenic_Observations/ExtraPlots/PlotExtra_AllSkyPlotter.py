# --- Plot1_AllSky.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: All Sky Imager data plot



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
from src.myImports import *
# 'C:\Users\cfelt\AppData\Local\Microsoft\Windows\Fonts'
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import cartopy.crs as ccrs
from netCDF4 import Dataset as netcdf_dataset
from spaceToolsLib.colorbars.matlab_parula import matlab_parula_cmap

print(color.UNDERLINE + f'Plot1_AllSky' + color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# -------------GENERAL PLOT TOGGLES-------------------
# cmapColor = 'viridis'
cmapColor = matlab_parula_cmap()
faceColorChoice = (156 / 255, 156 / 255, 156 / 255, 0.5)  # in normalize RGBA
timeTargetsUTC = [dt.datetime(2022,11,20,17,23,20,100000),
                      dt.datetime(2022,11,20,17,24,00,100000),
                      dt.datetime(2022,11,20,17,24,40,100000),
                      dt.datetime(2022,11,20,17,25,20,100000),
                      dt.datetime(2022,11,20,17,26,00,100000),
                      dt.datetime(2022,11,20,17,26,40,100000),
                      dt.datetime(2022,11,20,17,27,20,100000)] # find the UTC dates times of the specifically sampled labels
# ---------------BigAllSky-----------------
showNorway = True
showTrajectories = False
grayBackground = False
wImage = []
trajColors = ['tab:red', 'tab:orange']
lonW = 7.5
lonE = 27.5
latS = 68.5
latN = 73
# lonW = 0
# lonE = 180
# latS = 0
# latN = 180
res = '50m' # 50m or 10m
cbarVmin, cbarVmax = 0, 10 # in kRayleigh
BigAllSky_textSize = 35
BigAllSky_tickLabelSize = 0
BigAllSky_scatterSize = 550
BigAllSky_lineThickness = 8
BigAllSky_GridSize = 5
BigAllSky_TitleSize = 90
BigAllSky_costLineSize = 3
# --------------------------------
makeColorbarPlot = True


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
if wImage == []:
    wImage = [i for i in range(len(data_dict_allSky5577['AllSkyImages'][0]))]


for ImageIdx in wImage:
    # --- BIG ALLSKYIMAGER PLOT ---

    prgMsg('Plotting BigAllSky')

    BigAllSky_wScatterPoints = np.array([int((1E-9) * (pycdf.lib.datetime_to_tt2000(tme) - pycdf.lib.datetime_to_tt2000( dt.datetime(2022, 11, 20, 17, 20, 00, 000)))) for tme in timeTargetsUTC])

    # --- --- --- --- --- --
    # --- BigAllSky plot ---
    # --- --- --- --- --- --
    for i in range(2):
        # --- PLOT MAP OF NORWAY ---
        projProjection =ccrs.Orthographic(central_longitude=15, central_latitude=70)
        projTransform = ccrs.PlateCarree()

        fig, axBigAllSky = plt.subplots(1,subplot_kw=dict(projection=projProjection))
        figure_height = 20
        figure_width = 20
        fig.set_figwidth(figure_width)
        fig.set_figheight(figure_height)

        # get the elevation map data and plot it
        if grayBackground:
            fname = glob(r"C:\Users\cfelt\PycharmProjects\UIOWA_CDF_operator\ACESII_code\Papers\ACESII_Alfvenic_Observations\ElevationData\*.nc*")[0]
            dataset = netcdf_dataset(fname) # Load data into separate arrays,if variables are not known, print(dataset.variables) to check them
            elev = dataset.variables['elevation'][:]
            lats = dataset.variables['lat'][:]
            lons = dataset.variables['lon'][:]
            vmin, vmax = -8000, 3000
            v = np.linspace(vmin, vmax, 100, endpoint=True) # Set how many contour lines to display
            axBigAllSky.contourf(lons, lats, elev, v, cmap="gray", vmin=vmin, vmax=vmax, transform=projTransform)

        # gridlines
        gl = axBigAllSky.gridlines(draw_labels=True, linewidth=BigAllSky_GridSize,
                                   alpha=0.25,
                                   linestyle='--',
                                   color='black', zorder=2)

        gl.xlabel_style = {'size': BigAllSky_tickLabelSize, 'color': 'black', 'weight': 'bold'}
        gl.ylabel_style = {'size': BigAllSky_tickLabelSize, 'color': 'black', 'weight': 'bold'}
        gl.top_labels = False

        # extent of map
        axBigAllSky.set_extent([lonW, lonE, latS, latN])  # controls lat/long axes display

        # coastlines
        if showNorway:
            axBigAllSky.coastlines(resolution=res, color='black',  alpha=1,linewidth=BigAllSky_costLineSize,zorder=2)  # adds coastlines with resolution



        if i == 0:
            #--- Plot the Big AllSky image ---

            cmapBigAllSky = axBigAllSky.pcolormesh(data_dict_allSky5577['GLongs'][0], data_dict_allSky5577['GLats'][0], data_dict_allSky5577['AllSkyImages'][0][ImageIdx],
                                                   cmap=cmapColor,
                                                   transform=projTransform,
                                                   vmin=cbarVmin,
                                                   vmax=cbarVmax,
                                                   alpha=1,
                                                   zorder=1)
            axBigAllSky.text(0.5, 1.05, s=data_dict_allSky5577['Epoch'][0][ImageIdx].strftime("%H:%M:%S"), fontsize=BigAllSky_TitleSize, ha='center',transform=axBigAllSky.transAxes,weight='bold')
            BigAllSky_outputPath = rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\PlotExtra\AllSky\\BigAllSky_5570_{ImageIdx}.png'
            # fig.suptitle('Skibton 5577$\AA$ - 150 km\n' + data_dict_allSky5577['Epoch'][0][ImageIdx].strftime("%Y-%B-%d %H:%M:%S") + ' UTC',fontsize=BigAllSky_TitleSize,weight='bold')
        elif i == 1:
            cmapBigAllSky = axBigAllSky.pcolormesh(data_dict_allSky6300['GLongs'][0], data_dict_allSky6300['GLats'][0], data_dict_allSky6300['AllSkyImages'][0][ImageIdx],
                                                   cmap=cmapColor,
                                                   transform=projTransform,
                                                   vmin=cbarVmin,
                                                   vmax=cbarVmax,
                                                   alpha=1,
                                                   zorder=1)
            axBigAllSky.text(0.5,1.05,s=data_dict_allSky6300['Epoch'][0][ImageIdx].strftime("%H:%M:%S"),fontsize=BigAllSky_TitleSize, ha='center',transform=axBigAllSky.transAxes,weight='bold')
            BigAllSky_outputPath = rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\PlotExtra\AllSky\\BigAllSky_6300_{ImageIdx}.png'
            # fig.suptitle('Skibton 6300$\AA$ - 250 km\n' + data_dict_allSky6300['Epoch'][0][ImageIdx].strftime("%Y-%B-%d %H:%M:%S")+ ' UTC',fontsize=BigAllSky_TitleSize,weight='bold')
        axBigAllSky.set_facecolor(faceColorChoice)


        if showTrajectories:
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
                    axBigAllSky.scatter(x=xPos,y=yPos,s=BigAllSky_scatterSize,marker="o",color=trajColors[i],transform=projTransform,zorder=2)

                    # plot the text itself
                    deltaX = 0.0075*xPos if i == 0 else -1*0.0075*xPos
                    # if tt==0 or tt==6:
                    #     axBigAllSky.text(x=xPos + deltaX, y=yPos, s=label.strftime("%H:%M:%S"), color='white', weight='bold',
                    #                      ha=alignment[i], transform=projTransform, size=BigAllSky_textSize+4)

        plt.tight_layout()
        plt.savefig(BigAllSky_outputPath)
        plt.close()
        Done(start_time)






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
    cmapBigAllSky = axColorbarPlot.pcolormesh(data_dict_allSky5577['GLongs'][0], data_dict_allSky5577['GLats'][0],
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
    plt.savefig(r'C:\Users\cfelt\OneDrive\Desktop\Paper_Photos\Plot1\colorbar.png')
    # plt.show()
    Done(start_time)
