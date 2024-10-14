# --- ESAmovie.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Produce a movie of ESA data in either (1) DiffEFlux (2) DiffNFlux or (3) Distribution Function
# for the EEPAA, IEPAA or LEESA (Depending on the rocket).


# --- --- --- --- ---
from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
frame_skips = 1 # 1 is no skip i.e. all frames used, 5 is use only every 5th frame, 10 etc...
wInstr_movie = 'eepaa'

# plot specifc locations
plotSpecificLocations = False
specific_locations = [6000 + i for i in range(0, 100)]

# --- reduce dataset ---
plotSpecificTimeRange = True
targetTimes = [pycdf.lib.datetime_to_tt2000(dt.datetime(2022,11,20,17,24,54,000)), pycdf.lib.datetime_to_tt2000(dt.datetime(2022,11,20,17,25,11,000))]

# -- parameters involved in plotting the interpolated distribution function/DiffNFlux --
wDataPlot = [2] # 0 - Counts, 1 - DiffEFlux, 2 - Dist. Func
N = 200 # number of points for NxN interpolation grid
fillvalue = 0 # interpolated plots fillvalue

normalizeToThermal = True # normalize ESA plots to the electron thermal velocity, Te = 0.1


# --- --- --- ---
# --- import ---
# --- --- --- ---

from matplotlib import pyplot as plt, animation
from matplotlib.patches import Wedge
from scipy.interpolate import LinearNDInterpolator
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from ACESII_code.class_var_func import color, loadDictFromFile

print(color.BOLD + color.CYAN + 'ACESII_trajectory_plotting.py' + color.END + color.END)


def storeData(inputDict, dataPaths):
    parentDict = deepcopy(inputDict)
    for pkey, pval in parentDict.items():  # for each instrument, 3 loops

        size = 1 if pkey == 'leesa' else 2  # adjust for leesa only being on High Flyer

        for i in range(size):  # one for each flyer
            data_dict = {}

            try: # Case for if the datafile doesn't exist yet
                with pycdf.CDF(dataPaths[pkey][i][0]) as dataFile:
                    for key, val in dataFile.items():
                        data_dict = {**data_dict, **{key: [dataFile[key][...], {key: val for key, val in dataFile[key].attrs.items()}]}}
                try:
                    data_dict['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_esa'][0][i]) for i in range(len(data_dict['Epoch_esa'][0]))])
                except:
                    data_dict['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch'][0][i]) for i in range(len(data_dict['Epoch'][0]))])

                parentDict[pkey].append(data_dict)
            except:
                parentDict[pkey].append(data_dict)

    return parentDict

def ACESIIplotting():

    rocketAttrs, b, c = ACES_mission_dicts()


    # Trajectory Data
    trajFolderPath = f'{ACES_data_folder}trajectories\\'
    dataPath_traj = [glob(trajFolderPath + rf'{fliers[0]}\\*_ILat_ILong*'),
                     glob(trajFolderPath + rf'{fliers[1]}\\\\*_ILat_ILong*')]

    # L1 ESA data
    FolderPath = f'{ACES_data_folder}L1\\'
    dataPath_L1ESA = {'eepaa':[glob(FolderPath + rf'{fliers[0]}\\*eepaa_*'), glob(FolderPath + rf'{fliers[1]}\\\\*eepaa_*')],
                      'iepaa':[glob(FolderPath + rf'{fliers[0]}\\*iepaa*'),glob(FolderPath + rf'{fliers[1]}\\\\*iepaa*')],
                      'leesa':[glob(FolderPath + rf'{fliers[0]}\\*leesa*'),glob(FolderPath + rf'{fliers[1]}\\\\*leesa*')]}
    # L2 ESA data
    FolderPath = f'{ACES_data_folder}L2\\'
    dataPath_L2ESA = {'eepaa': [glob(FolderPath + rf'{fliers[0]}\\*eepaa_fullCal*'), glob(FolderPath + rf'{fliers[1]}\\\\*eepaa_fullCal*')],
                      'iepaa': [glob(FolderPath + rf'{fliers[0]}\\*iepaa*'), glob(FolderPath + rf'{fliers[1]}\\\\*iepaa*')],
                      'leesa': [glob(FolderPath + rf'{fliers[0]}\\*leesa*'), glob(FolderPath + rf'{fliers[1]}\\\\*leesa*')]}

    # Distribution Function Data
    FolderPath = f'{ACES_data_folder}\science\DistFunc\\'
    dataPath_Dist = {'eepaa':[glob(FolderPath + rf'{fliers[0]}\\*eepaa_fullCal*'), glob(FolderPath + rf'{fliers[1]}\\\\*eepaa_fullCal*')],
                      'iepaa':[glob(FolderPath + rf'{fliers[0]}\\*iepaa*'),glob(FolderPath + rf'{fliers[1]}\\\\*iepaa*')],
                      'leesa':[glob(FolderPath + rf'{fliers[0]}\\*leesa*'),glob(FolderPath + rf'{fliers[1]}\\\\*leesa*')]}

    prgMsg('Collecting Trajectory data')

    # --- Create TRAJECTORY Data Dicts ---
    data_dicts_traj = []

    for i in range(len(dataPath_traj)):
        data_dict_traj = {}
        data_dict_traj = loadDictFromFile(dataPath_traj[i][0],data_dict_traj)
        data_dict_traj['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_traj['Epoch_esa'][0][i]) for i in range(len(data_dict_traj['Epoch_esa'][0]))])
        data_dicts_traj.append(data_dict_traj)

    Done(start_time)

    # --- Create ESA L1 Data Dicts ---
    prgMsg('Collecting ESA L1 data')
    data_dicts_template = {'eepaa':[],'iepaa':[],'leesa':[]}
    data_dicts_ESA_L1 = storeData(data_dicts_template,dataPath_L1ESA)
    Done(start_time)

    # --- Create ESA L2 Data Dicts ---
    prgMsg('Collecting ESA L2 data')
    data_dicts_template = {'eepaa': [], 'iepaa': [], 'leesa': []}
    data_dicts_ESA_L2 = storeData(data_dicts_template, dataPath_L2ESA)
    Done(start_time)

    # --- Create Distribution_Func Data Dicts ---
    prgMsg('Collecting Distribution Function data')
    data_dicts_template = {'eepaa': [], 'iepaa': [], 'leesa': []}
    data_dicts_dist = storeData(data_dicts_template,dataPath_Dist)
    Done(start_time)

    ########################
    # --- MAKE THE MOVIE ---
    ########################

    wInstr_movie = 'eepaa'

    prgMsg('Creating ESA movie (This may take awhile) \n')

    if wDataPlot[0] == 2: # Distribution Function
        data_dicts_movie = data_dicts_dist[wInstr_movie] # DONT CHANGE THE ORDER OF THIS LINE AND THE NEXT LINE
        wInstr_movie = 'Distribution_Function'
    elif wDataPlot[0] == 1:
        data_dicts_movie = data_dicts_ESA_L2[wInstr_movie]  # DONT CHANGE THE ORDER OF THIS LINE AND THE NEXT LINE
        wInstr_movie = 'Differential_Energy_Flux'
    else:
        wInstr_movie = wInstr_movie
        data_dicts_movie = data_dicts_ESA_L1[wInstr_movie]

    from plottingStyles import ESAmovie

    # Get the data
    Energy = data_dicts_movie[0]['Energy'][0]
    Pitch = data_dicts_movie[0]['Pitch_Angle'][0]
    geoAlt = [data_dicts_traj[0]['geoAlt'][0], data_dicts_traj[1]['geoAlt'][0]]
    geoLat = [data_dicts_traj[0]['geoLat'][0], data_dicts_traj[1]['geoLat'][0]]
    geoLong = [data_dicts_traj[0]['geoLong'][0], data_dicts_traj[1]['geoLong'][0]]
    EpochData = [np.array(data_dicts_movie[0]['Epoch_esa'][0]), np.array(data_dicts_movie[1]['Epoch_esa'][0])]
    locMin = np.abs(np.array(data_dicts_movie[0]['Epoch_esa'][0]) - targetTimes[0]).argmin()
    locMax = np.abs(np.array(data_dicts_movie[0]['Epoch_esa'][0]) - targetTimes[1]).argmin()
    locations = [i for i in range(locMin, locMax, frame_skips)]

    EpochData_dates = [np.array([pycdf.lib.tt2000_to_datetime(EpochData[0][i]).strftime("%H:%M:%S:%f") for i in range(len(EpochData[0]))]), np.array([pycdf.lib.tt2000_to_datetime(EpochData[1][i]).strftime("%H:%M:%S:%f") for i in range(len(EpochData[1]))]),]
    ESAData = [data_dicts_movie[0][wInstr_movie][0], data_dicts_movie[1][wInstr_movie][0]]

    ##################################################
    # --- RESIZE LOW FLYER DATA TO MATCH HIGHFLYER ---
    ##################################################
    highFlyerSize, lowFlyerSize = len(geoAlt[0]),len(geoAlt[1])

    # --- Append start Values ---
    no_of_points_start = int((EpochData[1][0] - EpochData[0][0]) / (rocketAttrs.MinorFrameTime))
    newAlt = [geoAlt[1][0] for i in range(no_of_points_start)]
    newLat = [geoLat[1][0] for i in range(no_of_points_start)]
    newLong = [geoLong[1][0] for i in range(no_of_points_start)]
    ESAfillVal = [[0 for engy in range(len(Energy))] for ptch in range(len(Pitch)) ]
    newESA = [ESAfillVal for i in range(no_of_points_start)]

    for i in range(len(geoAlt[1])):
        newAlt.append(geoAlt[1][i])
        newLat.append(geoLat[1][i])
        newESA.append(ESAData[1][i])
        newLong.append(geoLong[1][i])

    # --- Append the ending values ---
    remainingIndicies = highFlyerSize - (lowFlyerSize + no_of_points_start)

    for i in range(remainingIndicies):
        newAlt.append(geoAlt[1][-1])
        newLat.append(geoLat[1][-1])
        newLong.append(geoLong[1][-1])
        newESA.append(ESAfillVal)

    geoAlt[1], geoLat[1], geoLong[1] = np.array(newAlt), np.array(newLat), np.array(newLong)
    ESAData[1] = np.array(newESA)

    ##################################
    # --- FORMAT DATA FOR PLOTTING ---
    ##################################

    if wDataPlot[0] in [1, 2]: # Handles Distribution Function and DiffEFlux case

        # --- Velocity Space coordinate transformation ---
        modelData = ESAData[0]
        scienceVar = [[], []]

        # The vperp and vparallel coordinates never change, so pre-calculate them here for use later
        Vpars = np.array([np.cos(np.radians(Pitch[ptch])) * np.sqrt(2 * q0 * Energy[engy] / m_e) for engy in
                    range(ESAmovie.EnergyStart, ESAmovie.EnergyEnd) for ptch in range(len(modelData[0]))])
        Vperps = np.array([np.sin(np.radians(Pitch[ptch])) * np.sqrt(2 * q0 * Energy[engy] / m_e) for engy in
                     range(ESAmovie.EnergyStart, ESAmovie.EnergyEnd) for ptch in range(len(modelData[0]))])

        for i in range(2):

            if plotSpecificLocations:
                iterateThis = specific_locations
            else:
                iterateThis = range(len(modelData))

            for tme in iterateThis:
                scienceVar[i].append([ESAData[i][tme][ptch][engy] for engy in range(ESAmovie.EnergyStart, ESAmovie.EnergyEnd) for ptch in range(len(Pitch))])

        scienceVar = np.array(scienceVar)

        # --- interpolate and store the Data for each flyer ----

        # Note: I'm not interested in interpolating over the entire Vspace range. Much of the interesting stuff happens <2000eV
        # Reduce the X,Y linspace to only relvant data

        if normalizeToThermal:
            if wDataPlot[0] == 2:
                xlimitBot = ESAmovie.normalizedXLimits_dist[0] * ESAmovie.vth_e
                xlimitTop = ESAmovie.normalizedXLimits_dist[1] * ESAmovie.vth_e
                ylimitBot = ESAmovie.normalizedYLimits_dist[0] * ESAmovie.vth_e
                ylimitTop = ESAmovie.normalizedYLimits_dist[1] * ESAmovie.vth_e
            elif wDataPlot[0] == 1:
                xlimitBot = ESAmovie.normalizedXLimits_diff[0] * ESAmovie.vth_e
                xlimitTop = ESAmovie.normalizedXLimits_diff[1] * ESAmovie.vth_e
                ylimitBot = ESAmovie.normalizedYLimits_diff[0] * ESAmovie.vth_e
                ylimitTop = ESAmovie.normalizedYLimits_diff[1] * ESAmovie.vth_e
        else:
            if wDataPlot[0] == 2:
                ylimitBot = ESAmovie.YLimits_dist[0]
                ylimitTop = ESAmovie.YLimits_dist[1]
                xlimitBot = ESAmovie.XLimits_dist[0]
                xlimitTop = ESAmovie.XLimits_dist[1]
            elif wDataPlot[0] == 1:
                ylimitBot = ESAmovie.YLimits_diff[0]
                ylimitTop = ESAmovie.YLimits_diff[1]
                xlimitBot = ESAmovie.XLimits_diff[0]
                xlimitTop = ESAmovie.XLimits_diff[1]

        X = np.linspace(xlimitBot, xlimitTop, int(N/2))
        Y = np.linspace(ylimitBot, ylimitTop, int(N/2))

        X, Y = np.meshgrid(X, Y, indexing='xy')
        all_Z_data =[[],[]]

        for i in range(2):
            prgMsg(f'Interpolating for {rocketAttrs.rocketID[i]}\n')

            if plotSpecificLocations:
                iterateThis = range(len(scienceVar[i]))
            else:
                iterateThis = range(len(ESAData[i]))

            for tme in tqdm(iterateThis):
                interp = LinearNDInterpolator(list(zip(Vpars, Vperps)), scienceVar[i][tme], fill_value=fillvalue)
                all_Z_data[i].append(interp(X, Y))

            Done(start_time)

    elif wDataPlot[0] == 0:
        all_X_data = Pitch
        all_Y_data = Energy[ESAmovie.EnergyStart:ESAmovie.EnergyEnd]
        all_Z_data = [[], []]
        for i in range(2):
            for tme in range(1, len(ESAData[i])):
                all_Z_data[i].append(np.array([ESAData[i][tme][ptch][ESAmovie.EnergyStart:ESAmovie.EnergyEnd] for ptch in range(len(Pitch))]).transpose())


    ##########################
    # --- SETUP THE FIGURE ---
    ##########################
    prgMsg('Initializing Figure')
    fig = plt.figure()
    fig.set_figwidth(ESAmovie.figure_width)
    fig.set_figheight(ESAmovie.figure_height)

    # The overall plot shapeAlt/m_to_km
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1, 5])

    if wDataPlot[0] in [1,2]:
        # Shape of altitude and latvslong plots
        gsAlt = gs[0].subgridspec(nrows=1, ncols=3, width_ratios=[60, 60,1])
        ax1 = fig.add_subplot(gsAlt[0, 0])
        ax_lat = fig.add_subplot(gsAlt[0, 1])

        # Shape of ESA plots
        gsESA = gs[1].subgridspec(nrows=1, ncols=3, width_ratios=[60, 60, 1])
        ax2 = fig.add_subplot(gsESA[0])
        ax3 = fig.add_subplot(gsESA[1])

        # Give colorbar its own axis
        ax_colorbar = fig.add_subplot(gsESA[2])
    else:
        # Shape of altitude and latvslong plots
        gsAlt = gs[0].subgridspec(nrows=1, ncols=2, width_ratios=[1, 1])
        ax1 = fig.add_subplot(gsAlt[0, 0])
        ax_lat = fig.add_subplot(gsAlt[0, 1])

        # Shape of ESA plots
        gsESA = gs[1].subgridspec(nrows=2, ncols=2, hspace=0.3, height_ratios=[1, 1], width_ratios=[90, 1])

        # ESA plots
        ax2 = fig.add_subplot(gsESA[0, 0])
        ax3 = fig.add_subplot(gsESA[1, 0], sharex=ax2,sharey=ax2)

        # ColorBars
        ax_colorbar_high = fig.add_subplot(gsESA[0, 1])
        ax_colorbar_low = fig.add_subplot(gsESA[1, 1])

    ##############################
    # --- Initialize the plots ---
    ##############################

    title = fig.suptitle(f'{wInstr_movie.upper()} (T_th,e ={ESAmovie.T_e} eV)\n{EpochData_dates[0][0]}', **ESAmovie.title_style)

    # --- INITIALIZE TRAJECTORY DATA ---
    latAltMarker_high = ax1.scatter(geoLat[0][0], geoAlt[1][0], **ESAmovie.marker_style[0])
    latAltMarker_low = ax1.scatter(geoLat[0][0], geoAlt[1][0], **ESAmovie.marker_style[1])

    # Plot the altitude text
    latAltText_high = ax1.text(geoLat[0][0] - ESAmovie.textOffset_alt, geoAlt[0][0], f'{round(geoAlt[0][0], ESAmovie.altrounding)} km', ha=ESAmovie.text_alignment[1], **ESAmovie.textUTC_style_alt[0])
    latAltText_low = ax1.text(geoLat[1][0] - ESAmovie.textOffset_alt, geoAlt[1][0], f'{round(geoAlt[1][0], ESAmovie.altrounding)} km', ha=ESAmovie.text_alignment[1], **ESAmovie.textUTC_style_alt[1])

    latLongMarker_high = ax_lat.scatter(geoLat[0][0], geoLong[0][0], **ESAmovie.marker_style[0])
    latLongMarker_low = ax_lat.scatter(geoLat[1][0], geoLong[1][0], **ESAmovie.marker_style[1])

    latlongText_high = ax_lat.text(geoLat[0][0] - ESAmovie.textOffset_latlong[0][0], geoLong[0][0] - ESAmovie.textOffset_latlong[0][1], f'({round(geoLat[0][0], ESAmovie.altrounding)}' + '$^{\circ}$' + f', {round(geoLong[0][0], ESAmovie.altrounding)}' + '$^{\circ}$)', ha=ESAmovie.text_alignment[0], **ESAmovie.textUTC_style_lat[0])
    latlongText_low = ax_lat.text(geoLat[1][0] - ESAmovie.textOffset_latlong[1][0], geoLong[1][0] - ESAmovie.textOffset_latlong[1][1], f'({round(geoLat[1][0], ESAmovie.altrounding)}' + '$^{\circ}$' + f', {round(geoLong[1][0], ESAmovie.altrounding)}' + '$^{\circ}$)', ha=ESAmovie.text_alignment[1], **ESAmovie.textUTC_style_lat[1])

    for i in range(2):
        ax1.plot(geoLat[i], geoAlt[i], **ESAmovie.trajectory_style[i])  # Alttitude plot trajectory lines
        ax_lat.plot(geoLat[i], geoLong[i],**ESAmovie.trajectory_style[i]) # Lat vs Long plot trajectory lines

    # --- SET TRAJECTORY PLOTS LABEL PARAMETERS ---
    ax1.set_xlabel('Geographic Latitude', **ESAmovie.axes_label_traj_style)
    ax1.set_ylabel('Altitude [km]', **ESAmovie.axes_label_traj_style)
    ax_lat.set_xlabel('Geographic Latitude', **ESAmovie.axes_label_traj_style)
    ax_lat.set_ylabel('Geographic Longitude', **ESAmovie.axes_label_traj_style)
    trajAxes = [ax1,ax_lat]

    for i in range(2):
        trajAxes[i].tick_params(**ESAmovie.majorTicks_style)
        trajAxes[i].tick_params(**ESAmovie.minorTicks_style)

    # --- INITIALIZE ESA DATA ---
    subtitles = ['ACESII 36359', 'ACESII 36364']

    if wDataPlot[0] in [1, 2]:

        axes = [ax2, ax3]

        # Plot background color first before
        if normalizeToThermal:
            blackRadius = (2.5 * np.sin(np.radians(90)) * np.sqrt(2 * q0 * Energy[12] / m_e)) / ESAmovie.vth_e
        else:
            blackRadius = (2.5 * np.sin(np.radians(90)) * np.sqrt(2 * q0 * Energy[12] / m_e))

        center = (0, 0)

        for i in range(2):
            background_black_circle = Wedge(center, blackRadius, 340, 200, fc=(0.18995, 0.07176, 0.23217), edgecolor=None)
            axes[i].add_artist(background_black_circle)

        # --- INITIALIZE INTERPOLATED PLOTS ---
        if normalizeToThermal:
            if wDataPlot[0] == 2:
                xlimits = ESAmovie.normalizedXLimits_dist
                ylimits = ESAmovie.normalizedYLimits_dist

            elif wDataPlot[0] == 1:
                xlimits = ESAmovie.normalizedXLimits_diff
                ylimits = ESAmovie.normalizedYLimits_diff

            X = X / ESAmovie.vth_e
            Y = Y / ESAmovie.vth_e
        else:
            if wDataPlot[0] == 2:
                xlimits = ESAmovie.XLimits_dist
                ylimits = ESAmovie.YLimits_dist
            elif wDataPlot[0] == 1:
                xlimits = ESAmovie.XLimits_diff
                ylimits = ESAmovie.YLimits_diff

        xlimBot, xlimTop = xlimits[0], xlimits[1]  # values chosen based on Te = 1000eV
        ylimBot, ylimTop = ylimits[0], ylimits[1]

        if wDataPlot[0] == 2:
            # PLOT THE DISTRIBUTION DATA
            cmap1 = axes[0].pcolormesh(X, Y, all_Z_data[0][0], **ESAmovie.dist_cmap_style)
            cmap2 = axes[1].pcolormesh(X, Y, all_Z_data[1][0], **ESAmovie.dist_cmap_style)
            cmap_style = dict(cmap='turbo', vmin=ESAmovie.Vmin_dist, vmax=ESAmovie.Vmax_dist, norm='log')
        elif wDataPlot[0] == 1:
            # PLOT THE DIFFERENTIAL FLUX DATA
            cmap1 = axes[0].pcolormesh(X, Y, all_Z_data[0][0], **ESAmovie.diff_cmap_style)
            cmap2 = axes[1].pcolormesh(X, Y, all_Z_data[1][0], **ESAmovie.diff_cmap_style)
            cmap_style = dict(cmap='turbo', vmin=ESAmovie.Vmin_diff, vmax=ESAmovie.Vmax_diff, norm='log')

        # --- ADJUST PLOT LABELS ---
        for i in range(2):
            axes[i].set_xlim(xlimBot, xlimTop)
            axes[i].set_ylim(ylimBot, ylimTop)

            if normalizeToThermal:
                axes[i].set_xlabel('$v_{\parallel}/v_{th,e}$',**ESAmovie.axes_label_ESA_style)
                axes[i].set_ylabel('$v_{\perp}/v_{th,e}$',**ESAmovie.axes_label_ESA_style)
            else:
                axes[i].set_xlabel('$v_{\parallel}$', **ESAmovie.axes_label_ESA_style)
                axes[i].set_ylabel('$v_{\perp}$', **ESAmovie.axes_label_ESA_style)

            axes[i].set_title(subtitles[i], **ESAmovie.subtitle_style)
            axes[i].tick_params(**ESAmovie.majorTicks_style)
            axes[i].tick_params(**ESAmovie.minorTicks_style)

            # Get the "order of magnitude" text object
            if not normalizeToThermal:
                text = axes[i].yaxis.get_offset_text()
                text.set_size(ESAmovie.size_of_axis_orderOfMagnitude)
                text = axes[i].xaxis.get_offset_text()
                text.set_size(ESAmovie.size_of_axis_orderOfMagnitude)

        # --- INITIALIZE COLORBAR ---
        cbar = plt.colorbar(mappable=cmap1,cax=ax_colorbar, **ESAmovie.colorbar_style)

        if wDataPlot[0] == 2:
            cbar.set_label('Distribution Function [$m^{-6}s^{3}$]', **ESAmovie.cbar_label_style)
        elif wDataPlot[0] == 1:
            cbar.set_label('Differential Energy Flux [$cm^{-2}str^{-1}s^{-1} eV/eV$]', **ESAmovie.cbar_label_style)

        # --- COLORBAR LABEL PARAMS ---
        ax_colorbar.tick_params(**ESAmovie.majorTicks_colorbar_style)
        ax_colorbar.tick_params(**ESAmovie.minorTicks_colorbar_style)

        # --- INITIALIZED BLOCKING SHAPES ---
        center = (0, 0)
        theta1, theta2 = 0, 360
        radius = ESAmovie.radius_modifier * np.sin(np.radians(90)) * np.sqrt(2 * q0 * Energy[ESAmovie.EnergyEnd - 1] / m_e)

        if normalizeToThermal:
            radius = radius/ESAmovie.vth_e

        for i in range(2):

            # add sector lines
            for j in range(21):
                sector_circle = Wedge(center, blackRadius, 10*j + 345, 10*(j+1) + 345, fc=(1,1,1,0), edgecolor='black')
                axes[i].add_artist(sector_circle)

            semicircle = Wedge(center, radius, theta1, theta2, **ESAmovie.shapes_style)
            leftTriangle = plt.Polygon([[0.0, -0.0 * radius], [0.0, ylimBot], [-1 * ylimBot / np.tan(np.radians(15.5)), ylimBot]],**ESAmovie.shapes_style)
            rightTriangle = plt.Polygon([[0.0, -0.0 * radius], [0.0, ylimBot], [ylimBot / np.tan(np.radians(15.5)), ylimBot]],**ESAmovie.shapes_style)
            axes[i].add_artist(semicircle)
            axes[i].add_artist(leftTriangle)
            axes[i].add_artist(rightTriangle)

    else:
        vExtremes = [[1E-17, 1E-13], [1E-17, 1E-13]]
        cmap1 = ax2.pcolormesh(all_X_data, all_Y_data, all_Z_data[0][0],vmin=vExtremes[0][0], vmax=vExtremes[0][1], **ESAmovie.pcolormesh_style)
        cmap2 = ax3.pcolormesh(all_X_data, all_Y_data, all_Z_data[1][0],vmin=vExtremes[1][0], vmax=vExtremes[1][1], **ESAmovie.pcolormesh_style)

        axes = [ax2, ax3]

        # INITIALIZE COLORBARS
        axes_colorbars = [ax_colorbar_high, ax_colorbar_low]
        for i in range(2):
            cbar = plt.colorbar(**ESAmovie.colorbar_style, cax=axes_colorbars[i])
            cbar.minorticks_on()
            cbar.ax.tick_params(**ESAmovie.majorTicks_style)
            cbar.set_label('Counts', **ESAmovie.cbar_label_style)

        # --- ADJUST PLOT LABELS ---
        for i in range(2):
            axes[i].set_title(subtitles[i],**ESAmovie.subtitle_style)
            axes[i].tick_params(**ESAmovie.majorTicks_style)
            axes[i].tick_params(**ESAmovie.minorTicks_style)
            axes[i].set_ylabel('Energy [eV]', **ESAmovie.axes_label_ESA_style)
            axes[i].set_yscale('log')
            axes[i].xaxis.set_major_locator(MultipleLocator(ESAmovie.spacing_of_majorticks_ESA))
            axes[i].xaxis.set_major_formatter('{x:.0f}')

    plt.tight_layout()
    Done(start_time)
    prgMsg('Compiling the Movie')
    def animatePlot(i):

        # update Epoch title
        title.set_text(f'{wInstr_movie.upper()} (T_th,e ={ESAmovie.T_e} eV)\n{EpochData_dates[0][i]}')

        # update ESA DATA
        cmap1.set_array(all_Z_data[0][i])
        cmap2.set_array(all_Z_data[1][i])

        ##### UPDATE ATTITUDE DATA #####
        # --- high ---
        latAltMarker_high.set_offsets([geoLat[0][i], geoAlt[0][i]])
        latAltText_high.set_x(geoLat[0][i])
        latAltText_high.set_y(geoAlt[0][i])
        latAltText_high.set_text(f'{round(geoAlt[0][i], ESAmovie.altrounding)} km')

        latLongMarker_high.set_offsets([geoLat[0][i], geoLong[0][i]])
        latlongText_high.set_x(geoLat[0][i])
        latlongText_high.set_y(geoLong[0][i])
        latlongText_high.set_text(f'({round(geoLat[0][i], ESAmovie.altrounding)}' + '$^{\circ}$' + f', {round(geoLong[0][i], ESAmovie.altrounding)}' + '$^{\circ}$)')

        # --- low ---
        latAltMarker_low.set_offsets([geoLat[1][i],geoAlt[1][i]])
        latAltText_low.set_x(geoLat[1][i])
        latAltText_low.set_y(geoAlt[1][i])
        latAltText_low.set_text(f'{round(geoAlt[1][i], ESAmovie.altrounding)} km')

        latLongMarker_low.set_offsets([geoLat[1][i],geoLong[1][i]])
        latlongText_low.set_x(geoLat[1][i])
        latlongText_low.set_y(geoLong[1][i])
        latlongText_low.set_text(f'({round(geoLat[1][i], ESAmovie.altrounding)}' + '$^{\circ}$' + f', {round(geoLong[1][i], ESAmovie.altrounding)}' + '$^{\circ}$)')


    if plotSpecificLocations:
        locations = specific_locations
    elif plotSpecificTimeRange:
        # find the closest match in the high flyer data
        locMin = np.abs(np.array(data_dicts_movie[0]['Epoch_esa'][0]) - targetTimes[0]).argmin()
        locMax = np.abs(np.array(data_dicts_movie[0]['Epoch_esa'][0]) - targetTimes[1]).argmin()
        locations = [i for i in range(locMin, locMax, frame_skips)]
    else:
        locations = [i for i in range(0, int(len(EpochData_dates[0])), frame_skips)]  # NEEDS TO BE THE HIGH FLYER LENGTH

    anim = animation.FuncAnimation(fig=fig, func=animatePlot, interval= 1000/ESAmovie.fps, frames=locations)
    anim.save(f'C:\Data\ACESII\\trajectories\\trajectory_plots\\movies\\ACESII_{wInstr_movie}.mp4')
    Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
ACESIIplotting()
