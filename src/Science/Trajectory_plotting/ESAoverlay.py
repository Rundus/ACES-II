# --- csv_to_cdf_attitude.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Turn the .cdf files of the TRICE attitude data into cdf files


# --- --- --- --- ---
from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

plotESAoverlay = False
wRocket_overlay = 4
wPitches_overlay = [1, 10, 19] # 0deg, 90deg, 180deg
wInstr_overlay = 'eepaa'


# --- --- --- ---
# --- import ---
# --- --- --- ---
from matplotlib import pyplot as plt
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

    rocketAttrs,b,c = ACES_mission_dicts()

    # Trajectory Data
    trajFolderPath = f'{ACES_data_folder}trajectories\\'
    dataPath_traj = [glob(trajFolderPath + rf'{fliers[0]}\\*_ILat_ILong*'),
                     glob(trajFolderPath + rf'{fliers[1]}\\\\*_ILat_ILong*')]

    # L1 ESA data
    FolderPath = f'{ACES_data_folder}L1\\'
    dataPath_L1ESA = {'eepaa':[glob(FolderPath + rf'{fliers[0]}\\*eepaa_*'), glob(FolderPath + rf'{fliers[1]}\\\\*eepaa_*')],
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

    ####################################################################
    # --- find the points 100s, 200s ... 600s in the Trajectory data ---
    ####################################################################
    prgMsg('Finding Target Times')
    timeTargets = [[100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600], [100, 150, 200, 250, 300, 350, 400]]

    timeTargetsIndices = [[], []]
    timeTargetsEpoch = []
    for i in range(len(timeTargets)):
        for timeTarg in timeTargets[i]:
            timeTargetsIndices[i].append(np.abs(data_dicts_traj[i]['Epoch_esa'][0] - (rocketAttrs.Launch_Times[i] + timeTarg * (10 ** (9)))).argmin())

        target_in_Epoch = [pycdf.lib.tt2000_to_datetime(data_dicts_traj[i]['Epoch_esa'][0][index]) for index in timeTargetsIndices[i]]
        timeTargetsEpoch.append([targetTime.strftime("%H:%M:%S") for targetTime in target_in_Epoch])

    timeTargetsData = {
        'Epoch': [[pycdf.lib.tt2000_to_datetime(data_dicts_traj[0]['Epoch_esa'][0][index]).time().strftime("%H:%M:%S") for index in timeTargetsIndices[0]],
                  [pycdf.lib.tt2000_to_datetime(data_dicts_traj[1]['Epoch_esa'][0][index]).time().strftime("%H:%M:%S") for index in timeTargetsIndices[1]]],
        'geoAlt': [[data_dicts_traj[0]['geoAlt'][0][index] for index in timeTargetsIndices[0]],
                   [data_dicts_traj[1]['geoAlt'][0][index] for index in timeTargetsIndices[1]]],
        'geoLat': [[data_dicts_traj[0]['geoLat'][0][index] for index in timeTargetsIndices[0]],
                   [data_dicts_traj[1]['geoLat'][0][index] for index in timeTargetsIndices[1]]],
        'geoLong': [[data_dicts_traj[0]['geoLong'][0][index] for index in timeTargetsIndices[0]],
                    [data_dicts_traj[1]['geoLong'][0][index] for index in timeTargetsIndices[1]]],
        'geomagLat': [[data_dicts_traj[0]['geomagLat'][0][index] for index in timeTargetsIndices[0]],
                      [data_dicts_traj[1]['geomagLat'][0][index] for index in timeTargetsIndices[1]]],
        'geomagLong': [[data_dicts_traj[0]['geomagLong'][0][index] for index in timeTargetsIndices[0]],
                       [data_dicts_traj[1]['geomagLong'][0][index] for index in timeTargetsIndices[1]]],
        'geoLat_km': [[data_dicts_traj[0]['geoLat_km'][0][index] for index in timeTargetsIndices[0]],
                      [data_dicts_traj[1]['geoLat_km'][0][index] for index in timeTargetsIndices[1]]],
        'geoLong_km': [[data_dicts_traj[0]['geoLong_km'][0][index] for index in timeTargetsIndices[0]],
                       [data_dicts_traj[1]['geoLong_km'][0][index] for index in timeTargetsIndices[1]]],
        'geoILat': [[data_dicts_traj[0]['geoILat'][0][index] for index in timeTargetsIndices[0]],
                    [data_dicts_traj[1]['geoILat'][0][index] for index in timeTargetsIndices[1]]],
        'geoILong': [[data_dicts_traj[0]['geoILong'][0][index] for index in timeTargetsIndices[0]],
                     [data_dicts_traj[1]['geoILong'][0][index] for index in timeTargetsIndices[1]]],
        'geoILat_km': [[data_dicts_traj[0]['geoILat_km'][0][index] for index in timeTargetsIndices[0]],
                       [data_dicts_traj[1]['geoILat_km'][0][index] for index in timeTargetsIndices[1]]],
        'geoILong_km': [[data_dicts_traj[0]['geoILong_km'][0][index] for index in timeTargetsIndices[0]],
                        [data_dicts_traj[1]['geoILong_km'][0][index] for index in timeTargetsIndices[1]]]}
    Done(start_time)


    #######################
    # --- MAKE THE PLOT ---
    #######################


    from plottingStyles import ESAoverlay

    data_dicts_ESA = data_dicts_ESA_L1[wInstr_overlay]

    prgMsg('Plotting ESA overlay')

    #############################
    # --- "Zoom" into Dataset ---
    #############################

    # index locations of specific latitude
    start = np.abs(data_dicts_traj[wRocket_overlay-4]['geoLat'][0] -ESAoverlay.start[wRocket_overlay-4]).argmin()
    end = np.abs(data_dicts_traj[wRocket_overlay-4]['geoLat'][0] -ESAoverlay.end[wRocket_overlay-4]).argmin()

    # Get a subset of the data
    Energy = data_dicts_ESA[0]['Energy'][0]
    Pitch = data_dicts_ESA[0]['Pitch_Angle'][0]
    EpochData = data_dicts_ESA[wRocket_overlay-4]['Epoch_esa'][0][start:end]
    ESAData = data_dicts_ESA[wRocket_overlay-4][wInstr_overlay][0][start:end]

    ##########################
    # --- PLOT Spectrogram ---
    ##########################

    fig, ax = plt.subplots(3)

    fig.set_figheight(ESAoverlay.figure_height)
    fig.set_figwidth(ESAoverlay.figure_width)

    for i in range(len(wPitches_overlay)):

        # Trajectory and ESA data
        X = data_dicts_traj[wRocket_overlay-4]['geoLat'][0][start:end]
        Y = Energy
        Z = np.array([ESAData[j][wPitches_overlay[i]] for j in range(len(ESAData))]).transpose() # Get a slice in pitch data

        # PCOLORMESH Plot
        cmap = ax[i].pcolormesh(X, Y, Z, **ESAoverlay.pcolormesh_style)
        cbar = plt.colorbar(cmap, ax = ax[i])
        cbar.ax.tick_params(**ESAoverlay.majorTicks_style_cbar)
        cbar.set_label('Counts', **ESAoverlay.cbar_label_style)
        ax[i].set_ylabel('Energy [eV]', **ESAoverlay.axes_label_style)

        ax[i].minorticks_on()
        ax[i].yaxis.set_tick_params(**ESAoverlay.tick_style)
        ax[i].xaxis.set_tick_params(**ESAoverlay.tick_style)

        if i == 2:
            ax[i].set_xlabel('Geodetic Lat [deg]', **ESAoverlay.axes_label_style)

        ax[i].set_ylim(**ESAoverlay.yaxis_lim)

        # title
        ax[i].set_title(f'{wInstr_overlay}, Pitch = {Pitch[wPitches_overlay[i]]}$^\circ$', **ESAoverlay.title_style)

        ##################
        # --- Alt AXIS ---
        ##################
        ax1 = ax[i].twinx()
        ax1.plot(data_dicts_traj[wRocket_overlay - 4]['geoLat'][0][start:end], data_dicts_traj[wRocket_overlay - 4]['geoAlt'][0][start:end], **ESAoverlay.trajectory_style_twin[wRocket_overlay - 4])
        ax1.set_ylabel('Alt [km]', **ESAoverlay.axes_label_style_twin)
        ax1.minorticks_on()
        ax1.yaxis.set_tick_params(**ESAoverlay.tick_style)
        ax1.xaxis.set_tick_params(**ESAoverlay.tick_style)
        plt.tight_layout()

    fig.savefig(rf'D:\Data\ACESII\trajectories\trajectory_plots\{fliers[wRocket_overlay - 4]}\{wInstr_overlay}_overlay_ALTvsLAT.png')
    Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
ACESIIplotting()
