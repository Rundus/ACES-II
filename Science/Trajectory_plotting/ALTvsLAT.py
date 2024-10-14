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

# Use connecting lines in all the relevant plots?
connectingLines = True



# --- --- --- ---
# --- import ---
# --- --- --- ---
from matplotlib import pyplot as plt
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

    rocketAttrs,b,c = ACES_mission_dicts()

    # Trajectory Data
    trajFolderPath = f'{ACES_data_folder}trajectories\\'
    dataPath_traj = [glob(trajFolderPath + rf'{fliers[0]}\\*_ILat_ILong*'),
                     glob(trajFolderPath + rf'{fliers[1]}\\\\*_ILat_ILong*')]

    prgMsg('Collecting Trajectory data')

    # --- Create TRAJECTORY Data Dicts ---
    data_dicts_traj = []

    for i in range(len(dataPath_traj)):
        data_dict_traj = {}
        data_dict_traj = loadDictFromFile(dataPath_traj[i][0],data_dict_traj)
        data_dict_traj['Epoch_esa'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_traj['Epoch_esa'][0][i]) for i in range(len(data_dict_traj['Epoch_esa'][0]))])
        data_dicts_traj.append(data_dict_traj)

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

    from plottingStyles import ALTvsLAT

    # --- Plot the Data ---

    prgMsg('Plotting ALTvsLAT Data')
    fig = plt.figure(figsize=ALTvsLAT.figure_size)

    #########################
    # --- geographic data ---
    #########################

    ax1 = fig.add_subplot(111)

    # Plot trajectories and UTC labels
    for i in range(len(dataPath_traj)):
        ax1.plot(data_dicts_traj[i]['geoLat'][0], data_dicts_traj[i]['geoAlt'][0], **ALTvsLAT.trajectory_style[i]) # Trajectories
        ax1.scatter(timeTargetsData['geoLat'][i], timeTargetsData['geoAlt'][i], **ALTvsLAT.marker_style[i]) # UTC Text

    # ticks major/minor
    ax1.minorticks_on()
    ax1.set_ylabel('Altitude (km)', **ALTvsLAT.axes_label_style[0])
    ax1.set_xlabel('Geographic Latitude', **ALTvsLAT.axes_label_style[0])
    ax1.xaxis.set_major_locator(MultipleLocator(ALTvsLAT.spacing_of_majorticks))
    ax1.xaxis.set_major_formatter('{x:.0f}')
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.tick_params(**ALTvsLAT.majorTicks_style[0])
    ax1.tick_params(**ALTvsLAT.minorTicks_style[0])

    # Plot the text
    for j in range(2):
        for i, txt in enumerate(timeTargetsEpoch[j]):
            ax1.text(timeTargetsData['geoLat'][j][i] + ALTvsLAT.scatter_text_alignment_offsets[j][i],timeTargetsData['geoAlt'][j][i], txt, ha=ALTvsLAT.scatter_text_alignment[j][i], **ALTvsLAT.textUTC_style[j])

    # plot connecting lines
    if connectingLines:
        highFlyerOffset = 2
        for i in range(len(timeTargetsData['geoLat'][1])):
            x_values = [timeTargetsData['geoLat'][1][i], timeTargetsData['geoLat'][0][i+highFlyerOffset]]
            y_values = [timeTargetsData['geoAlt'][1][i], timeTargetsData['geoAlt'][0][i+highFlyerOffset]]
            ax1.plot(x_values, y_values, **ALTvsLAT.connecting_lines_style)

    ##########################
    # --- geomagnetic data ---
    ##########################
    ax2 = ax1.twiny()
    ax2.plot(data_dicts_traj[0]['geomagLat'][0], data_dicts_traj[0]['geoAlt'][0], **ALTvsLAT.plot_style_MAG) # High Flyer
    ax2.plot(data_dicts_traj[1]['geomagLat'][0], data_dicts_traj[1]['geoAlt'][0], **ALTvsLAT.plot_style_MAG)  # Low Flyer

    # ticks major/minor
    ax2.minorticks_on()
    ax2.xaxis.set_major_locator(MultipleLocator(ALTvsLAT.spacing_of_majorticks))
    ax2.xaxis.set_major_formatter('{x:.0f}')
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.tick_params(**ALTvsLAT.majorTicks_style[1])
    ax2.tick_params(**ALTvsLAT.minorTicks_style[1])
    ax2.set_xlabel('Geomagnetic Latitude',**ALTvsLAT.axes_label_style[1])

    ax3 = ax1.twinx()
    ax3.get_yaxis().set_visible(False)
    slope = -2.378842439667985558915 # corresponds to line with -78.13deg inclination
    for i in range(11):
        ax3.axline(xy1=(69+i*0.5,-1*slope), xy2=(69 + (0.5*(i+1)), 0),**ALTvsLAT.mag_background_style)

    fig.savefig(r'D:\Data\ACESII\trajectories\trajectory_plots\Altitude_vs_lat_vs_maglat.png')
    Done(start_time)



