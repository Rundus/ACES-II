# --- csv_to_cdf_attitude.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Turn the .cdf files of the TRICE attitude data into cdf files


# --- --- --- --- ---
from src.my_imports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# --- Alt vs Lat/geoLat ---
plotALTvsLAT = False

# Use connecting lines in all the relevant plots?
connectingLines = True

useKM = False


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

    from plottingStyles import ILATvsILONG

    fig = plt.figure(figsize=ILATvsILONG.figureSize)
    ax1 = fig.add_subplot(111)
    ax1.minorticks_on()
    ax1.xaxis.set_major_formatter('{x:.0f}')

    ax1.tick_params(**ILATvsILONG.majorTicks_style[0])
    ax1.tick_params(**ILATvsILONG.minorTicks_style[0])

    plt.grid(True, **ILATvsILONG.gridMajor_style)
    plt.grid(True, **ILATvsILONG.gridMinor_style )

    # Plot the Trajectory Data
    if useKM:
        mod_file_ID = 'km'

        # AXIS LABELS
        ax1.set_xlabel(f'(+E|-W) [km]', **ILATvsILONG.axes_label_style)
        ax1.set_ylabel(f'(+N|-S) [km]', **ILATvsILONG.axes_label_style)

        # TICKS
        major_xticks = np.arange(-50, 50, ILATvsILONG.spacing_of_majorticks_km)
        major_yticks = np.arange(0, 550, ILATvsILONG.spacing_of_majorticks_km)
        ax1.set_xticks(major_xticks)
        ax1.set_yticks(major_yticks)
        ax1.xaxis.set_minor_locator(AutoMinorLocator())

        # GROUND TRACK TRAJECTS
        for i in range(len(dataPath_traj)):
            ax1.plot(data_dicts_traj[i]['geoLong_km'][0], data_dicts_traj[i]['geoLat_km'][0], **ILATvsILONG.trajectory_style[i])
            ax1.plot(data_dicts_traj[i]['geoILong_km'][0], data_dicts_traj[i]['geoILat_km'][0],**ILATvsILONG.intersections_style[i])

            # Plot X marks
            ax1.scatter(timeTargetsData['geoLong_km'][i], timeTargetsData['geoLat_km'][i], **ILATvsILONG.marker_style[i])

            # UTC Text
            for k in range(len(timeTargetsEpoch[i])):
                ax1.text(timeTargetsData['geoLong_km'][i][k] + ILATvsILONG.adjusts_km[i][0],
                         timeTargetsData['geoLat_km'][i][k] + ILATvsILONG.adjusts_km[i][1],
                         timeTargetsEpoch[i][k], **ILATvsILONG.textUTC_style[i])

        if connectingLines:

            # plot the time tags for the High Flyer Only but plot the lines on the low flyer only
            for k in range(len(timeTargetsEpoch[1])):
                ax1.text(timeTargetsData['geoILong_km'][0][k] + ILATvsILONG.adjusts_km[0][0],timeTargetsData['geoILat_km'][0][k] + ILATvsILONG.adjusts_km[0][1],timeTargetsEpoch[1][k], **ILATvsILONG.textUTC_intersect_style)

                # Plot X marks
                ax1.scatter(timeTargetsData['geoILong_km'][0][k], timeTargetsData['geoILat_km'][0][k], **ILATvsILONG.marker_intersect_style)

                # Plot the connecting lines
                x_values = [timeTargetsData['geoILong_km'][0][k],timeTargetsData['geoILong_km'][1][k]]
                y_values = [timeTargetsData['geoILat_km'][0][k], timeTargetsData['geoILat_km'][1][k]]

                ax1.plot(x_values, y_values, **ILATvsILONG.connecting_lines_style)

    # TYPICALLY LAT/LONG
    else:
        mod_file_ID = 'latlong'
        ax1.set_xlabel('Geographic Long', **ILATvsILONG.axes_label_style)
        ax1.set_ylabel('Geographic Lat', **ILATvsILONG.axes_label_style)

        # TICKS
        major_xticks = np.arange(-50, 50, ILATvsILONG.spacing_of_majorticks)
        major_yticks = np.arange(0, 550, ILATvsILONG.spacing_of_majorticks)
        ax1.set_xticks(major_xticks)
        ax1.set_yticks(major_yticks)
        ax1.xaxis.set_minor_locator(AutoMinorLocator())

        # GROUND TRACK TRAJECTS
        for i in range(len(dataPath_traj)):
            ax1.plot(data_dicts_traj[i]['geoLong'][0], data_dicts_traj[i]['geoLat'][0],
                     **ILATvsILONG.trajectory_style[i])
            ax1.plot(data_dicts_traj[i]['geoILong'][0], data_dicts_traj[i]['geoILat'][0],
                     **ILATvsILONG.intersections_style[i])

            # Plot X marks
            ax1.scatter(timeTargetsData['geoLong'][i], timeTargetsData['geoLat'][i],
                        **ILATvsILONG.marker_style[i])

            # UTC Text
            for k in range(len(timeTargetsEpoch[i])):
                ax1.text(timeTargetsData['geoLong'][i][k] + ILATvsILONG.adjusts[i][0],
                         timeTargetsData['geoLat'][i][k] + ILATvsILONG.adjusts[i][1],
                         timeTargetsEpoch[i][k], **ILATvsILONG.textUTC_style[i])

            # GEOMAGNETIC AXIS
            ax2 = ax1.twinx()
            ax2.plot(data_dicts_traj[0]['geoLong'][0], data_dicts_traj[0]['geomagLat'][0],**ILATvsILONG.plot_style_MAG)
            # ticks major/minor
            ax2.minorticks_on()
            ax2.yaxis.set_major_locator(MultipleLocator(ILATvsILONG.spacing_of_majorticks))
            ax2.yaxis.set_major_formatter('{x:.0f}')
            ax2.yaxis.set_minor_locator(AutoMinorLocator())
            ax2.tick_params(**ILATvsILONG.majorTicks_style[1])
            ax2.tick_params(**ILATvsILONG.minorTicks_style[1])
            ax2.set_ylabel('Geomagnetic Lat', **ILATvsILONG.axes_label_style)

        if connectingLines:

            # plot the time tags for the High Flyer Only but plot the lines on the low flyer only
            for k in range(len(timeTargetsEpoch[1])):
                ax1.text(timeTargetsData['geoILong'][0][k] + ILATvsILONG.adjusts[0][0],
                         timeTargetsData['geoILat'][0][k] + ILATvsILONG.adjusts[0][1], timeTargetsEpoch[1][k],
                         **ILATvsILONG.textUTC_intersect_style)

                # Plot X marks
                ax1.scatter(timeTargetsData['geoILong'][0][k], timeTargetsData['geoILat'][0][k],
                            **ILATvsILONG.marker_intersect_style)

                # Plot the connecting lines
                x_values = [timeTargetsData['geoILong'][0][k], timeTargetsData['geoILong'][1][k]]
                y_values = [timeTargetsData['geoILat'][0][k], timeTargetsData['geoILat'][1][k]]

                ax1.plot(x_values, y_values, **ILATvsILONG.connecting_lines_style)

    ax1.legend()
    plt.title(f'EEPAA \n B-Field Projection Altitude: {ILATvsILONG.targetProjectionAltitude} km \n Origin: {ILATvsILONG.refLat}$^\circ$N {ILATvsILONG.refLong}$^\circ$E')
    fig.savefig(rf'D:\Data\ACESII\trajectories\trajectory_plots\BProjection_{mod_file_ID}.png')

