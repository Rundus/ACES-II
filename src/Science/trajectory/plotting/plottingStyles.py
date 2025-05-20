# --- plottingStyles.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Location to store the settings for the plots in ACESII_trajectory_plotting
import numpy as np


###############################################################
# --- PLOT 1: BOTH FLYERS, ALT VS LAT (w/ CONNECTING LINES) ---
###############################################################

class ALTvsLAT:

    figure_size = (13, 10)

    # high/low flyer plot
    trajectory_style = [dict(color='red', linewidth=1),dict(color='blue', linewidth=1)]
    plot_style_MAG = dict(color='blue', linewidth=1, alpha=0.0)
    axes_label_style = [dict(size=18, labelpad=10),dict(size=18, labelpad=10)]

    # ticks
    spacing_of_majorticks = 1
    majorTicks_style = [dict(labelsize=15, which='major', size=15, pad=5),dict(labelsize=15, which='major', size=15, pad=5)]
    minorTicks_style = [dict(labelsize=15, which='minor', size=10),dict(labelsize=15, which='minor', size=10)]

    # x-mark lines
    marker_style = [dict(marker='x', color='red', s=150), dict(marker='x', color='blue', s=150)]
    textUTC_style = [dict(size=14, color='red'),dict(size=14, color='blue')]

    # text alignement
    scatter_text_alignment = [['right', 'right', 'right', 'right', 'right', 'left', 'left', 'left', 'left','left', 'left'],['right', 'right', 'right', 'left', 'left', 'left', 'left']]
    scatter_text_alignment_offsets = [[-0.1, -0.1, -0.1, -0.1, -0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],[-0.1, -0.1, -0.1, 0.1, 0.1, 0.1, 0.1]]

    # magnetic lines
    mag_background_style = dict(color='black', linewidth=2, linestyle='-.', alpha=0.2)

    # Connecting lines
    connecting_lines_style = dict(color='green', linestyle="--", alpha=0.7)



################################################################
# --- PLOT 2: BOTH FLYERS, LAT VS LONG (w/ CONNECTING LINES) ---
################################################################
class LATvsLONG:

    figure_size = (13, 10)

    # Trajectory Plot
    trajectory_style = [dict(color='red', linewidth=1),dict(color='blue', linewidth=1)]
    axes_label_style = [dict(size=18, labelpad=10),dict(size=18, labelpad=10)]

    # invisible geomag plots
    plot_style_MAG = dict(color='blue', linewidth=1, alpha=0.0)

    # ticks
    spacing_of_majorticks = 1
    majorTicks_style = [dict(labelsize=15, which='major', size=15, pad=5),
                        dict(labelsize=15, which='major', size=15, pad=5)]
    minorTicks_style = [dict(labelsize=15, which='minor', size=10), dict(labelsize=15, which='minor', size=10)]

    # x-mark lines
    marker_style = [dict(marker='x', color='red', s=150), dict(marker='x', color='blue', s=150)]
    textUTC_style = [dict(size=14, color='red'), dict(size=14, color='blue')]

    # text alignement
    scatter_text_alignment = [['left', 'left', 'left', 'left', 'left', 'left', 'left', 'left', 'left', 'left', 'left'],
                              ['right', 'right', 'right', 'right', 'right', 'right', 'right']]
    scatter_text_alignment_offsets = [[0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15],
                                      [-0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, ]]


###########################################################################
# --- PLOT 3: BOTH FLYERS, ILAT,LAT VS ILONG,LONG (w/ CONNECTING LINES) ---
###########################################################################
class ILATvsILONG:

    targetProjectionAltitude = 100
    refLong = 16.020833
    refLat = 69.294167

    figureSize = (9, 9)

    trajectory_style = [dict(color='red', linewidth=1,label='High Flyer Trajectory'), dict(color='blue', linewidth=1,label='Low Flyer Trajectory')]
    intersections_style = [dict(color='salmon', linewidth=1,label='High Flyer - Geograhpic Projection'),
                           dict(color='deepskyblue', linewidth=1,label='Low Flyer - Geographic Projection')]
    axes_label_style = dict(size=18, labelpad=10)

    # invisible geomag plots
    plot_style_MAG = dict(color='blue', linewidth=1, alpha=0.0)

    # Grid
    gridMajor_style = dict(which='Major', alpha=0.5)
    gridMinor_style = dict(which='Minor', alpha=0.2)

    # ticks
    spacing_of_majorticks = 1
    spacing_of_majorticks_km = 50
    majorTicks_style = [dict(labelsize=15, which='major', size=15, pad=5), dict(labelsize=15, which='major', size=15, pad=5)]
    minorTicks_style = [dict(labelsize=15, which='minor', size=10), dict(labelsize=15, which='minor', size=10)]

    # x-mark lines
    marker_style = [dict(marker='x', color='red', s=150), dict(marker='x', color='blue', s=150)]
    marker_intersect_style = dict(marker='x', color='salmon', s=150)
    ha = ['left', 'right']
    textUTC_style = [dict(size=10, color='red',ha = ha[0]), dict(size=10, color='blue',ha=ha[1])]
    textUTC_intersect_style = dict(size=10, color='salmon', ha=ha[0])
    adjusts_km = [[1, 0], [-1 * 2, -1 * 8]]  # km adjustments to the timestamp labels
    adjusts = [[0.05, 0.05], [-0.05, -0.05]] # deg adjustments to the timestamp labels

    # Connecting lines
    connecting_lines_style = dict(color='green', linestyle="--", alpha=0.7)

############################################
# --- PLOT 4: ALT vs EPOCH (ESA OVERLAY) ---
############################################

class ESAoverlay:

    figure_height = 15
    figure_width = 23

    # "Zoom" into data
    lowFlyer_adjust = 1000

    ####################
    # --- FIRST AXIS ---
    ####################

    # pcolormesh
    pcolormesh_style = dict(vmin=3, vmax=80, cmap="turbo", alpha=0.8)

    # colorbar
    colorbar_style = dict(pad=0.15)
    cbar_label_style = dict(size=20)
    majorTicks_style_cbar = dict(labelsize=20)
    axes_label_style = dict(size=30, labelpad=5)
    tick_style = dict(labelsize=20)

    #####################
    # --- SECOND AXIS ---
    #####################

    # given in geodetic Latitude
    start = [71,70.8]
    end = [72.4,73.5]

    # high/low flyer plot
    trajectory_style = [dict(color='red', linewidth=2), dict(color='cyan', linewidth=2)]
    plot_style_MAG = dict(color='blue', linewidth=1, alpha=0.0)
    axes_label_style_twin = dict(size=20, labelpad=5)
    yaxis_lim = dict(bottom=8, top=1435)
    title_style = dict(size=40)
    trajectory_style_twin = [dict(color='white', linewidth=3), dict(color='cyan', linewidth=3)]


    # ticks
    spacing_of_majorticks_twin = 0.2

##########################
# --- PLOT 5: ESAMOVIE ---
##########################

class ESAmovie:

    # fps of the video
    fps = 5

    # figure size
    figure_height = 20
    figure_width = 36

    # "Zoom" into energy range
    EnergyStart = 12
    EnergyEnd = 41

    ########################
    # --- TRAJECT PLOT --- #
    ########################

    # high/low flyer plot
    trajectory_style = [dict(color='red', linewidth=2), dict(color='blue', linewidth=2)]
    axes_label_traj_style = dict(size=20, labelpad=5)
    title_style = dict(size=40)
    subtitle_style = dict(size=25)

    # ticks
    spacing_of_majorticks_Traj = 0.5
    majorTicks_style = dict(labelsize=20, which='major', size=15, pad=2)
    minorTicks_style = dict(labelsize=15, which='minor', size=10)
    polarTicks_degrees_style = dict(size=15 )
    polarTicks_radii_style = dict(size=50)

    # Altitude Text
    textOffset_alt = 0.15
    textOffset_latlong = [[-0.15,-0.15],[0.15, 0.25]]
    text_alignment = ['left','right']
    textUTC_style_lat = [dict(size=22, color='red'), dict(size=22, color='blue')]
    textUTC_style_alt = [dict(size=22, color='red'), dict(size=22, color='blue')]
    altrounding = 1 # how many decimals to round to

    # x-mark lines
    marker_style = [dict(marker='o', color='red', s=175), dict(marker='o', color='blue', s=175)]

    ##########################
    # --- ESA DATA PLOTS --- #
    ##########################

    spacing_of_majorticks_ESA = 10
    size_of_axis_orderOfMagnitude = 30

    #### if plotDistFunc ###
    T_e = 1000
    vth_e = (2*(T_e)*(1.602176565 * 10**(-19)) / (9.11 * 10**(-31)))**(1/2)
    Vmin_dist, Vmax_dist = 1E-18, 1E-14
    dist_cmap_style = dict(cmap='turbo', vmin=Vmin_dist, vmax=Vmax_dist, norm='log')
    normalizedXLimits_dist = [-1.25, 1.25] # values chosen based on Te = 1000eV for vthermal
    normalizedYLimits_dist = [-0.25, 1.25] # values chosen based on Te = 1000eV for vthermal
    XLimits_dist = [-3E7, 3E7]
    YLimits_dist = [-1E7, 3E7]

    #### if plotDiffEFlux ###
    Vmin_diff, Vmax_diff = 1E6, 1E9
    diff_cmap_style = dict(cmap='turbo', vmin=Vmin_diff, vmax=Vmax_diff, norm='log')
    normalizedXLimits_diff = [-1.5, 1.5]  # values chosen based on Te = 1000eV for vthermal
    normalizedYLimits_diff = [-0.25, 1.5]  # values chosen based on Te = 1000eV for vthermal
    XLimits_diff = [-3E7, 3E7]
    YLimits_diff = [-1E7, 3E7]

    # Plot elements
    shapes_style = dict(fc='white', edgecolor=None)
    radius_modifier = 0.95

    # --- if ESA counts data ---
    vmin = [1E-19, 1E-19]
    vmax = [1E-13, 1E-13]

    # --- pcolormesh ---
    pcolormesh_style = dict(cmap='turbo', alpha=1)

    # axes label spectrogram
    axes_label_ESA_style = dict(size=38, labelpad=5)


    ##################
    # --- COLORBAR ---
    ##################
    colorbar_style = dict(pad=-4,cmap='turbo')
    cbar_label_style = dict(size=35, labelpad=50,rotation=270)
    majorTicks_colorbar_style = dict(labelsize=35, which='major', size=15, pad=2)
    minorTicks_colorbar_style = dict(labelsize=25, which='minor', size=10)



class AttitudeMovie:

    fps = 30

    # Initial Condition data
    launcherSettings = np.array([[1, 352.4, 78.3, 136.3], [1, 353.2, 76.2, 138]])  # vector mag, azimuth, elevation, altitutde (in ft)





