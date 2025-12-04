

class L1toL2_FixedLPToggles:

    # --- FILE I/O ---
    modifier = ''
    inputPath_modifier = 'L1'  # e.g. 'L1' or 'L1'. It's the name of the broader input folder
    outputPath_modifier = 'L2'  # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder


class L1toL2_SweptLPToggles:

    # --- FILE I/O ---
    modifier = ''
    inputPath_modifier = 'L1'  # e.g. 'L1' or 'L1'. It's the name of the broader input folder
    outputPath_modifier = 'L2'  # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

    # --- RC effect downsample ---
    down_sample_RC_effect_bool = False
    keep_this_many_points = 1 # how many points to keep when you downsample RC effect

    # --- break into curves toggles ---
    breakIntoCurves = True
    targetVoltage_min = -1  # only care about voltage sweeps above this voltage value. Nominally -1
    digitalVariance = 5  # how much the digitized step point can vary when looking for the top and bottom of curves. nominally = 5
    indvEpochThresh = 15000000  # Value, in tt2000, that determines the time diff needed between epoch points to identify particular sweeps


class L2toL3_FixedLPToggles:

    # --- FILE I/O ---
    modifier = ''
    inputPath_modifier = 'L2'  # e.g. 'L1' or 'L1'. It's the name of the broader input folder
    outputPath_modifier = 'L3\Langmuir'  # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
    errorPath_modifier = 'calibration\LP_calibration'



class L2toL3_SweptLPToggles:

    # --- FILE I/O ---
    modifier = ''
    inputPath_modifier = 'L2'  # e.g. 'L1' or 'L1'. It's the name of the broader input folder
    outputPath_modifier = 'L3\Langmuir'  # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
    errorPath_modifier = 'calibration\LP_calibration'

    ###################
    ### SWEPT PROBE ###
    ###################
    # bad data: High Flyer 36
    # auroral case: 260
    # calm case:
    # Exponential case: 40, 100, 70,90, 150
    # Linear Cases: 130

    # --- TRADITIONAL FIT TOGGLES ---
    SECTION_SweptProbeniTe = False
    wSweeps = []  # [] --> all sweeps, [#1,#2,...] specific sweeps
    IgnoreTheseSweeps = [[22, 23], []]  # some particular sweeps in the data are bad, ignore them
    plotSpecificSweeps = False  # plot the sweeps in wSweeps
    traditionalFit = False
    showTeFittingMethod = False  # shows how Te is derived. Using the lowest Chisquare
    fitsStartHere = 0.7  # in volts. Start searching for the linear region here
    fitsEndHere = 1.8  # in volts. Finish searching for the linear region here
    showVspFittingMethod = False
    alternateData = False
    wSet = 1  # nominally == 0 or 1, determines keeping every odd sweep or even sweep

    # ^^^ I noticed a pattern in the data: every other sweep follows its own pattern,
    # almost like two datasets are being sampled. Splitting the data in 2 may reveal something

    # --- GRID SEARCH Vsp FIT TOGGLES ---
    gridSearchFit = False  # the toggle parameters are in LP_gridSearch_toggles
    showGridSearchFitting = False  # shows the best gridsearch fit