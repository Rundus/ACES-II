class FixedLPToggles:
    ###################
    ### FIXED PROBE ###
    ###################
    SECTION_calculateFixedni = True
    fixed_Ti_assumed = True  # IF FALSE use IRI model
    tromsoScales = [1 / 50, 1 / 50]  # values used to make LP density match ~ 5.7E4 cm^-4 at the E-Region
    Ti_assumed = 0.1  # assuming an ion temperature (in eV)
    unit_conversion = 1E9  # 1 for Amps 10**9 for nano amps, etc


class SweptLPToggles:
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