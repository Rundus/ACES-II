# --- LP_gridSearch_toggles.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Place to put the toggles for the LP GridSearch technique

import numpy as np

#############################
# --- GRID SEARCH TOGGLES ---
#############################

# --- Vsp ---
searchStartsHere = 0.4 # use this as the minimum transition xData value for the fits. a1_range minimum must be larger than this

# number of points used to perform grid search
N_nTe = 10
N_Vsp = 10
N_Te = 100


a0_range = np.linspace(10E6, 10E12, N_nTe)
a1_range = np.linspace(0.7, 2, N_Vsp)
a2_range = np.linspace(0.001, 5, N_Te)

transGuess = [1E15, 1, 0.2]
transBounds = ([1E10, 0, 0], [np.Inf, 2.5, 100]) # Bounds for: a0 = nT_e, a1 = Vsp, a2 = Te
satGuess = [1E15, 1, 0.2]
satBounds = ([1E10, 0, 0], [np.Inf, 2.5, 100]) # Bounds for: a0 = nT_e, a1 = Vsp, a2 = Te

# number of points the electron saturation region will loop through to check for best ChiSquare
satNumPoints = [3]

# Transition region fitting params
transFitParameters = dict(maxfev=1000000, bounds=transBounds, p0=transGuess)

# saturation region Fitting params
satFitParameters = dict(maxfev=1000000, bounds=satBounds, p0=satGuess)
