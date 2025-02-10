# --- myImports.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: There are some common imports that every file uses. In order to de-clutter my code
# I can place these imports here. Only the imports which EVERY file uses will go here.

#########################
# --- IMPORTS IMPORTS ---
#########################
######################
# --- FROM IMPORTS ---
######################


# from spaceToolsLib.tools.epochTime import dateTimetoTT2000
from spaceToolsLib.setupFuncs.setupSpacepy import setupPYCDF

#####################
# --- SETUP PYCDF ---
#####################
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)