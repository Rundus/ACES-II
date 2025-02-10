# --- myImports.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: There are some common imports that every file uses. In order to de-clutter my code
# I can place these imports here. Only the imports which EVERY file uses will go here.

#########################
# --- IMPORTS IMPORTS ---
#########################
import scipy
import itertools
import time
import matplotlib.pyplot as plt
import numpy as np
import os
import datetime as dt
######################
# --- FROM IMPORTS ---
######################


from spaceToolsLib.variables import u0,q0,m_e,cm_to_m,kB,ep0,Re, ion_dict
from spaceToolsLib.tools.coordinates import getCoordinateKeys
# from spaceToolsLib.tools.epochTime import dateTimetoTT2000
from spaceToolsLib.tools.interpolate import InterpolateDataDict
from spaceToolsLib.tools.diagnoistics import Done,prgMsg
from spaceToolsLib.setupFuncs.setupSpacepy import setupPYCDF
from spaceToolsLib.tools.colors import color
from spaceToolsLib.tools.CDF_load import loadDictFromFile,getInputFiles
from spaceToolsLib.tools.CDF_output import outputCDFdata

from tqdm import tqdm
from glob import glob
from os.path import getsize
from scipy.optimize import curve_fit
from copy import deepcopy

from data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from class_var_func import  L1_TRICE_Quick,L2_TRICE_Quick,L0_TRICE_Quick,L0_ACES_Quick,L1_ACES_Quick
from missionAttributes import ACES_mission_dicts, TRICE_mission_dicts

#####################
# --- SETUP PYCDF ---
#####################
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)