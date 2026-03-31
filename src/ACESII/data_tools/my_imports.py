# --- my_imports.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: There are some common imports that every file uses. In order to de-clutter my code
# I can place these imports here. Only the imports which EVERY file uses will go here.

#################
# --- IMPORTS ---
#################
from spaceToolsLib.setupFuncs.setupSpacepy import setupPYCDF

#####################
# --- SETUP PYCDF ---
#####################
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)


########################
# --- COMMON IMPORTS ---
########################
import time
import numpy as np
from src.ACESII.data_tools.data_classes import DataClasses
import spaceToolsLib as stl
from tqdm import tqdm
from glob import glob
from src.ACESII.data_tools.data_paths import DataPaths
from src.ACESII.mission_attributes import ACESII
import os
from copy import deepcopy
