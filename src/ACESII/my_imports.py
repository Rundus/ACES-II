# --- my_imports.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: There are some common imports that every file uses. In order to de-clutter my code
# I can place these imports here. Only the imports which EVERY file uses will go here.

#################
# --- IMPORTS ---
#################
from spaceToolsLib.setupFuncs.setupSpacepy import setupPYCDF
from src.ACESII.data_paths import DataPaths
from src.ACESII.mission_attributes import ACESII
from copy import deepcopy
from glob import glob
from tqdm import tqdm
import datetime as dt
import os
import numpy as np
import spaceToolsLib as stl
import time


#####################
# --- SETUP PYCDF ---
#####################
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)