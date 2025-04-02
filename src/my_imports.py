# --- my_imports.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: There are some common imports that every file uses. In order to de-clutter my code
# I can place these imports here. Only the imports which EVERY file uses will go here.

#################
# --- IMPORTS ---
#################
from spaceToolsLib.setupFuncs.setupSpacepy import setupPYCDF
from copy import deepcopy
from glob import glob
from src.data_paths import DataPaths
from src.mission_attributes import ACESII
from tqdm import tqdm

import datetime as dt
import os
from os.path import getsize
import numpy as np
import spaceToolsLib as stl
import time


#####################
# --- SETUP PYCDF ---
#####################
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)