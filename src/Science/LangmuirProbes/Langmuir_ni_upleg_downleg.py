# --- Langmuir_ni_upleg_downleg---
# Use the ACES-II LP data to analyze upleg/downleg n_i data from the
# payloads and try to estimate a nominal n_e vs altitude for the whole flight

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal

from src.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- OutputData ---
outputData = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.Science.LangmuirProbes.toggles import FloatingPotentialToggles as fToggles


def langmuir_upleg_downleg():