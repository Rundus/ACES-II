# --- WavePacketTimes.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Place to store the time windows to define (roughly) where
# the wave packets in the deltaB/deltaE data start and stop. The provide a function
# to provide that data conveniently. The time targets were chosen based on field-aligned data USING THE POYNTING FLUX DATA
# for the low flyer

import numpy as np
import datetime as dt

wavePacketTimes = np.array([
    # High Flyer
    [[],  # s1
     []], # s2

    # Low Flyer
    [[dt.datetime(2022, 11, 20, 17, 24, 45, 000),dt.datetime(2022, 11, 20, 17, 24, 46, 000)], # s1
     [dt.datetime(2022, 11, 20, 17, 24, 47, 000),dt.datetime(2022, 11, 20, 17, 24, 47, 500)], # s2
     [dt.datetime(2022, 11, 20, 17, 24, 47, 000),dt.datetime(2022, 11, 20, 17, 24, 47, 500)], # s2
     ] # s2
])