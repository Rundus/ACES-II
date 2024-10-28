# --- dispersionAttributes.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Location to store the information for the dispersion times
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt


###############################################################
# --- PLOT 1: BOTH FLYERS, ALT VS LAT (w/ CONNECTING LINES) ---
###############################################################
def boxRemove(Data, EnergyMin, EnergyMax, TimeMin, TimeMax, Energy, Time):

    # perform the removal
    newData = np.array(Data)
    for tme in range(len(Time)):
        for engy in range(len(Energy)):
            if Energy[engy] >= EnergyMin and Energy[engy] <= EnergyMax and Time[tme] >= TimeMin and Time[tme] <= TimeMax:
                for ptch in range(len(newData[0])):
                    newData[tme][ptch][engy] = 0

    return newData

def diagonalRemove(Data, EnergyStartPoint ,TimeStart, TimeEnd, Energy, Time, upper):

    # find the indices of energy cor
    # responding to EnergyMin,EnergyMax
    Esp = np.abs(Energy - EnergyStartPoint).argmin()

    # find time indicies of Tmin, Tmax
    Tmin = np.abs(Time - TimeStart).argmin()
    Tmax = np.abs(Time - TimeEnd).argmin()
    k = (Tmin - Esp)

    # remove any super large negative fillvals
    newData = np.array(Data).clip(min=0)

    # Construct the Diagonal Mask
    from copy import deepcopy
    zeroDiagonal_mask = np.array(deepcopy(newData)).T

    for ptch in range(21): # loop through all pitch angles
        pitchSlice = np.array(Data[:, ptch, :]).T

        # get the upper diagonal elements that will mask out the data later
        if upper:  # zero the upper section of the array
            zeroDiagonal_mask[:, ptch, :] = np.triu(pitchSlice, k)
        else: # zero the lower section of the array
            zeroDiagonal_mask[:, ptch, :] = pitchSlice - np.triu(pitchSlice, k+1)

    # Construct the adjusted mask - LEFT
    # description: remove vertical chunks from left or right of the mask based on Tmin and Tmax
    zeroDiagonal_mask_left = boxRemove(Data=zeroDiagonal_mask.T,
                                            EnergyMin=0,
                                            EnergyMax=20000,
                                            TimeMin=0,
                                            TimeMax=TimeStart,
                                            Energy=Energy,
                                            Time=Time)

    # right
    zeroDiagonal_mask_right = boxRemove(Data=zeroDiagonal_mask_left,
                                            EnergyMin=0,
                                            EnergyMax=20000,
                                            TimeMin=TimeEnd,
                                            TimeMax=Time[-1],
                                            Energy=Energy,
                                            Time=Time)
    # subtract mask from original data
    newData = np.array(newData - zeroDiagonal_mask_right)

    return newData

class dispersionAttributes(object):

    # taken from the 10deg bin on eepaa_fullCal. These are the most noticable
    keyDispersionTimes = [
        [dt.datetime(2022, 11, 20, 17, 24, 56, 000000), dt.datetime(2022, 11, 20, 17, 24, 57, 600000)],# s1
        [dt.datetime(2022, 11, 20, 17, 24, 57, 630000), dt.datetime(2022, 11, 20, 17, 24, 59, 186000)],# s2
        [dt.datetime(2022, 11, 20, 17, 24, 58, 800000), dt.datetime(2022, 11, 20, 17, 25, 00, 50000)],# s3
        [dt.datetime(2022, 11, 20, 17, 24, 59, 965000), dt.datetime(2022, 11, 20, 17, 25, 00, 758000)],# s4  right before Inverted V
        [dt.datetime(2022, 11, 20, 17, 25, 00, 501000), dt.datetime(2022, 11, 20, 17, 25, 1, 350000)], # s5 inside inverted V's left edge
        [dt.datetime(2022, 11, 20, 17, 25, 0, 850000),  dt.datetime(2022, 11, 20, 17, 25, 2, 0000)], # s6 under V, it's faint
        [dt.datetime(2022, 11, 20, 17, 25, 3, 800000), dt.datetime(2022, 11, 20, 17, 25, 5, 200000)], # s7 in choatic region
        [dt.datetime(2022, 11, 20, 17, 25, 5, 450000), dt.datetime(2022, 11, 20, 17, 25, 5, 920000)], # s8 VERTICAL BAR that doesn't look like it's in an inverted V
        [dt.datetime(2022, 11, 20, 17, 25, 5, 806000), dt.datetime(2022, 11, 20, 17, 25, 6, 506000)], # s9 Last dispersive feature in the chaotic region
        [dt.datetime(2022, 11, 20, 17, 25, 6, 760000), dt.datetime(2022, 11, 20, 17, 25, 7, 606000)], # s10 Very Strong Wave-like feature on north end of Inverted V
        [dt.datetime(2022, 11, 20, 17, 25, 7, 641000), dt.datetime(2022, 11, 20, 17, 25, 8, 406000)], # s11 Smaller feature, similar to s1
        [dt.datetime(2022, 11, 20, 17, 24, 18, 000), dt.datetime(2022, 11, 20, 17, 24, 19, 000)], # s12 first STEB observed. Dispersion underneath inverted V
        [dt.datetime(2022, 11, 20, 17, 24, 21, 000), dt.datetime(2022, 11, 20, 17, 24, 22, 000)], # s13 STEB northward of first aurora observed. Faint
        [dt.datetime(2022, 11, 20, 17, 24, 31, 00000), dt.datetime(2022, 11, 20, 17, 24, 32, 271000)], # s14 edge-type STEB on 2nd auroral form observed
        [dt.datetime(2022, 11, 20, 17, 24, 32, 600000), dt.datetime(2022, 11, 20, 17, 24, 33, 461000)], # s15 Undearath aurora
        [dt.datetime(2022, 11, 20, 17, 24, 45, 200000), dt.datetime(2022, 11, 20, 17, 24, 46, 000)], # s16 near the south edge of inverted V, almost undernearth it.
        [dt.datetime(2022, 11, 20, 17, 25, 23, 700000), dt.datetime(2022, 11, 20, 17, 25, 24, 300000)], # s17 Strong edge-type STEB near big aurora
        [dt.datetime(2022, 11, 20, 17, 26, 10, 000000), dt.datetime(2022, 11, 20, 17, 26, 11, 000000)], # s18 small STEB above aurora
    ]

    # the overhead noise on the dispersion features shouldn't be included. Here we have the energy indicies for the ~max energy the feature appears to be at when looking at pitch 0deg
    # NOTE: these are in Energy units (eV) where [min, Max]
    # This data was collected on the COUNTS graphs
    keyDispersionEnergyLimits = np.array([
        [28.22, 400],  # s1
        [28.22, 300],  # s2
        [28.22, 330],  # s3
        [28.22, 1000],  # s4
        [28.22, 1200],  # s5
        [114, 900],  # s6
        [83, 618],  # s7
        [28.22, 460],  # s8
        [28.22, 618],  # s9
        [28.22, 722],  # s10
        [28.22, 330],  # s11
        [28.22, 2000],  # s12
        [28.22, 1000],  # s13
        [28.22, 1000],  # s14
        [28.22, 1100],  # s15
        [28.22, 1000],  # s16
        [28.22, 600],  # s17
        [28.22, 360]  # s18
        ])

    # For the results of the study we need to know the EXACT deltaE and deltaT that the STEBS had, not just
    # limits useful for isolating data. this is hard to do by programming, so data was collected by hand
    # NOTE: Energy data is in eV
    keyDispersionDeltaE = np.array([
        [28.22, 340],  # s1
        [28.22, 256],  # s2
        [28.22, 293],  # s3
        [80, 750],  # s4
        [28.22, 845],  # s5
        [114, 630],  # s6
        [83, 288],  # s7
        [28.22, 460],  # s8
        [28.22, 540],  # s9
        [28.22, 540],  # s10
        [28.22, 245],  # s11
        [155, 990],  # s12
        [28.22, 250],  # s13
        [28.22, 470],  # s14
        [85, 714],  # s15
        [28.22, 115],  # s16
        [40, 290],  # s17
        [28.22, 250]  # s18
    ]
    )

    keyDispersionDeltaT = np.array([
        [dt.datetime(2022, 11, 20, 17, 24, 56, 205000), dt.datetime(2022, 11, 20, 17, 24, 57, 610000)],  # s1
        [dt.datetime(2022, 11, 20, 17, 24, 57, 710000), dt.datetime(2022, 11, 20, 17, 24, 59, 160000)],  # s2
        [dt.datetime(2022, 11, 20, 17, 24, 59, 55000), dt.datetime(2022, 11, 20, 17, 25, 0, 10000)],  # s3
        # [dt.datetime(2022, 11, 20, 17, 25, 00, 60000), dt.datetime(2022, 11, 20, 17, 25, 00, 552000)], # s4  right before Inverted V
        [dt.datetime(2022, 11, 20, 17, 25, 00, 60000), dt.datetime(2022, 11, 20, 17, 25, 0, 700000)], # s4  right before Inverted V
        [dt.datetime(2022, 11, 20, 17, 25, 00, 563000), dt.datetime(2022, 11, 20, 17, 25, 1, 312000)], # s5 inside inverted V's left edge
        [dt.datetime(2022, 11, 20, 17, 25, 1, 10000), dt.datetime(2022, 11, 20, 17, 25, 1, 462000)], # s6 under V, it's faint
        [dt.datetime(2022, 11, 20, 17, 25, 4, 8000), dt.datetime(2022, 11, 20, 17, 25, 4, 412000)], # s7 in choatic region
        [dt.datetime(2022, 11, 20, 17, 25, 5, 514000), dt.datetime(2022, 11, 20, 17, 25, 5, 870000)], # s8 VERTICAL BAR that doesn't look like it's in an inverted V
        [dt.datetime(2022, 11, 20, 17, 25, 5, 906000), dt.datetime(2022, 11, 20, 17, 25, 6, 506000)], # s9 Last dispersive feature in the chaotic region
        [dt.datetime(2022, 11, 20, 17, 25, 6, 858000), dt.datetime(2022, 11, 20, 17, 25, 7, 362000)], # s10 Very Strong Wave-like feature on north end of Inverted V
        [dt.datetime(2022, 11, 20, 17, 25, 7, 714000), dt.datetime(2022, 11, 20, 17, 25, 8, 310000)], # s11 Smaller feature, similar to s1
        [dt.datetime(2022, 11, 20, 17, 24, 18, 37000), dt.datetime(2022, 11, 20, 17, 24, 18, 666000)], # s12 first STEB observed. Dispersion underneath inverted V
        [dt.datetime(2022, 11, 20, 17, 24, 21, 152000), dt.datetime(2022, 11, 20, 17, 24, 22, 000)], # s13 STEB northward of first aurora observed. Faint
        [dt.datetime(2022, 11, 20, 17, 24, 31, 164000), dt.datetime(2022, 11, 20, 17, 24, 32, 00000)], # s14 edge-type STEB on 2nd auroral form observed
        [dt.datetime(2022, 11, 20, 17, 24, 32, 616000), dt.datetime(2022, 11, 20, 17, 24, 33, 356000)], # s15 Undearath aurora
        [dt.datetime(2022, 11, 20, 17, 24, 45, 558000), dt.datetime(2022, 11, 20, 17, 24, 45, 800000)], # s16 near the south edge of inverted V, almost undernearth it.
        [dt.datetime(2022, 11, 20, 17, 25, 23, 758000), dt.datetime(2022, 11, 20, 17, 25, 24, 166000)], # s17 Strong edge-type STEB near big aurora
        [dt.datetime(2022, 11, 20, 17, 26, 10, 262000), dt.datetime(2022, 11, 20, 17, 26, 10, 860000)], # s18 small STEB above aurora
    ]
    )

    # type options: (1) Auroral (2) isolated (3) edgetype
    keyDispersionType = [
        'isolated',  # s1
        'isolated',  # s2
        'isolated',  # s3
        'edgetype',  # s4
        'Auroral',  # s5
        'Auroral',  # s6
        'edgetype',  # s7
        'isolated',  # s8
        'isolated',  # s9
        'isolated',  # s10
        'isolated',  # s11
        'Auroral',  # s12
        'isolated',  # s13
        'edgetype',  # s14
        'Auroral',  # s15
        'edgetype',  # s16
        'edgetype',  # s17
        'edgetype'  # s18
    ]

    ##############################
    # --- CLEAN THE SIGNATURES ---
    ##############################

    # I need to isolate the dispersion signatures so I can T.O.F. fit them. To do this,
    # have to manually isolate them. Lets create a special set of functions for each dispersion to do just this.

    additionalRemoval = True

    @staticmethod
    def cleanS1(inputData,Energy,Time): # cleans up the corners

        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove upper Right stuff (COARSE)
        upper = True
        EnergyStartPoint, TimeStart, TimeEnd = 600, 0.8, 1.6
        newData = diagonalRemove(newData, EnergyStartPoint, TimeStart, TimeEnd, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS2(inputData,Energy,Time): # remove data in top right
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove upper Right stuff (COARSE)
        upper = True
        EnergyStartPoint, TimeStart, TimeEnd = 80, 1.1, 1.6
        newData = diagonalRemove(newData, EnergyStartPoint, TimeStart, TimeEnd, Energy, Time, upper)


        # # --- DIAGONAL REMOVE ---
        # remove upper Right stuff (COARSE)
        upper = True
        EnergyStartPoint, TimeStart, TimeEnd = 70, 0.9, 1.6
        newData = diagonalRemove(newData, EnergyStartPoint, TimeStart, TimeEnd, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS3(inputData, Energy, Time):  # doesn't need anything
        newData = inputData

        # --- BOX REMOVE ---
        # remove lower left
        EnergyMin, EnergyMax, TimeMin, TimeMax = 236, 10000, -0.1, 0.6
        newData = boxRemove(newData, EnergyMin, EnergyMax, TimeMin, TimeMax,Energy,Time)

        # --- DIAGONAL REMOVE ---
        # remove lower left stuff (COARSE)
        upper = False
        EnergyStartPoint, TimeStart, TimeEnd = 155, -0.05, 0.2
        newData = diagonalRemove(newData, EnergyStartPoint, TimeStart, TimeEnd, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left stuff (COARSE)
        upper = False
        EnergyStartPoint, TimeStart, TimeEnd = 80, -0.05, 0.401
        newData = diagonalRemove(newData, EnergyStartPoint, TimeStart, TimeEnd, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove upper Right stuff (COARSE)
        upper = True
        EnergyStartPoint, TimeStart, TimeEnd = 85, 0.6, 1.6
        newData = diagonalRemove(newData, EnergyStartPoint, TimeStart, TimeEnd, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS4(inputData, Energy, Time):  # doesn't need anything
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove bottom left stuff (COARSE)
        upper = False
        EnergyStartPoint, TimeStart, TimeEnd = 150, -0.1, 0.47
        newData = diagonalRemove(newData, EnergyStartPoint, TimeStart, TimeEnd, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove bottom left stuff (COARSE)
        upper = True
        EnergyStartPoint, TimeStart, TimeEnd = 150, 0.4, 1
        newData = diagonalRemove(newData, EnergyStartPoint, TimeStart, TimeEnd, Energy, Time, upper)


        # # --- DIAGONAL REMOVE ---
        # # remove bottom left stuff (COARSE)
        # upper = True
        # EnergyStartPoint, TimeStart, TimeEnd = 116, 0.45, 1
        # newData = diagonalRemove(newData, EnergyStartPoint, TimeStart, TimeEnd, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove bottom left stuff (COARSE)
        # upper = False
        # EnergyStartPoint, TimeStart, TimeEnd = 55, 0.5, 1
        # newData = diagonalRemove(newData, EnergyStartPoint, TimeStart, TimeEnd, Energy, Time, upper)


        return newData

    @staticmethod
    def cleanS5(inputData, Energy, Time):  # doesn't need anything
        newData = inputData

        if dispersionAttributes.additionalRemoval:
            # --- DIAGONAL REMOVE ---
            # remove lower left (COARSE)
            upper = False
            EnergyStartPoint, TimeStart, TimeEnd = 210, -0.1, 0.25
            newData = diagonalRemove(newData, EnergyStartPoint, TimeStart, TimeEnd, Energy, Time, upper)


            # --- DIAGONAL REMOVE ---
            upper = True
            EnergyStartPoint, TimeStart, TimeEnd = 250, 0.4, 1.5
            newData = diagonalRemove(newData, EnergyStartPoint, TimeStart, TimeEnd, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS6(inputData, Energy, Time):
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 540, -0.1, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 240, 0.15, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # # --- DIAGONAL REMOVE ---
        # # remove lower left (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = 520, 0.3, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- BOX REMOVE ---
        EnergyMin, EnergyMax, TimeMin, TimeMax = 0, 1000, 0.62, 2
        newData = boxRemove(newData, EnergyMin, EnergyMax, TimeMin, TimeMax, Energy, Time)

        return newData

    @staticmethod
    def cleanS7(inputData, Energy, Time):
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 250, -0.1, 0.8
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove upper right (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = 400, 0.4, 1.6
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove upper right (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = 180, 0.55, 1.6
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS8(inputData, Energy, Time):
        newData = inputData

        # # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 180, -0.1, 0.2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)
        #
        # # --- DIAGONAL REMOVE ---
        # remove upper right (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = 450, 0.2, 1
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)
        #
        upper = False
        diagonalCounter, TimeMin, TimeMax = 60, 0.14, 0.25
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)
        #
        upper = True
        diagonalCounter, TimeMin, TimeMax = 615, 0.0, 0.7
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = True
        diagonalCounter, TimeMin, TimeMax = 280, 0.17, 0.7
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = False
        diagonalCounter, TimeMin, TimeMax = 70, 0.13, 0.25
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = True
        diagonalCounter, TimeMin, TimeMax = 180, 0.24, 0.6
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = False
        diagonalCounter, TimeMin, TimeMax = 335, -0.1, 0.05
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData


    @staticmethod
    def cleanS9(inputData, Energy, Time):
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 280, -0.1, 0.05
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 112, -0.1, 0.4
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # remove lower left (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = 390, 0.35, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)


        return newData

    @staticmethod
    def cleanS10(inputData, Energy, Time):
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 240, -0.1, 0.15
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 60, 0.1, 0.5
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = 280, 0.5, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = 94, 0.646, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS11(inputData, Energy, Time):
        newData = inputData

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 180, -0.1, 0.3
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 60, 0.26, 0.5
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = 120,  0.6, 1
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        upper = True
        diagonalCounter, TimeMin, TimeMax = 380, 0.14, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS12(inputData, Energy, Time):
        newData = inputData



        EnergyMin, EnergyMax, TimeMin, TimeMax = 0, 2000, -0.1, 0.35
        newData = boxRemove(newData, EnergyMin, EnergyMax, TimeMin, TimeMax, Energy, Time)

        EnergyMin, EnergyMax, TimeMin, TimeMax = 0, 2000, 0.7, 2
        newData = boxRemove(newData, EnergyMin, EnergyMax, TimeMin, TimeMax, Energy, Time)

        EnergyMin, EnergyMax, TimeMin, TimeMax = 1404, 10000, 0.0, 2
        newData = boxRemove(newData, EnergyMin, EnergyMax, TimeMin, TimeMax, Energy, Time)

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = True
        diagonalCounter, TimeMin, TimeMax = 900, 0.44, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # --- DIAGONAL REMOVE ---
        # remove lower left (COARSE)
        upper = False
        diagonalCounter, TimeMin, TimeMax = 350, 0.3455, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        # upper = False
        # diagonalCounter, TimeMin, TimeMax = 500, 0.34, 1
        # newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)
        #
        # upper = True
        # diagonalCounter, TimeMin, TimeMax = 1000, 0.4,2
        # newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)


        return newData

    @staticmethod
    def cleanS13(inputData, Energy, Time):
        newData = inputData

        upper = False
        diagonalCounter, TimeMin, TimeMax = 350, -0.1, 0.65
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = True
        diagonalCounter, TimeMin, TimeMax = 400, 0.48, 1.5
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS14(inputData, Energy, Time):
        newData = inputData

        upper = False
        diagonalCounter, TimeMin, TimeMax = 240, -0.1, 0.4
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = True
        diagonalCounter, TimeMin, TimeMax = 175, 0.78, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = False
        diagonalCounter, TimeMin, TimeMax = 972, 0.1, 0.55
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = False
        diagonalCounter, TimeMin, TimeMax = 95, 0.54, 0.855
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = True
        diagonalCounter, TimeMin, TimeMax = 96, 0.948, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = True
        diagonalCounter, TimeMin, TimeMax = 1182, 0.447, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        EnergyMin, EnergyMax, TimeMin, TimeMax = 0, 10000, -0.1, 0.43
        newData = boxRemove(newData, EnergyMin, EnergyMax, TimeMin, TimeMax, Energy, Time)

        return newData

    @staticmethod
    def cleanS15(inputData, Energy, Time):
        newData = inputData

        upper = False
        diagonalCounter, TimeMin, TimeMax = 400, -0.1, 0.24
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = False
        diagonalCounter, TimeMin, TimeMax = 110, 0.4, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = True
        diagonalCounter, TimeMin, TimeMax = 400, 0.45, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        EnergyMin, EnergyMax, TimeMin, TimeMax = 0, 60, -0.1, 2
        newData = boxRemove(newData, EnergyMin, EnergyMax, TimeMin, TimeMax, Energy, Time)

        EnergyMin, EnergyMax, TimeMin, TimeMax = 0, 10000, -0.1, 0.21
        newData = boxRemove(newData, EnergyMin, EnergyMax, TimeMin, TimeMax, Energy, Time)

        return newData

    @staticmethod
    def cleanS16(inputData, Energy, Time):
        newData = inputData

        upper = True
        diagonalCounter, TimeMin, TimeMax = 550, 0.2, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = True
        diagonalCounter, TimeMin, TimeMax = 250, 0.38, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = True
        diagonalCounter, TimeMin, TimeMax = 130, 0.56, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = False
        diagonalCounter, TimeMin, TimeMax = 150, -0.1,0.45
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = False
        diagonalCounter, TimeMin, TimeMax = 200, 0.09, 0.3
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS17(inputData, Energy, Time):
        newData = inputData

        upper = True
        diagonalCounter, TimeMin, TimeMax = 400, 0.25, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData

    @staticmethod
    def cleanS18(inputData, Energy, Time):
        newData = inputData

        upper = True
        diagonalCounter, TimeMin, TimeMax = 250, 0.6, 2
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        upper = False
        diagonalCounter, TimeMin, TimeMax = 120, -0.1, 0.4
        newData = diagonalRemove(newData, diagonalCounter, TimeMin, TimeMax, Energy, Time, upper)

        return newData


    # Store the functions in a datavariable
    isolationFunctions = {'s1': cleanS1,
                          's2': cleanS2,
                          's3': cleanS3,
                          's4': cleanS4,
                          's5': cleanS5,
                          's6': cleanS6,
                          's7': cleanS7,
                          's8': cleanS8,
                          's9': cleanS9,
                          's10': cleanS10,
                          's11': cleanS11,
                          's12': cleanS12,
                          's13': cleanS13,
                          's14': cleanS14,
                          's15': cleanS15,
                          's16': cleanS16,
                          's17': cleanS17,
                          's18': cleanS18}