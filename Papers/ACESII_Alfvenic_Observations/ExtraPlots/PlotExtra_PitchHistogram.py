# --- Plots9_PoyntingFlux.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
import matplotlib.pyplot as plt
from myImports import *
import matplotlib.gridspec as gridspec
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import math
import matplotlib as mpl
import spaceToolsLib as stl

print(color.UNDERLINE + f'Plot8_keyObservations' + color.END)

#################
# --- TOGGLES ---
#################
figure_height = (16)
figure_width = (12)

cmap = 'turbo'
cmap = stl.apl_rainbow_black0_cmap()


plot_LineWidth = 3
plot_textFontSize = 20
plot_MarkerSize = 20
plot_Colors = ['tab:blue', 'tab:red', 'tab:orange', 'tab:green', 'tab:purple', 'tab:olive', 'black']
text_FontSize = 25
title_FontSize = 14
labels_FontSize = 22
labels_subplot_fontsize = 19
legend_FontSize = 15
legend_SubAxes_FontSize = 15
tick_LabelSize = 15
tick_SubplotLabelSize = 15
tick_Width = 2
tick_Length = 4
cbar_FontSize = 15
dpi = 100

dispersiveRegionTargetTime_HF = [dt.datetime(2022,11,20,17,24,55,900000),
                              dt.datetime(2022,11,20,17,25,2,000000)]

# plot toggles - Show STEB itself ----------
cbarLow_counts, cbarHigh_counts = 1, 100
diffEFlux_limit_Low, diffEFlux_limit_High = 1E3, 1E7
wPitch_Engy_vs_Time = 2 # the pitch angle index to plot for the Energy vs time plot
Energy_yLimit = 1350

# plot toggles - Histogram ---------------------------
Histogram_countsthresh = 4
# consider only the pitch angles between -10 and 90
# [  0 1  2  3  4  5  6  7  8  9 10 ...]
# [-10 0 10 20 30 40 50 60 70 80 90 ...]
pitchAngleWidthAcceptance_lowerlimit = 2
pitchAngleWidthAcceptance_upperlimit = 10 +1 # NEEDS +1 for array indexing
tempCmap = plt.cm.turbo_r # define new colormap
cmaplist = [tempCmap(i) for i in range(tempCmap.N)] # extract all colors from the colormap
cmap_hist = mpl.colors.LinearSegmentedColormap.from_list('turboCustom',cmaplist,tempCmap.N) # create the new map
bounds = 10*np.array([-1.5+i for i in range(11)])
histNorm = mpl.colors.BoundaryNorm(bounds, tempCmap.N)



##########################
# --- --- --- --- --- ---
# --- LOADING THE DATA ---
# --- --- --- --- --- ---
##########################

prgMsg('Loading Data')
targetVarName = 'Epoch'
targetVar = dispersiveRegionTargetTime

# Attitude data (for geomagnetic lat/long info)
inputFiles_Attitude = [glob(r'C:\Data\ACESII\attitude\high\*Attitude_Solution*')[0], glob(r'C:\Data\ACESII\attitude\low\*Attitude_Solution*')[0]]
data_dict_attitude_high = loadDictFromFile(inputFilePath=inputFiles_Attitude[0])
data_dict_attitude_low = loadDictFromFile(inputFilePath=inputFiles_Attitude[1])

# Magnetometer Data
inputFile_B = 'C:\Data\ACESII\L2\high\ACESII_36359_l2_RingCore_ENU.cdf'  # get the B data
data_dict_B = loadDictFromFile(inputFile_B,targetVar=[targetVar, targetVarName])

inputFile_deltaB = glob('C:\Data\ACESII\L3\deltaB\high\ACESII_36359_RingCore_Field_Aligned_WL250_stitchedFlight.cdf')[0] # get the deltaB data
data_dict_deltaB = deepcopy(loadDictFromFile(inputFile_deltaB,targetVar=[targetVar, targetVarName]))

# Langmuir Data
inputFile_Langmuir = 'C:\Data\ACESII\L3\Langmuir\high\ACESII_36359_langmuir_fixed_LowPass_low0.3_high0.3.cdf'
data_dict_langmuir = deepcopy(loadDictFromFile(inputFile_Langmuir,targetVar=[targetVar, targetVarName],wKeys_Load=['ni', 'Epoch', 'ILat']))
indexVals = [np.abs(data_dict_langmuir['Epoch'][0] - tme).argmin() for i,tme in enumerate(data_dict_B['Epoch'][0])]
data_dict_langmuir['ni'][0] = deepcopy(data_dict_langmuir['ni'][0][indexVals])
data_dict_langmuir['Epoch'][0] = deepcopy(data_dict_langmuir['Epoch'][0][indexVals])

# EISCAT Data (# Up-sample the EISCAT data (it's fine since the EISCAT variables are very slowly varying))
inputFile_EISCAT = 'C:\Data\ACESII\science\EISCAT_ACESII_Slice\high\ACESII_36359_EISCAT_Tromso_rktSlice.cdf'
data_dict_EISCAT = deepcopy(loadDictFromFile(inputFile_EISCAT, targetVar=[targetVar, targetVarName],wKeys_Load=['Ion_Comp', 'Op_Comp','ILat','Epoch']))
data_dict_EISCAT_interp = InterpolateDataDict(InputDataDict=data_dict_EISCAT,InputEpochArray=data_dict_EISCAT['Epoch'][0],targetEpochArray=data_dict_B['Epoch'][0],wKeys=[])
data_dict_EISCAT = deepcopy(data_dict_EISCAT_interp)

Done(start_time)


# --- --- --- --- --- --- --- ---
# --- MEDIAN PITCH ANGLE PLOT ---
# --- --- --- --- --- --- --- ---
X, Y = np.meshgrid(rktTime_counts, Energy)
Z = deepcopy(X)

# populate the new data
for tme in range(len(rktTime_counts)):
    for engy in range(len(Energy)):

        # get the data across all pitch angles here
        pitchData = counts[tme, :, engy]

        # Find errors in data and eliminate them from the median calculation
        for k, val in enumerate(pitchData):
            if val <= Histogram_countsthresh: # theshold away small COUNT values
                pitchData[k] = 0
            elif val > 1E14: # for extremely large outlier values
                pitchData[k] = 0
            elif math.isnan(val): # if is a nan
                pitchData[k] = 0

        pitchData_useThis = pitchData[pitchAngleWidthAcceptance_lowerlimit:pitchAngleWidthAcceptance_upperlimit]
        medianVal_pitchVal = np.nan

        halfDataVal = sum(pitchData_useThis) / 2
        for h in range(len(pitchData_useThis)):
            if halfDataVal - sum(pitchData_useThis[:h + 1]) < 0:
                medianVal_pitchVal = Pitch[h+pitchAngleWidthAcceptance_lowerlimit]
                break

        Z[engy][tme] = medianVal_pitchVal

# adjust the plot

# cmapHist = axPitchHist.pcolormesh(X, Y, Z, cmap=cmap_hist, norm=histNorm)
cmapHist = axPitchHist.pcolormesh(X, Y, np.transpose(data_dict_eepaa_high['Differential_Energy_Flux'][0][:,2,:]), cmap=cmap, vmin=5E7,vmax=1E9 ,norm='log')
axPitchHist.set_yscale('log')
axPitchHist.set_ylim(Energy[-1], Energy_yLimit)
axPitchHist.set_ylabel('Energy  [eV]', fontsize=labels_FontSize)
axPitchHist.tick_params(axis='x', which='major', labelbottom=False)
axPitchHist.set_xmargin(0)

