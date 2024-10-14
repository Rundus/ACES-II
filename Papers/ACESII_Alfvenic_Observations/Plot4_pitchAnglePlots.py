# --- Plots4_pitchAnglePlots.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Recreation of DiffEFlux pitch angle Plots, focusing on
# a few particle signatures


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import matplotlib.gridspec as gridspec
from myspaceToolsLib.time import EpochTo_T0_Rocket
from ACESII_code.Science.AlfvenSingatureAnalysis.Particles.dispersionAttributes import dispersionAttributes




print(color.UNDERLINE + f'Plot4_pitchAnglePlots' + color.END)

#################
# --- TOGGLES ---
#################
useDiffNFlux = False
# plot toggles - Overview -----------------------------
sliceEpochIndicies = {
    's1':[5934, 5940, 5946],
    's2':[5959, 5966, 5974],
    's3':[5987 - 3, 5990 - 1, 5995],
    # 's4':[6002 , 6005 , 6008 ],
    's4':[6003, 6007, 6011],
    's5':[6014, 6018 + 1, 6021 + 3],
    's10':[6139, 6142, 6145]  # s10 The Big One on the poleward side of the aurora
}
# dispersiveRegionTargetTime = [dt.datetime(2022,11,20,17,24,57,600000),
#                               dt.datetime(2022,11,20,17,25,2,750000)]
dispersiveRegionTargetTime = [dt.datetime(2022,11,20,17,24,56,500000),
                              dt.datetime(2022,11,20,17,25,2,000000)]

figure_height = (15)
figure_width = (12.5)

from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
cmap = apl_rainbow_black0_cmap()

labelsFontSize = 20
smallPlots_labelsFontSize = 17
labelPadding = -0.25
textFontSize = 18
titleFontSize = 18
tickFontSize = 18
tickWidth = 2
tickLength = 2.5
lineWidth = 4
cbarFont = 18
dpi = 200


# plot toggles - Show STEB itself ----------
if useDiffNFlux:
    cbarLow, cbarHigh = 1E5, 1E7
else:
    cbarLow, cbarHigh = 4E7, 1E9

wDispersions = np.array([2,3,4,5])-1 # [s1, s2, s3, s4, etc] <-- Index
wPitch_Engy_vs_Time = [0,1,2] # the pitch angle index to plot for the Energy vs time plot
Energy_yLimit = 1350

# plot toggles - Slices pitch angle ------------------
X_Velocity_limits, Y_Velocity_limit = [-0.5, 1.6], [-1.6, 1.6]
NoOfSlices = 3



# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---

prgMsg('Loading Data')

# Attitude data (for geomagnetic lat/long info)
inputFiles_Attitude = [glob(r'C:\Data\ACESII\attitude\high\*Attitude_Solution*')[0], glob(r'C:\Data\ACESII\attitude\low\*Attitude_Solution*')[0]]
data_dict_attitude_high = loadDictFromFile(inputFilePath=inputFiles_Attitude[0])
data_dict_attitude_low = loadDictFromFile(inputFilePath=inputFiles_Attitude[1])

# EEPAA Particle Data
inputFiles_eepaa = [glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')[0], glob('C:\Data\ACESII\L2\low\*eepaa_fullCal*')[0]]
data_dict_eepaa_high = loadDictFromFile(inputFilePath=inputFiles_eepaa[0])

Done(start_time)


# --- --- --- --- --- --- --
# --- Calc Vpara & Vperp ---
# --- --- --- --- --- --- --
prgMsg('Calculating Vperp and Vparallel')
Energy = data_dict_eepaa_high['Energy'][0]
Pitch = data_dict_eepaa_high['Pitch_Angle'][0]
if useDiffNFlux:
    countsTemp = data_dict_eepaa_high['Differential_Number_Flux'][0][6000]
else:
    countsTemp = data_dict_eepaa_high['Differential_Energy_Flux'][0][6000]
Vperp = deepcopy(countsTemp)
Vpara = deepcopy(countsTemp)

for ptch in range(len(countsTemp)):
    for engy in range(len(countsTemp[0])):
        Vmag = np.sqrt(2*q0*Energy[engy]/(m_e))
        Vperp[ptch][engy] = np.sin(np.radians(Pitch[ptch]))*Vmag
        Vpara[ptch][engy] = np.cos(np.radians(Pitch[ptch]))*Vmag

Vpara, Vperp = np.array(Vpara)/(10000*1000), np.array(Vperp)/(10000*1000)
Done(start_time)


# define (in time) where the STEBS occur, this has been done already

dispersionTimes = dispersionAttributes.keyDispersionTimes
# dispersionTimes = dispersionAttributes.keyDispersionDeltaT

# the slices in time for each dispersion used
sliceTimes = {key:[data_dict_eepaa_high['Epoch'][0][val] for val in sliceEpochIndicies[key]] for key,val in sliceEpochIndicies.items()}


# --- --- --- --- ---
# --- GET THE DATA ---
# --- --- --- --- ---
IndexLow, IndexHigh = np.abs(data_dict_eepaa_high['Epoch'][0] - dispersiveRegionTargetTime[0]).argmin(), np.abs(data_dict_eepaa_high['Epoch'][0] - dispersiveRegionTargetTime[1]).argmin()
Epoch = EpochTo_T0_Rocket(InputEpoch=data_dict_eepaa_high['Epoch'][0][IndexLow:IndexHigh], T0=data_dict_eepaa_high['Epoch'][0][0])
if useDiffNFlux:
    dataArray = np.array(data_dict_eepaa_high['Differential_Number_Flux'][0][IndexLow:IndexHigh])
else:
    dataArray = np.array(data_dict_eepaa_high['Differential_Energy_Flux'][0][IndexLow:IndexHigh])




############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

prgMsg('Beginning Plot')
fig = plt.figure()
fig.set_size_inches(figure_width,figure_height)

gs0 = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1,4], hspace=0.18)

# DEFINE THE TOP PLOTS AND BOTTOM PLOTS
gs01 = gs0[1].subgridspec(3, 4, hspace=0.08,wspace=0.08)



###################
# --- TOP PLOTS ---
###################
# --- EEPAA Pitch Slice ---
ax00 = fig.add_subplot(gs0[0,:])

# Sum and average the lowest 3 bins
pitchSlices = [dataArray[:,slice,:] for slice in wPitch_Engy_vs_Time]

summedData = np.zeros(shape=(len(pitchSlices[0]),len(pitchSlices[0][0])))

for tme in range(len(pitchSlices[0])):
    for engy in range(len(pitchSlices[0][0])):
        values = np.array([slice[tme][engy] for slice in pitchSlices])
        rmvFILL = values[values>0]
        if len(rmvFILL) == 0:
            summedData[tme][engy]= 0
        else:
            summedData[tme][engy] = sum(rmvFILL)/len(rmvFILL)


eepaaPitchSlice = ax00.pcolormesh(Epoch, Energy, summedData.T, cmap=cmap,norm='log', vmin=cbarLow, vmax=cbarHigh)


ax00.set_yscale('log')
ax00.set_ylabel('Energy [eV]',fontsize=labelsFontSize, weight='bold')
ax00.set_ylim(28, Energy_yLimit)

ax00.set_xlabel('Time Since Launch [s]', fontsize=labelsFontSize, weight='bold',labelpad=labelPadding)
ax00.tick_params(axis='y', which='major', labelsize=tickFontSize, width=tickWidth, length=tickLength)
ax00.tick_params(axis='y', which='minor', labelsize=tickFontSize, width=tickWidth, length=tickLength/2)
ax00.tick_params(axis='x', which='major', labelsize=tickFontSize, width=tickWidth, length=tickLength)
ax00.tick_params(axis='x', which='minor', labelsize=tickFontSize, width=tickWidth, length=tickLength/2)

# plot the black vertical lines
vertical_lineStyles = ['dotted','dashdot','dotted','dashdot']

# for k,disIdx in enumerate(wDispersions):
#     for i in range(NoOfSlices):
#         timeTag = round(EpochTo_T0_Rocket(InputEpoch=[sliceTimes[f's{disIdx + 1}'][i]], T0=data_dict_eepaa_high['Epoch'][0][0])[0], 2)
#         ax00.axvline(x=timeTag, color='black', linewidth=lineWidth, linestyle=vertical_lineStyles[k],alpha=0.6)




######################
# --- BOTTOM PLOTS ---
######################
for rowIdx in range(NoOfSlices):
    for colIdx in range(len(wDispersions)):

        ax = fig.add_subplot(gs01[rowIdx,colIdx])

        # colored textbox
        timeTag = round(EpochTo_T0_Rocket(InputEpoch=[sliceTimes[f's{wDispersions[colIdx] + 1}'][rowIdx]], T0=data_dict_eepaa_high['Epoch'][0][0])[0], 2)
        props = dict(boxstyle='round', facecolor='white', alpha=1)
        ax.text(0.5, -1.3, f'$t_{rowIdx}$=' + f'{timeTag} s', fontsize=textFontSize, weight='bold', color='black', bbox=props, ha='center')

        # dataToPlot
        if useDiffNFlux:
            dataArray_Slice = data_dict_eepaa_high['Differential_Number_Flux'][0][np.abs(data_dict_eepaa_high['Epoch'][0] - sliceTimes[f's{wDispersions[colIdx] + 1}'][rowIdx]).argmin()]
        else:
            dataArray_Slice = data_dict_eepaa_high['Differential_Energy_Flux'][0][np.abs(data_dict_eepaa_high['Epoch'][0] - sliceTimes[f's{wDispersions[colIdx] + 1}'][rowIdx]).argmin()]



        # Set the background black by turning all 0 values into 1's (just for display purposes)
        # for tme in range(len(dataArray_Slice)):
        #     for engy in range(len(dataArray_Slice[0])):
        #         if not dataArray_Slice[tme][engy] >= 1:
        #             dataArray_Slice[tme][engy] = 1

        ax.pcolormesh(Vperp, Vpara, dataArray_Slice, cmap=cmap, shading='nearest', norm='log', vmin=cbarLow, vmax=cbarHigh)
        ax.set_xlim(X_Velocity_limits[0], X_Velocity_limits[1])
        ax.set_ylim(Y_Velocity_limit[0], Y_Velocity_limit[1])
        ax.invert_yaxis()

        ax.tick_params(axis='y', which='major', labelsize=tickFontSize, width=tickWidth, length=tickLength)
        ax.tick_params(axis='y', which='minor', labelsize=tickFontSize, width=tickWidth, length=tickLength / 2)
        ax.tick_params(axis='x', which='major', labelsize=tickFontSize, width=tickWidth, length=tickLength)
        ax.tick_params(axis='x', which='minor', labelsize=tickFontSize, width=tickWidth, length=tickLength / 2)

        Eval = [0.55, 1.05, 1.545]
        for eval in Eval:
            xVals = [eval * np.sin(np.radians(ptch)) for ptch in [-15 + i * 10 for i in range(21)]]
            yVals = [eval * np.cos(np.radians(ptch)) for ptch in [-15 + i*10 for i in range(21)]]
            ax.plot(xVals, yVals, label=f'{eval}', color='black', linewidth=1.5,alpha=0.5)

            if colIdx == 0:
                ax.text(x=eval * np.sin(np.radians(130)),y=eval * np.cos(np.radians(130)),s=f'{eval}', fontsize=13)

        if colIdx == 0:
            ax.set_ylabel('V$_{\parallel}$ [10$^{4}$ km/s]', fontsize=smallPlots_labelsFontSize , labelpad=labelPadding + 5, weight='bold')

        else:
            ax.set_yticklabels([])

        if rowIdx == NoOfSlices-1:
            ax.set_xlabel('V$_{\perp}$ [10$^{4}$ km/s]', fontsize=smallPlots_labelsFontSize, labelpad=labelPadding + 5, weight='bold')
        else:
            ax.set_xticklabels([])




##################
# --- Colorbar ---
##################
cax = fig.add_axes([0.93, 0.05, 0.025, 0.94])
cbar = plt.colorbar(eepaaPitchSlice, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=tickFontSize,length=tickLength)

cbar.ax.get_yaxis().labelpad = 17
for l in cbar.ax.yaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(cbarFont)


plt.subplots_adjust(left=0.07, bottom=0.05, right=0.92, top=0.99, wspace=None, hspace=None)

plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot4\\Plot4_pitchAngle_base.png', dpi=dpi)

