# --- Plots6_pitchAnglePlots.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Recreation of DiffEFlux pitch angle Plots, focusing on
# a few particle signatures


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
from src.my_imports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import spaceToolsLib as stl
from src.Science.AlfvenSingatureAnalysis.Particles.dispersionAttributes import dispersionAttributes

print(stl.color.UNDERLINE + f'Plot6_pitchAnglePlots' + stl.color.END)

#################
# --- TOGGLES ---
#################
useDiffNFlux = False
useCounts = False

sliceEpochIndicies = {
    's1':[5934, 5940, 5946],
    's2':[5959, 5966, 5974],
    's3':[5987 - 3, 5990 - 1, 5995],
    's4':[6003, 6007, 6011],
    's5':[6014, 6018 + 1, 6021 + 3],
    's10':[6139, 6142, 6145]  # s10 The Big One on the poleward side of the aurora
}
dispersiveRegionTargetTime = [dt.datetime(2022,11,20,17,24,56,500000),
                              dt.datetime(2022,11,20,17,25,2,000000)]


wDispersions = np.array([2, 3, 4, 5])-1 # [s1, s2, s3, s4, etc] <-- Index
wPitch_Engy_vs_Time_STEBS = [2, 4] # the pitch angle index to plot for the Energy vs time plot
wPitch_Engy_vs_Time_InvertedV = [5, 11] # the pitch angle index to plot for the Energy vs time plot
Energy_yLimit = 2000

# plot toggles - Slices pitch angle ------------------
X_Velocity_limits, Y_Velocity_limit = [-0.5, 1.6], [-1.6, 1.6]
NoOfSlices = 3

######################
# --- PLOT TOGGLES ---
######################
plt.rcParams["font.family"] = "Arial"
figure_height = (18.5)
figure_width = (16)

dpi = 200

Label_FontSize = 24
Label_Padding = 0.5
velSpace_Label_FontSize = 25

Text_FontSize = 18
Title_FontSize = 18
Tick_FontSize = 25
Tick_Width = 3
Tick_Length = 7
Plot_Linewidth = 6.5
vertical_lineStyles = ['dashdotdotted', 'dashdot', 'dashdotdotted', 'dashdot']
cbar_Fontsize = 24

my_cmap = stl.apl_rainbow_black0_cmap()
my_cmap.set_extremes(bad='white', under='black')

if useDiffNFlux:
    cbarLow, cbarHigh = 1E4, 5E6
elif useCounts:
    cbarLow, cbarHigh = 1, 100
else:
    # cbarLow, cbarHigh = 3E6, 1E9
    cbarLow, cbarHigh = 5E6, 1E9

normVal = 'log'

# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---
stl.prgMsg('Loading Data')

# Attitude data (for geomagnetic lat/long info)
inputFiles_Attitude = [glob(r'C:\Data\ACESII\attitude\high\*Attitude_Solution*')[0], glob(r'C:\Data\ACESII\attitude\low\*Attitude_Solution*')[0]]
data_dict_attitude_high = stl.loadDictFromFile(inputFilePath=inputFiles_Attitude[0])
data_dict_attitude_low = stl.loadDictFromFile(inputFilePath=inputFiles_Attitude[1])

# EEPAA Particle Data
inputFiles_eepaa = [glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')[0], glob('C:\Data\ACESII\L2\low\*eepaa_fullCal*')[0]]
data_dict_eepaa_high = stl.loadDictFromFile(inputFilePath=inputFiles_eepaa[0])
data_dict_counts_high = stl.loadDictFromFile(inputFilePath='C:\Data\ACESII\L1\high\ACESII_36359_l1_eepaa_fullCal.cdf')

stl.Done(start_time)

# --- --- --- --- --- --- --
# --- Calc Vpara & Vperp ---
# --- --- --- --- --- --- --
stl.prgMsg('Calculating Vperp and Vparallel')
Energy = data_dict_eepaa_high['Energy'][0]
Pitch = data_dict_eepaa_high['Pitch_Angle'][0]
if useDiffNFlux:
    countsTemp = data_dict_eepaa_high['Differential_Number_Flux'][0][6000]
elif useCounts:
    countsTemp =data_dict_counts_high['eepaa'][0][6000]
else:
    countsTemp = data_dict_eepaa_high['Differential_Energy_Flux'][0][6000]
Vperp = deepcopy(countsTemp)
Vpara = deepcopy(countsTemp)

for ptch in range(len(countsTemp)):
    for engy in range(len(countsTemp[0])):
        Vmag = np.sqrt(2*stl.q0*Energy[engy]/(stl.m_e))
        Vperp[ptch][engy] = np.sin(np.radians(Pitch[ptch]))*Vmag
        Vpara[ptch][engy] = np.cos(np.radians(Pitch[ptch]))*Vmag

Vpara, Vperp = np.array(Vpara)/(10000*1000), np.array(Vperp)/(10000*1000)
stl.Done(start_time)


# define (in time) where the STEBS occur, this has been done already
dispersionTimes = dispersionAttributes.keyDispersionTimes

# the slices in time for each dispersion used
sliceTimes = {key:[data_dict_eepaa_high['Epoch'][0][val] for val in sliceEpochIndicies[key]] for key,val in sliceEpochIndicies.items()}


# --- --- --- --- ---
# --- GET THE DATA ---
# --- --- --- --- ---
IndexLow, IndexHigh = np.abs(data_dict_eepaa_high['Epoch'][0] - dispersiveRegionTargetTime[0]).argmin(), np.abs(data_dict_eepaa_high['Epoch'][0] - dispersiveRegionTargetTime[1]).argmin()
Epoch = stl.EpochTo_T0_Rocket(InputEpoch=data_dict_eepaa_high['Epoch'][0][IndexLow:IndexHigh], T0=data_dict_eepaa_high['Epoch'][0][0])
if useDiffNFlux:
    dataArray = np.array(data_dict_eepaa_high['Differential_Number_Flux'][0][IndexLow:IndexHigh])
elif useCounts:
    dataArray = np.array(data_dict_counts_high['eepaa'][0][IndexLow:IndexHigh])
else:
    dataArray = np.array(data_dict_eepaa_high['Differential_Energy_Flux'][0][IndexLow:IndexHigh])


############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

stl.prgMsg('Beginning Plot')
fig = plt.figure()
fig.set_size_inches(figure_width,figure_height)
gs0 = gridspec.GridSpec(4, 1, figure=fig, height_ratios=[1, 1, 0.5, 4], hspace=0.05)

ax00 = fig.add_subplot(gs0[0, :])
ax01 = fig.add_subplot(gs0[1, :], sharex=ax00)


###################
# --- TOP PLOTS ---
###################

# --- --- --- --- --- --- --- --- --- --- -
# --- EEPAA Pitch Slice (10deg to 40deg)---
# --- --- --- --- --- --- --- --- --- --- -


# Sum and average the STEB pitch bins
totalDirFlux_engy = np.array(dataArray)
# totalDirFlux_engy[totalDirFlux_engy<0] = 0 # THIS WAS HERE IN THE 2nd Revision
totalDirFlux_engy = totalDirFlux_engy[:, wPitch_Engy_vs_Time_STEBS[0]:wPitch_Engy_vs_Time_STEBS[1]+1, :]
summedData = np.sum(totalDirFlux_engy, axis=1).T/(wPitch_Engy_vs_Time_STEBS[1] - wPitch_Engy_vs_Time_STEBS[0] + 1)
summedData[summedData<0] = np.nan # do this to make the colorbar look right
summedData[summedData ==0] = 1 # do this to make the colorbar look right

eepaaPitchSlice = ax00.pcolormesh(Epoch, Energy, summedData, cmap=my_cmap, norm=normVal, vmin=cbarLow, vmax=cbarHigh)
ax00.set_yscale('log')
ax00.set_ylabel(rf'{data_dict_eepaa_high["Pitch_Angle"][0][wPitch_Engy_vs_Time_STEBS[0]]}$^\circ < \alpha <$  {data_dict_eepaa_high["Pitch_Angle"][0][wPitch_Engy_vs_Time_STEBS[1]]}$^\circ$' +"\nEnergy [eV]",fontsize=Label_FontSize, weight='bold')
ax00.set_ylim(28, Energy_yLimit)
ax00.tick_params(axis='y', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
ax00.tick_params(axis='y', which='minor', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length*0.6)
ax00.tick_params(axis='x', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length, labelbottom=False)
ax00.tick_params(axis='x', which='minor', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length*0.6, labelbottom=False)


# --- --- --- --- --- --- --- --- --- --- -
# --- EEPAA Pitch Slice (50deg to 130deg)---
# --- --- --- --- --- --- --- --- --- --- -

totalDirFlux_engy = np.array(dataArray)
totalDirFlux_engy = totalDirFlux_engy[:, wPitch_Engy_vs_Time_InvertedV[0]:wPitch_Engy_vs_Time_InvertedV[1]+1, :]
summedData = np.sum(totalDirFlux_engy, axis=1).T/(wPitch_Engy_vs_Time_InvertedV[1] - wPitch_Engy_vs_Time_InvertedV[0] + 1)
summedData[summedData<0] = np.nan # do this to make the colorbar look right
summedData[summedData ==0] = 1 # do this to make the colorbar look right

eepaaPitchSlice = ax01.pcolormesh(Epoch, Energy, summedData, cmap=my_cmap,norm=normVal, vmin=cbarLow, vmax=cbarHigh)

ax01.set_yscale('log')
ax01.set_ylabel(rf'{data_dict_eepaa_high["Pitch_Angle"][0][wPitch_Engy_vs_Time_InvertedV[0]]}$^\circ < \alpha <$  {data_dict_eepaa_high["Pitch_Angle"][0][wPitch_Engy_vs_Time_InvertedV[1]]}$^\circ$' +"\nEnergy [eV]",fontsize=Label_FontSize, weight='bold')
ax01.set_ylim(28, Energy_yLimit)

ax01.set_xlabel('Seconds Since 17:20:00 UTC [s]', fontsize=Label_FontSize, weight='bold',labelpad=Label_Padding)
ax01.tick_params(axis='y', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
ax01.tick_params(axis='y', which='minor', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length/2)
ax01.tick_params(axis='x', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
ax01.tick_params(axis='x', which='minor', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length/2)
ax01.minorticks_on()

# plot the black vertical lines
vertical_lineStyles = ['dotted','dashdot','dotted','dashdot']
# for k,disIdx in enumerate(wDispersions):
#     for i in range(NoOfSlices):
#         index = np.abs(data_dict_eepaa_high['Epoch'][0] - sliceTimes[f's{wDispersions[colIdx] + 1}'][rowIdx]).argmin()
#         timeTag = round(EpochTo_T0_Rocket(InputEpoch=[sliceTimes[f's{disIdx + 1}'][i]], T0=data_dict_eepaa_high['Epoch'][0][0])[0], 2)
#         ax01.axvline(x=timeTag, color='black', linewidth=lineWidth, linestyle=vertical_lineStyles[k],alpha=0.6)





######################
# --- BOTTOM PLOTS ---
######################

# DEFINE THE TOP PLOTS AND BOTTOM PLOTS
topLabels = [['(c)','S2'],['(d)','S3'],['(e)','S4'],['(f)','S5']]
gs01 = gs0[3].subgridspec(3, 4, hspace=0.08,wspace=0.08)
for rowIdx in range(NoOfSlices):
    for colIdx in range(len(wDispersions)):

        ax = fig.add_subplot(gs01[rowIdx, colIdx])

        # colored textbox
        timeTag = round(stl.EpochTo_T0_Rocket(InputEpoch=[sliceTimes[f's{wDispersions[colIdx] + 1}'][rowIdx]], T0=data_dict_eepaa_high['Epoch'][0][0])[0], 2)
        props = dict(boxstyle='round', facecolor='white', alpha=1,lw=4)
        ax.text(0.5, -1.5, f'$t_{rowIdx}$=' + f'{timeTag}$\pm$0.05 s', fontsize=Text_FontSize, weight='bold', color='black', bbox=props, ha='center')

        # dataToPlot
        if useDiffNFlux:
            dataArray_Slice = data_dict_eepaa_high['Differential_Number_Flux'][0][np.abs(data_dict_eepaa_high['Epoch'][0] - sliceTimes[f's{wDispersions[colIdx] + 1}'][rowIdx]).argmin()]
        elif useCounts:
            dataArray_Slice = data_dict_counts_high['eepaa'][0][np.abs(data_dict_counts_high['Epoch'][0] - sliceTimes[f's{wDispersions[colIdx] + 1}'][rowIdx]).argmin()]
        else:
            index = np.abs(data_dict_eepaa_high['Epoch'][0] - sliceTimes[f's{wDispersions[colIdx] + 1}'][rowIdx]).argmin()
            dataArray_Slice = data_dict_eepaa_high['Differential_Energy_Flux'][0][index-1:index+2]
            dataArray_Slice = np.array(dataArray_Slice)
            dataArray_Slice = np.sum(dataArray_Slice, axis=0) / (3)
            dataArray_Slice[dataArray_Slice == 0] = 1
            dataArray_Slice[dataArray_Slice < 0] = np.nan

            # demarcation lines on the spectrogram plot
            sliceTime_TSL = (pycdf.lib.datetime_to_tt2000(sliceTimes[f's{wDispersions[colIdx] + 1}'][rowIdx]) - pycdf.lib.datetime_to_tt2000(data_dict_eepaa_high['Epoch'][0][0]))/1E9
            ax00.axvline(x= sliceTime_TSL,color='white', linewidth=Plot_Linewidth, linestyle=vertical_lineStyles[colIdx], alpha=0.8)

            if rowIdx == 1:
                props = dict(boxstyle='round', facecolor='white', alpha=1, lw=4)
                ax00.text(sliceTime_TSL, Energy_yLimit-200, f'S{colIdx+2}', fontsize=Text_FontSize, weight='bold', color='black', bbox=props, ha='center')

        # if normVal == 'log':
            # dataArray_Slice[dataArray_Slice<0] = np.nan
            # dataArray_Slice[dataArray_Slice ==0] = 1

        ax.pcolormesh(Vperp, Vpara, dataArray_Slice, cmap=my_cmap, shading='nearest', norm=normVal, vmin=cbarLow, vmax=cbarHigh)
        ax.set_xlim(X_Velocity_limits[0], X_Velocity_limits[1])
        ax.set_ylim(Y_Velocity_limit[0], Y_Velocity_limit[1])
        ax.invert_yaxis()

        ax.tick_params(axis='y', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
        ax.tick_params(axis='y', which='minor', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length / 2)
        ax.tick_params(axis='x', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
        ax.tick_params(axis='x', which='minor', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length / 2)

        Eval = [0.55, 1.05, 1.545]
        for eval in Eval:
            xVals = [eval * np.sin(np.radians(ptch)) for ptch in [-15 + i * 10 for i in range(21)]]
            yVals = [eval * np.cos(np.radians(ptch)) for ptch in [-15 + i*10 for i in range(21)]]
            ax.plot(xVals, yVals, label=f'{eval}', color='white', linewidth=1.5,alpha=0.5)

            if colIdx == 0:
                ax.text(x=eval * np.sin(np.radians(130)),y=eval * np.cos(np.radians(130)),s=f'{eval}', fontsize=13)

        if colIdx == 0:
            ax.set_ylabel('V$_{\parallel}$ [10$^{4}$ km/s]', fontsize=velSpace_Label_FontSize , labelpad=Label_Padding+5, weight='bold')

        else:
            ax.set_yticklabels([])

        if rowIdx == NoOfSlices-1:
            ax.set_xlabel('V$_{\perp}$ [10$^{4}$ km/s]', fontsize=velSpace_Label_FontSize, labelpad=Label_Padding, weight='bold')
        else:
            ax.set_xticklabels([])

        if rowIdx ==0:
            ax.text(X_Velocity_limits[0], -1.9, topLabels[colIdx][0], color='black', weight='bold', fontsize=Text_FontSize+12)
            ax.text((X_Velocity_limits[0]+X_Velocity_limits[1])/2, -1.9, topLabels[colIdx][1], color='black', weight='bold', fontsize=Text_FontSize + 20,ha='center')

# add in the label for S1
props = dict(boxstyle='round', facecolor='white', alpha=1, lw=4)
sliceTime_TSL = (pycdf.lib.datetime_to_tt2000(sliceTimes[f's{1}'][rowIdx]) - pycdf.lib.datetime_to_tt2000(data_dict_eepaa_high['Epoch'][0][0]))/1E9
ax00.text(sliceTime_TSL, Energy_yLimit-200, f'S1', fontsize=Text_FontSize, weight='bold', color='black', bbox=props, ha='center')

##################
# --- Colorbar ---
##################
cax = fig.add_axes([0.88, 0.05, 0.025, 0.93])
cbar = plt.colorbar(eepaaPitchSlice, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=Tick_FontSize,length=Tick_Length)

cbar.ax.get_yaxis().labelpad = 17
for l in cbar.ax.yaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(cbar_Fontsize)

if useDiffNFlux:
    cbarLabel = '[1/cm$^{2}$-s-sr-eV]'
elif useCounts:
    cbarLabel = 'Counts'
else:
    cbarLabel = '[eV/cm$^{2}$-s-sr-eV]'

cbar.set_label(cbarLabel,fontsize=cbar_Fontsize+25,rotation=90,labelpad=-5)


# Add (a), (b), (c), etc labels
ax00.text(296.5, 0.75E3, '(a)',color='white',weight='bold',fontsize=Text_FontSize+12)
ax01.text(296.5, 0.75E3, '(b)',color='white',weight='bold',fontsize=Text_FontSize+12)
plt.subplots_adjust(left=0.1, bottom=0.05, right=0.86, top=0.98, wspace=None, hspace=None)
plt.savefig(rf'C:\Users\cfelt\Desktop\Research\ACESII\Feltman2025_ACESII_Alfven_Observations\PLOTS\Plot6\Plot6_pitchAngle_base.png', dpi=dpi)

