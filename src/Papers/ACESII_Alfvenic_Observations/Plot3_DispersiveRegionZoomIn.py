# --- Plot1_AllSky.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Stack Plot detailing field-aligned
# particle data along with electric and magnetic signatures


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import spaceToolsLib as stl
from src.my_imports import *
from scipy.signal import spectrogram
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
print(stl.color.UNDERLINE + f'Plot3_DispersiveregionZoomIn' + stl.color.END)
import matplotlib.pyplot as plt

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
dpi = 200


# --- Cbar ---
cbarMin, cbarMax = 1E6, 1E9
cbar_Tick_LabelSize = 14
my_cmap = stl.apl_rainbow_black0_cmap()
specCbarMin, specCbarMax = 1E-2, 1E1
spectrogramCmap = stl.blue_green_white_yellow_red_cmap()
my_cmap.set_bad(color=(0,0,0))

# Plot toggles
Figure_width = 7.5 # in inches
Figure_height =5.5*2.2# in inches
Text_FontSize = 20
Label_FontSize = 14.5
Label_Padding = 5
Tick_FontSize = 11.5
Tick_Length = 5
Tick_Width = 2
Tick_Padding = 10
Plot_LineWidth = 1
Legend_FontSize = 13

# Physics Toggles
targetILat = [71.91, 72.03]
targetEpoch = [dt.datetime(2022, 11, 20, 17, 24, 52, 0000000), dt.datetime(2022, 11, 20, 17, 25, 9, 500000)]

EField_scale = 1000 # what to scale the deltaE field by
EEPAAslice_wPitch = 2 # 10deg pitch
PAengyLimits = [17, 40] # determines in the P.A. panel which energies to count

Spectrogram_Freqlimits = [0, 12]



# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---
stl.prgMsg('Loading Data')
targetVar = [targetEpoch,'Epoch']


# delta B
inputMagFiles_high = glob('C:\Data\ACESII\L3\deltaB\high\*Field_Aligned*')[0]
data_dict_mag_high = stl.loadDictFromFile(inputFilePath=inputMagFiles_high, targetVar=targetVar, wKeys_Reduce=['B_e', 'B_r', 'B_p', 'ILat', 'Epoch', 'Alt'])
inputMagFiles_low = glob('C:\Data\ACESII\L3\deltaB\low\*Field_Aligned*')[0]
data_dict_mag_low = stl.loadDictFromFile(inputFilePath=inputMagFiles_low, targetVar=targetVar, wKeys_Reduce=['B_e', 'B_r', 'B_p', 'ILat', 'Epoch', 'Alt'])

# delta E
inputEFIFiles_low = glob('C:\Data\ACESII\L3\deltaE\low\*Field_Aligned*')[0]
data_dict_Efield_low = stl.loadDictFromFile(inputFilePath=inputEFIFiles_low, targetVar=targetVar, wKeys_Reduce=['E_e', 'E_r', 'E_p', 'ILat', 'Epoch', 'Alt'])

data_dict_Efield_low['E_e'][0] = EField_scale*data_dict_Efield_low['E_e'][0]
data_dict_Efield_low['E_p'][0] = EField_scale*data_dict_Efield_low['E_p'][0]
data_dict_Efield_low['E_r'][0] = EField_scale*data_dict_Efield_low['E_r'][0]

# EEPAA Particle Data
inputEEPAA_low = glob('C:\Data\ACESII\L2\low\*eepaa_fullCal*')[0]
data_dict_eepaa_low = stl.loadDictFromFile(inputFilePath=inputEEPAA_low, targetVar=targetVar, wKeys_Reduce=['Differential_Energy_Flux', 'ILat', 'Epoch', 'Alt'])
inputEEPAA_high = glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')[0]
data_dict_eepaa_high = stl.loadDictFromFile(inputFilePath=inputEEPAA_high, targetVar=targetVar, wKeys_Reduce=['Differential_Energy_Flux', 'ILat', 'Epoch', 'Alt'])

magDicts = [data_dict_mag_high, data_dict_mag_low]
eepaaDicts = [data_dict_eepaa_high, data_dict_eepaa_low]

stl.Done(start_time)

############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################
stl.prgMsg('Making Plot')
# --- PLOT EVERYTHING ---
fig, ax = plt.subplots(9, height_ratios=[1, 1, 1, 1,  0.6,   1, 1, 1, 1 ])
fig.set_figwidth(Figure_width)
fig.set_figheight(Figure_height)

for wRocket in [4, 5]:

    idxAdjust = 0 if wRocket == 4 else 5

    magDict = magDicts[wRocket-4]
    eepaaDict = eepaaDicts[wRocket-4]

    # --- Calculate Spectrogram ---
    spectrogramData = magDict['B_e'][0]
    windowType, npersegN, scalingType = 'hann', 128, 'spectrum'  # spectrogram toggles
    overlap = int(npersegN * (1 / 2))  # hanning filter overlap
    f, t, Sxx = spectrogram(spectrogramData,
                            fs=128,
                            window=windowType,
                            nperseg=npersegN,  # note: if ==None default size is 256
                            noverlap=overlap,
                            scaling=scalingType)  # scaling = density or scaling = spectrum

    # - determine a new ILat variable for the spectrogram -
    # first determine the times of the spectrogram
    startTime = pycdf.lib.datetime_to_tt2000(magDict['Epoch'][0][0])
    specTempTimes = [pycdf.lib.tt2000_to_datetime(int(startTime + 1E9*tme)) for tme in t]
    specIndicies = [np.abs(magDict['Epoch'][0] - specTempTimes[k]).argmin() for k in range(len(specTempTimes))]
    specILats = magDict['ILat'][0][specIndicies]
    specAlts = magDict['Alt'][0][specIndicies]
    specTimes = magDict['Epoch'][0][specIndicies]

    # --- get the averaged Differential Energy FLux over all Pitch Angles (Energy vs Time) ---
    totalDirFlux_engy = np.array(eepaaDict['Differential_Energy_Flux'][0])
    totalDirFlux_engy[totalDirFlux_engy<0] = 0
    totalDirFlux_engy = totalDirFlux_engy[:,2:19,:]
    totalDirFlux_engy = np.sum(totalDirFlux_engy, axis=1).T/(21-4)

    # --- get the averaged Differential Energy FLux over all energies (Pitch vs Time) ---
    totalDirFlux_ptch = np.array(eepaaDict['Differential_Energy_Flux'][0])
    totalDirFlux_ptch[totalDirFlux_ptch < 0] = 0
    totalDirFlux_ptch = np.sum(totalDirFlux_ptch, axis=2).T / len(data_dict_eepaa_high['Energy'][0])


    # --- EEPAA Epoch vs Time Data ---
    cmap = ax[0+idxAdjust].pcolormesh(eepaaDict['Epoch'][0], eepaaDict['Energy'][0], totalDirFlux_engy, cmap=my_cmap, vmin=cbarMin, vmax=cbarMax, norm='log')
    ax[0+idxAdjust].set_ylabel(rf'Energy [eV]', fontsize=Label_FontSize,labelpad=Label_Padding)
    ax[0+idxAdjust].set_yscale('log')
    ax[0+idxAdjust].set_ylim(33, 1150)
    ax[0+idxAdjust].tick_params(axis='both', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)

    # --- EEPAA Data All Pitchs ---
    ax[1+idxAdjust].pcolormesh(eepaaDict['Epoch'][0], eepaaDict['Pitch_Angle'][0], totalDirFlux_ptch, cmap=my_cmap, vmin=cbarMin, vmax=cbarMax, norm='log')
    ax[1+idxAdjust].set_ylabel(f'P. A. [deg]', fontsize=Label_FontSize,labelpad=Label_Padding)
    ax[1+idxAdjust].set_ylim(0, 180)
    ax[1+idxAdjust].margins(0)
    ax[1+idxAdjust].tick_params(axis='both', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
    yticks = [0, 60, 120, 180]
    ax[1+idxAdjust].set_yticks(yticks)
    ax[1+idxAdjust].set_yticklabels([str(tick) for tick in yticks])

    if wRocket == 4:
        # HF EEPAA colorbar
        cax = fig.add_axes([0.91, 0.774, 0.02, 0.212])
    elif wRocket == 5:
        # LF EEPAA colorbar
        cax = fig.add_axes([0.91, 0.276, 0.02, 0.211])

    cbar = plt.colorbar(cmap, cax=cax)
    cbar.ax.minorticks_on()
    cbar.ax.tick_params(labelsize=cbar_Tick_LabelSize+4)

    # --- Be/Er ---
    ax[2+idxAdjust].plot(magDict['Epoch'][0], magDict['B_e'][0], color='blue', linewidth=Plot_LineWidth)
    ax[2+idxAdjust].set_ylabel('$\delta B_{e}$\n[nT]', fontsize=Label_FontSize, color='blue',labelpad=Label_Padding)
    ax[2+idxAdjust].tick_params(axis='y', colors='blue')
    ax[2+idxAdjust].set_ylim(-8, 8)
    ax[2+idxAdjust].tick_params(axis='both', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
    ax[2+idxAdjust].margins(0)

    # Er
    if wRocket == 4:
        axEr = ax[2+idxAdjust].twinx()
        axEr.set_ylabel('E-Field\nN/A', fontsize=Label_FontSize, color='red', rotation=-90, labelpad=Label_Padding+30)
        axEr.tick_params(axis='y', colors='red',labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
        axEr.set_ylim(-8, 8)
    elif wRocket == 5:
        axEr = ax[2+idxAdjust].twinx()
        axEr.plot(data_dict_Efield_low['Epoch'][0],data_dict_Efield_low['E_r'][0], color='red', linewidth=Plot_LineWidth)
        axEr.set_ylabel('$\delta E_{r}$\n[mV/m]', fontsize=Label_FontSize, color='red', rotation=-90,labelpad=Label_Padding+30)
        axEr.tick_params(axis='y', colors='red',labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
        axEr.set_ylim(-8, 8)
        axEr.margins(0)

    # --- Be Spectrogram ---
    cmap = ax[3+idxAdjust].pcolormesh(specTimes, f, Sxx, shading='nearest', vmin=specCbarMin, vmax=specCbarMax, cmap=spectrogramCmap, norm='log')
    ax[3+idxAdjust].set_ylabel('$\delta B_{e}$\n[Hz]', fontsize=Label_FontSize, labelpad=Label_Padding, color='blue')
    ax[3+idxAdjust].set_xlabel('time [UTC]\nAlt [km]\nILat [deg]', fontsize=Tick_FontSize, weight='bold')
    ax[3+idxAdjust].xaxis.set_label_coords(-0.09, -0.09)
    ax[3+idxAdjust].set_ylim(Spectrogram_Freqlimits[0], Spectrogram_Freqlimits[1])
    ax[3 + idxAdjust].tick_params(axis='y', colors='blue')

    if wRocket ==4 :
        # spectrogram colorbar HF
        cax = fig.add_axes([0.91, 0.558, 0.02, 0.103])
    elif wRocket == 5:
        # spectrogram colorbar LF
        cax = fig.add_axes([0.91, 0.06, 0.02, 0.103])

    cbar = plt.colorbar(cmap, cax=cax)
    cbar.ax.minorticks_on()
    cbar.ax.tick_params(labelsize=cbar_Tick_LabelSize)


    # spectrogram Ticks
    yticks = [0, 4, 8, 12]
    ax[3+idxAdjust].set_yticks(yticks)
    ax[3+idxAdjust].set_yticklabels([str(tick) for tick in yticks])
    xtickTimes = [dt.datetime(2022,11,20,17,24,53),
                dt.datetime(2022,11,20,17,24,55,500000),
                dt.datetime(2022,11,20,17,24,58),
                dt.datetime(2022,11,20,17,25,00,500000),
                dt.datetime(2022,11,20,17,25,3),
                dt.datetime(2022,11,20,17,25,5,500000),
                dt.datetime(2022,11,20,17,25,8)]

    xtick_indicies = np.array([np.abs(magDict['Epoch'][0] - tick).argmin() for tick in xtickTimes])
    ILat_ticks = [str(round(tick, 2)) for tick in magDict['ILat'][0][xtick_indicies]]
    Alt_ticks = [str(round(tick, 1)) for tick in magDict['Alt'][0][xtick_indicies]]
    time_ticks = [tick.strftime("%H:%M:%S") for tick in magDict['Epoch'][0][xtick_indicies]]
    tickLabels = [f'{time_ticks[k]}\n{Alt_ticks[k]}\n{ILat_ticks[k]}' for k in range(len(xtick_indicies))]
    ax[3+idxAdjust].set_xticks(xtickTimes)
    ax[3+idxAdjust].set_xticklabels(tickLabels)
    ax[3+idxAdjust].tick_params(axis='y', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
    ax[3+idxAdjust].tick_params(axis='x', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)


# turn off center axis
ax[4].axis('off')

# turn off x-ticks for non-Emphemeris plots
for i in range(9):
    if i not in [3,8]:
        ax[i].tick_params(axis='x', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width,labelbottom=False,bottom=False,top=False)

# adjust plot size
fig.subplots_adjust(left=0.13, bottom=0.06, right=0.89, top=0.985, wspace=None,hspace=0.05)  # remove the space between plots
fig.align_ylabels(ax[:])

# output the figure
plt.savefig(rf'C:\Users\cfelt\Desktop\Research\ACESII\Feltman2025_ACESII_Alfven_Observations\PLOTS\Plot3\Plot3_base.png', dpi=dpi)
stl.Done(start_time)


