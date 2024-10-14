# --- Plot1_AllSky.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Stack Plot detailing field-aligned
# particle data along with electric and magnetic signatures


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import spaceToolsLib as stl
from myImports import *
from scipy.signal import spectrogram
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
print(color.UNDERLINE + f'Plot2_Conjugacy' + color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
dpi = 800




# --- Cbar ---
# cbarMin, cbarMax = 5E6, 3E9
cbarMin, cbarMax = 2E7, 3E9
cbarTickLabelSize = 14
my_cmap = stl.apl_rainbow_black0_cmap()
my_cmap.set_bad(color=(1, 1, 1))

Escale = 1000 # what to scale the deltaE field by


# --- HF/LF General ---
plot_General = False
General_figure_width = 7.5 # in inches
General_figure_height =10# in inches
General_targetILat = [71.35, 72.85]
General_targetEpoch = [dt.datetime(2022,11,20,17,24,54,000000), dt.datetime(2022,11,20,17,25,11,00000)]
General_LabelFontSize = 12
General_TickFontSize = 13
General_PlotLineWidth = 0.5
General_EBlimits = 9
General_LabelPadding = 8
GeneralCmap = my_cmap

# --- HF/LF Dispersive Region ---
plot_Dispersive = True
Dispersive_figure_width = 7.5
Dispersive_figure_height = 5.5*(2.2)
Disp_targetILat = [71.91, 72.03]
Dis_targetEpoch = [dt.datetime(2022,11,20,17,24,53,750000), dt.datetime(2022,11,20,17,25,8,000000)]
Dispersive_wPitch = 2 # 10deg pitch
PAengyLimits = [17, 40] # determines in the P.A. panel which energies to count
Dispersive_LabelFontSize = 14.5
Dispersive_TickFontSize = 11.5
Dispersive_TickLength = 5
Dispersive_TickWidth = 2
Dispersive_LabelPadding = 5
DispersiveCmap = my_cmap
specCbarMin, specCbarMax = 1E-2, 1E1
DispersiveFreqlimits = [0, 12]
Disp_PlotLineWidth = 1
spectrogramCmap = stl.blue_green_white_yellow_red_cmap()

if plot_General:
    alignByILat = True # Align the Data via ILat
elif plot_Dispersive:
    alignByILat = False  # Align the Data via ILat
else:
    alignByILat = True  # Align the Data via ILat

# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---
prgMsg('Loading Data')
if plot_General:
    targetILat = General_targetILat
    targetEpoch = General_targetEpoch
    plot_Dispersive = False
elif plot_Dispersive:
    targetILat = Disp_targetILat
    targetEpoch = Dis_targetEpoch

if alignByILat:
    targetVar = [targetILat,'ILat']
else:
    targetVar = [targetEpoch,'Epoch']



rocketAttrs, b, c = ACES_mission_dicts()

# delta B
inputMagFiles_high = glob('C:\Data\ACESII\L3\deltaB\high\*Field_Aligned*')[0]
data_dict_mag_high = loadDictFromFile(inputFilePath=inputMagFiles_high, targetVar=targetVar, wKeys_Reduce=['B_e', 'B_r', 'B_p', 'ILat', 'Epoch', 'Alt'])
inputMagFiles_low = glob('C:\Data\ACESII\L3\deltaB\low\*Field_Aligned*')[0]
data_dict_mag_low = loadDictFromFile(inputFilePath=inputMagFiles_low, targetVar=targetVar, wKeys_Reduce=['B_e', 'B_r', 'B_p', 'ILat', 'Epoch', 'Alt'])

# delta E
inputEFIFiles_low = glob('C:\Data\ACESII\L3\deltaE\low\*Field_Aligned*')[0]
data_dict_Efield_low = loadDictFromFile(inputFilePath=inputEFIFiles_low, targetVar=targetVar, wKeys_Reduce=['E_e', 'E_r', 'E_p', 'ILat', 'Epoch', 'Alt'])

data_dict_Efield_low['E_e'][0] = Escale*data_dict_Efield_low['E_e'][0]
data_dict_Efield_low['E_p'][0] = Escale*data_dict_Efield_low['E_p'][0]
data_dict_Efield_low['E_r'][0] = Escale*data_dict_Efield_low['E_r'][0]

# EEPAA Particle Data
inputEEPAA_low = glob('C:\Data\ACESII\L2\low\*eepaa_fullCal*')[0]
data_dict_eepaa_low = loadDictFromFile(inputFilePath=inputEEPAA_low, targetVar=targetVar, wKeys_Reduce=['Differential_Energy_Flux', 'ILat', 'Epoch', 'Alt'])
inputEEPAA_high = glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')[0]
data_dict_eepaa_high = loadDictFromFile(inputFilePath=inputEEPAA_high, targetVar=targetVar, wKeys_Reduce=['Differential_Energy_Flux', 'ILat', 'Epoch', 'Alt'])


# LP Particle Data
inputLP_low = glob('C:\Data\ACESII\L3\Langmuir\low\*langmuir_fixed*')[1]
data_dict_LP_low = loadDictFromFile(inputFilePath=inputLP_low, targetVar=targetVar, wKeys_Reduce=['ni', 'ILat', 'Epoch'])
inputLP_high = glob('C:\Data\ACESII\L3\Langmuir\high\*langmuir_fixed*')[1]
data_dict_LP_high = loadDictFromFile(inputFilePath=inputLP_high, targetVar=targetVar, wKeys_Reduce=['ni', 'ILat', 'Epoch'])
Done(start_time)




############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

# --- PLOT EVERYTHING ---
fig, ax = plt.subplots(9, height_ratios=[2, 2, 1, 1,     0.6,   2, 2, 1, 1 ])
fig.set_figwidth(Dispersive_figure_width)
fig.set_figheight(Dispersive_figure_height)


for wRocket in [4, 5]:

    idxAdjust = 0 if wRocket == 4 else 5



    magDicts = [data_dict_mag_high, data_dict_mag_low]
    magDict = magDicts[wRocket-4]
    eepaaDicts = [data_dict_eepaa_high, data_dict_eepaa_low]
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

    # --- get the total Differential Energy FLux over all energies ---
    totalDirFlux = np.zeros(shape=(len(eepaaDict['Differential_Energy_Flux'][0]), len(eepaaDict['Pitch_Angle'][0])))
    for tme in range(len(eepaaDict['Epoch'][0])):
        for ptch in range(len(eepaaDict['Pitch_Angle'][0])):
            sumVal = 0
            for engy in range(PAengyLimits[0],PAengyLimits[1]+1):
                val = eepaaDict['Differential_Energy_Flux'][0][tme][ptch][engy]
                if val >0:
                    sumVal += val
            totalDirFlux[tme][ptch] = sumVal

    # --- PLOT EVERYTHING ---
    # fig, ax = plt.subplots(4, sharex=True, height_ratios=[2,2,1,1])
    # fig.set_figwidth(Dispersive_figure_width)
    # fig.set_figheight(Dispersive_figure_height)
    # fig.subplots_adjust(top=0.97, hspace=0.1)  # remove the space between plots

    # --- EEPAA Data 10deg ---
    cmap = ax[0+idxAdjust].pcolormesh(eepaaDict['Epoch'][0], eepaaDict['Energy'][0], eepaaDict['Differential_Energy_Flux'][0][:,Dispersive_wPitch,:].T, cmap=GeneralCmap, vmin=cbarMin, vmax=cbarMax, norm='log')
    ax[0+idxAdjust].set_ylabel(rf'P.A. = {eepaaDict["Pitch_Angle"][0][Dispersive_wPitch]}$^\circ$' + '\n Energy [eV]', fontsize=Dispersive_LabelFontSize,labelpad=Dispersive_LabelPadding)
    ax[0+idxAdjust].set_yscale('log')
    ax[0+idxAdjust].set_ylim(33, 1150)
    ax[0+idxAdjust].tick_params(axis='both', labelsize=Dispersive_TickFontSize, length=Dispersive_TickLength, width=Dispersive_TickWidth)

    # --- EEPAA Data All Pitch ---
    ax[1+idxAdjust].pcolormesh(eepaaDict['Epoch'][0], eepaaDict['Pitch_Angle'][0], totalDirFlux.T, cmap=GeneralCmap, vmin=cbarMin, vmax=cbarMax, norm='log')
    ax[1+idxAdjust].set_ylabel(f'{round(eepaaDict["Energy"][0][PAengyLimits[1]])}-{round(eepaaDict["Energy"][0][PAengyLimits[0]])} eV \nP. A. [deg]', fontsize=Dispersive_LabelFontSize,labelpad=Dispersive_LabelPadding)
    ax[1+idxAdjust].set_ylim(0, 181)
    ax[1+idxAdjust].margins(0)
    ax[1+idxAdjust].tick_params(axis='both', labelsize=Dispersive_TickFontSize, length=Dispersive_TickLength, width=Dispersive_TickWidth)
    yticks = [0, 60, 120, 180]
    ax[1+idxAdjust].set_yticks(yticks)
    ax[1+idxAdjust].set_yticklabels([str(tick) for tick in yticks])

    if wRocket == 4:
        # HF EEPAA colorbar
        cax = fig.add_axes([0.91, 0.71, 0.02, 0.281])
    elif wRocket == 5:
        # LF EEPAA colorbar
        cax = fig.add_axes([0.91, 0.215, 0.02, 0.28])
    cbar = plt.colorbar(cmap, cax=cax)
    cbar.ax.minorticks_on()
    cbar.ax.tick_params(labelsize=cbarTickLabelSize+4)

    # --- Be/Er ---
    ax[2+idxAdjust].plot(magDict['Epoch'][0], magDict['B_e'][0], color='blue', linewidth=Disp_PlotLineWidth)
    ax[2+idxAdjust].set_ylabel('$\delta B_{e}$\n[nT]', fontsize=Dispersive_LabelFontSize, color='blue',labelpad=Dispersive_LabelPadding+10)
    ax[2+idxAdjust].tick_params(axis='y', colors='blue')
    ax[2+idxAdjust].set_ylim(-8, 8)
    ax[2+idxAdjust].tick_params(axis='both', labelsize=Dispersive_TickFontSize, length=Dispersive_TickLength, width=Dispersive_TickWidth)
    ax[2+idxAdjust].margins(0)

    # Er
    if wRocket == 4:
        axEr = ax[2+idxAdjust].twinx()
        axEr.set_ylabel('E-Field\nN/A', fontsize=Dispersive_LabelFontSize, color='red', rotation=-90, labelpad=Dispersive_LabelPadding+25)
        axEr.tick_params(axis='y', colors='red',labelsize=Dispersive_TickFontSize, length=Dispersive_TickLength, width=Dispersive_TickWidth)
        axEr.set_ylim(-8, 8)
    elif wRocket == 5:
        axEr = ax[2+idxAdjust].twinx()
        axEr.plot(data_dict_Efield_low['Epoch'][0],data_dict_Efield_low['E_r'][0], color='red', linewidth=Disp_PlotLineWidth)
        axEr.set_ylabel('$\delta E_{r}$\n[mV/m]', fontsize=Dispersive_LabelFontSize, color='red', rotation=-90,labelpad=Dispersive_LabelPadding+25)
        axEr.tick_params(axis='y', colors='red',labelsize=Dispersive_TickFontSize, length=Dispersive_TickLength, width=Dispersive_TickWidth)
        axEr.set_ylim(-8, 8)
        axEr.margins(0)

    # --- Be Spectrogram ---
    cmap = ax[3+idxAdjust].pcolormesh(specTimes, f, Sxx, shading='nearest', vmin=specCbarMin, vmax=specCbarMax, cmap=spectrogramCmap, norm='log')
    ax[3+idxAdjust].set_ylabel('$\delta B_{e}$ Freq.\n[Hz]', fontsize=Dispersive_LabelFontSize, labelpad=Dispersive_LabelPadding+5)
    ax[3+idxAdjust].set_xlabel('time [UTC]\nAlt [km]\nILat [deg]', fontsize=Dispersive_TickFontSize, weight='bold')
    ax[3+idxAdjust].xaxis.set_label_coords(-0.09, -0.14)
    ax[3+idxAdjust].set_ylim(DispersiveFreqlimits[0], DispersiveFreqlimits[1])


    if wRocket ==4 :
        # spectrogram colorbar HF
        cax = fig.add_axes([0.91, 0.554, 0.02, 0.0682])
    elif wRocket == 5:
        # spectrogram colorbar LF
        cax = fig.add_axes([0.91, 0.06, 0.02, 0.0682])

    cbar = plt.colorbar(cmap, cax=cax)
    cbar.ax.minorticks_on()
    cbar.ax.tick_params(labelsize=cbarTickLabelSize)


    # spectrogram Ticks
    yticks = [0, 4, 8, 12]
    ax[3+idxAdjust].set_yticks(yticks)
    ax[3+idxAdjust].set_yticklabels([str(tick) for tick in yticks])
    xtickTimes = [dt.datetime(2022,11,20,17,24,55),
                dt.datetime(2022,11,20,17,24,57),
                dt.datetime(2022,11,20,17,24,59),
                dt.datetime(2022,11,20,17,25,1),
                dt.datetime(2022,11,20,17,25,3),
                dt.datetime(2022,11,20,17,25,5),
                  dt.datetime(2022,11,20,17,25,7)]

    xtick_indicies = np.array([np.abs(magDict['Epoch'][0] - tick).argmin() for tick in xtickTimes])
    ILat_ticks = [str(round(tick, 2)) for tick in magDict['ILat'][0][xtick_indicies]]
    Alt_ticks = [str(round(tick, 1)) for tick in magDict['Alt'][0][xtick_indicies]]
    time_ticks = [tick.strftime("%H:%M:%S") for tick in magDict['Epoch'][0][xtick_indicies]]
    tickLabels = [f'{time_ticks[k]}\n{Alt_ticks[k]}\n{ILat_ticks[k]}' for k in range(len(xtick_indicies))]
    ax[3+idxAdjust].set_xticks(xtickTimes)
    ax[3+idxAdjust].set_xticklabels(tickLabels)
    ax[3+idxAdjust].tick_params(axis='y', labelsize=Dispersive_TickFontSize, length=Dispersive_TickLength, width=Dispersive_TickWidth)
    ax[3+idxAdjust].tick_params(axis='x', labelsize=Dispersive_TickFontSize, length=Dispersive_TickLength, width=Dispersive_TickWidth)


# turn off center axis
ax[4].axis('off')

fig.subplots_adjust(top=0.97, hspace=0.1)  # remove the space between plots
fig.subplots_adjust(left=0.13, bottom=0.06, right=0.90, top=0.99, wspace=None,hspace=0.1)  # remove the space between plots

# output the figure
outputMod = 'HF' if wRocket==4 else 'LF'
timeSpace = 'space' if alignByILat else 'time'
plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot2\Plot2_dispersive_{timeSpace}.png', dpi=dpi)



