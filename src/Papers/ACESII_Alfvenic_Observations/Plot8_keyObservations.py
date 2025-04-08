# --- Plots4_pitchAnglePlots.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Recreation of DiffEFlux pitch angle Plots, focusing on
# a few particle signatures


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

from src.my_imports import *
import matplotlib.gridspec as gridspec
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import spaceToolsLib as stl
import matplotlib.pyplot as plt
print(stl.color.UNDERLINE + f'Plot8_keyObservations' + stl.color.END)

#################
# --- TOGGLES ---
#################
figure_height = (12)
# figure_width = (12) # OLD
figure_width = (14)

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

dispersiveRegionTargetTime = [dt.datetime(2022,11,20,17,24,55,900000),
                              dt.datetime(2022,11,20,17,25,2,000000)]
# --- plot COUNTS and ENERGY  toggles ---
diffNFlux_limit_Low, diffNFlux_limit_High = 5E3, 3E7
wSTEBtoPlot = [1, 2, 3, 4, 5] # STEB number NOT index. Don't -1
wPitchs_to_plot = [2, 3, 4, 5] # decide which pitch angles to get the peak energy for
countsMask = 5
fluxMask = 5E7
Energy_yLimit_Idx = 14 # ONLY consider energies above this index -> energies BELOW ~ 1345 eV


##########################
# --- --- --- --- --- ---
# --- LOADING THE DATA ---
# --- --- --- --- --- ---
##########################

stl.prgMsg('Loading Data')
targetVarName = 'Epoch'
targetVar = dispersiveRegionTargetTime

# EEPAA Particle Data
inputFiles_eepaa = [glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')[0], glob('C:\Data\ACESII\L2\low\*eepaa_fullCal*')[0]]
data_dict_eepaa_high = stl.loadDictFromFile(inputFilePath=inputFiles_eepaa[0], targetVar=[targetVar, targetVarName])

inputFiles_eepaa_counts = [glob('C:\Data\ACESII\L1\high\*eepaa_fullCal*')[0], glob('C:\Data\ACESII\L1\low\*eepaa_fullCal*')[0]]
data_dict_counts_high = stl.loadDictFromFile(inputFilePath=inputFiles_eepaa_counts[0], targetVar=[targetVar, targetVarName],wKeys_Reduce=['counts','Epoch'])
stl.Done(start_time)


##########################
# --- --- --- --- --- ---
# --- PREPARE THE DATA ---
# --- --- --- --- --- ---
##########################
# --- time ---
from src.Science.AlfvenSingatureAnalysis.Particles.dispersionAttributes import dispersionAttributes
STEBtimes = [dispersionAttributes.keyDispersionDeltaT[idx-1] for idx in wSTEBtoPlot]
LaunchDateTime = pycdf.lib.datetime_to_tt2000(dt.datetime(2022,11,20,17,20,00,000000))
STEBtimes_rkt = [[(pycdf.lib.datetime_to_tt2000(tme[0]) - LaunchDateTime)/1E9, (pycdf.lib.datetime_to_tt2000(tme[1]) - LaunchDateTime)/1E9] for tme in STEBtimes]
rktTime_counts = stl.EpochTo_T0_Rocket(InputEpoch=data_dict_counts_high['Epoch'][0], T0=LaunchDateTime)

# --- particles ---
Pitch = data_dict_counts_high['Pitch_Angle'][0]
Energy = data_dict_counts_high['Energy'][0]
counts = data_dict_counts_high['counts'][0]
diffNFlux = data_dict_eepaa_high['Differential_Number_Flux'][0]


############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

stl.prgMsg('Beginning Plot')
fig = plt.figure()
fig.set_size_inches(figure_width, figure_height)

gs0 = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[0.3, 0.7])

gs00 = gs0[0].subgridspec(1, 1)
axPeakE = fig.add_subplot(gs00[:])

gs01 = gs0[1].subgridspec(2, 3, hspace=0.1, wspace=0.1)
axS2 = fig.add_subplot(gs01[0, 0])
axS3 = fig.add_subplot(gs01[0, 1])
axS4 = fig.add_subplot(gs01[1, 0])
axS5 = fig.add_subplot(gs01[1, 1])
ax_inV1 = fig.add_subplot(gs01[0, 2])
ax_inV2 = fig.add_subplot(gs01[1, 2])
subAxes = [axS2, axS3, axS4, axS5, ax_inV1, ax_inV2]



# --- --- --- --- --- --- --- --- ---
# --- PEAK COUNTS AND ENERGY PLOTS ---
# --- --- --- --- --- --- --- --- ---
STEB_peak_energies = []

for idx, ptchVal in enumerate(wPitchs_to_plot):

    peakFlux = []
    peakEnergy = []

    for t, tme in enumerate(STEBtimes_rkt):

        # Isolate the dispersion Feature
        wDis = wSTEBtoPlot[t]
        wDispersion_key = f's{wDis}'
        lowCut, highCut = np.abs(data_dict_eepaa_high['Epoch'][0] - dispersionAttributes.keyDispersionDeltaT[wDis - 1][0]).argmin(), np.abs(data_dict_eepaa_high['Epoch'][0] - dispersionAttributes.keyDispersionDeltaT[wDis - 1][1]).argmin()
        Epoch_dis = deepcopy(stl.dateTimetoTT2000(data_dict_eepaa_high['Epoch'][0][lowCut:highCut + 1], inverse=False))
        eepaa_dis_pre = deepcopy(counts[lowCut:highCut + 1])
        eepaa_dis = deepcopy(np.array(dispersionAttributes.isolationFunctions[wDispersion_key](eepaa_dis_pre, Energy, Epoch_dis)))  # apply the isolation functions found in dispersionAttributes.py

        # prepare the data by removing fillvals and applying the mask
        STEBdata = deepcopy(eepaa_dis)
        STEBdata[STEBdata < countsMask] = 0  # set anything below 0 = 0

        # --- PEAK ENERGY (from COUNTS data) ---
        peakEval = 0
        for timeIdx in range(len(STEBdata)):

            engySubSlice = STEBdata[timeIdx, ptchVal, Energy_yLimit_Idx:] # Contains Thresholded data of either non-zeros or zeros
            peakE_Idx = next((i for i, x in enumerate(engySubSlice) if x), None) # finds the first non-zero element in the array (highest energy)

            if peakE_Idx == None: # if there's no data at all
                EpVal = 0
            else: # if there's at least one element
                EpVal = Energy[Energy_yLimit_Idx+peakE_Idx]

            if EpVal > peakEval: # see if that element is higher than my peakEVal
                peakEval = Energy[Energy_yLimit_Idx+peakE_Idx]

        peakEnergy.append(peakEval)

    # plot the results
    axPeakE.plot(wSTEBtoPlot, peakEnergy, color=plot_Colors[idx], label=rf'$\alpha = {Pitch[ptchVal]}^\circ$', marker='.',ms=plot_MarkerSize)
    STEB_peak_energies.append(peakEnergy)


STEB_peak_energies = np.array(STEB_peak_energies).T

axPeakE.set_xmargin(0)
axPeakE.set_ylabel('STEB Peak Energy\n [eV]',fontsize=labels_FontSize,weight='bold')
axPeakE.legend(fontsize=legend_FontSize, loc='upper left')
axPeakE.set_ylim(-50, 1200)
axPeakE.set_xlabel('STEB #',fontsize=labels_FontSize, weight = 'bold')
axPeakE.grid(alpha=0.7)
axPeakE.set_xlim(0.1,5.9)
axPeakE.tick_params(axis='y', which='major', labelsize=tick_LabelSize, width=tick_Width, length=tick_Length)
axPeakE.tick_params(axis='y', which='minor', labelsize=tick_LabelSize, width=tick_Width, length=tick_Length / 2)
axPeakE.tick_params(axis='x', which='major', labelsize=tick_LabelSize+5, width=tick_Width, length=tick_Length)
axPeakE.tick_params(axis='x', which='minor', labelsize=tick_LabelSize+5, width=tick_Width, length=tick_Length / 2)



#########################################
# --- Calculate Flux vs energy graphs ---
#########################################
# --- Calculate Flux vs energy graphs ---
# peakEnergies = [245.74, 725.16, 245.74, 987.89]
for t, tme in enumerate(STEBtimes[1:]):

    if t == 0:
        wDis = 2
        peakEnergies =STEB_peak_energies[1]
    elif t == 1:
        wDis = 4
        peakEnergies = STEB_peak_energies[3]
    elif t == 2:
        wDis = 3
        peakEnergies = STEB_peak_energies[2]
    elif t == 3:
        wDis = 5
        peakEnergies = STEB_peak_energies[4]

    wDispersion_key = f's{wDis}'
    lowCut, highCut = np.abs(data_dict_eepaa_high['Epoch'][0] - dispersionAttributes.keyDispersionDeltaT[wDis - 1][0]).argmin(), np.abs(data_dict_eepaa_high['Epoch'][0] - dispersionAttributes.keyDispersionDeltaT[wDis - 1][1]).argmin()
    Epoch_dis = deepcopy(stl.dateTimetoTT2000(data_dict_eepaa_high['Epoch'][0][lowCut:highCut + 1], inverse=False))
    eepaa_dis_pre = deepcopy(diffNFlux[lowCut:highCut + 1])
    eepaa_dis = deepcopy(np.array(dispersionAttributes.isolationFunctions[wDispersion_key](eepaa_dis_pre, Energy, Epoch_dis)))  # pply the isolation functions found in dispersionAttributes.py

    # prepare the data by removing fillvals and applying the mask
    STEBdata = deepcopy(eepaa_dis)

    for idx, ptchVal in enumerate(wPitchs_to_plot):

        STEBdata[STEBdata <= 0] = np.nan  # set anything <= 0 to a nan for the averaging
        engyIdx = np.abs(Energy - peakEnergies[idx]).argmin()
        ptchSlice = STEBdata[:, ptchVal, engyIdx:].T
        Energies = Energy[engyIdx:]
        fluxSum = np.array([np.nanmean(ptchSlice[engy]) for engy in range(len(Energies))])

        # if t == 1:
        #     if ptchVal in [2,3,4,5]:
        #         fluxSum[np.where(Energies>1000)] = 0
        #     else:
        #         fluxSum[np.where(Energies > peakEnergy)] = 0

        subAxes[t].plot(Energies, fluxSum, color=plot_Colors[idx], label=rf'$\alpha = {Pitch[ptchVal]}^\circ$', marker='.', ms=plot_MarkerSize-7)

    props = dict(boxstyle='round', facecolor='white', alpha=1)
    subAxes[t].text(500, 2E7, s=f'S{wDis}', fontsize=text_FontSize, weight='bold', color='black', bbox=props, ha='center')
    subAxes[t].set_yscale('log')
    subAxes[t].set_ylim(diffNFlux_limit_Low, diffNFlux_limit_High)
    subAxes[t].set_xlim(-10, 1000)
    subAxes[t].tick_params(axis='y', which='major', labelsize=tick_SubplotLabelSize+4, width=tick_Width, length=tick_Length)
    subAxes[t].tick_params(axis='y', which='minor', labelsize=tick_SubplotLabelSize, width=tick_Width, length=tick_Length / 2)
    subAxes[t].tick_params(axis='x', which='major', labelsize=tick_SubplotLabelSize, width=tick_Width, length=tick_Length)
    subAxes[t].tick_params(axis='x', which='minor', labelsize=tick_SubplotLabelSize, width=tick_Width, length=tick_Length / 2)
    subAxes[t].grid(alpha=0.7)

    # set the labels and such
    if t in [0]:
        # subAxes[t].legend(prop={'size': legend_SubAxes_FontSize},loc='upper right') # OLD
        subAxes[t].legend(prop={'size': legend_SubAxes_FontSize}, loc='lower right')
    if t in [0, 2]:
        subAxes[t].set_ylabel('Avg. Diff. N. Flux \n [cm$^{-2}$-sr$^{-1}$-s$^{-1}$-eV$^{-1}$]', fontsize=labels_subplot_fontsize, weight='bold')
    if t in [2, 3]:
        subAxes[t].set_xlabel('Energy [eV]', fontsize=labels_FontSize, weight='bold')
    if t in [0, 1]:
        subAxes[t].set_xticklabels([])
    if t in [1, 3]:
        subAxes[t].set_yticklabels([])



##############################
# --- PLOT THE INVERTED-Vs ---
##############################

# --- --- --- --- -
# The S4 Inverted-V
# --- --- --- --- -
axID = 4
lowCut, highCut = np.abs(data_dict_eepaa_high['Epoch'][0] - dt.datetime(2022,11,20,17,25,0,512000)).argmin(), np.abs(data_dict_eepaa_high['Epoch'][0] - dt.datetime(2022,11,20,17,25,0,612000)).argmin()
inV_flux = deepcopy(diffNFlux[lowCut:highCut + 1])
# inV_flux[inV_flux < 0] = np.nan  # set anything below 0 = 0
# peakEnergy = 322 # inverted-V peak enregy

inV_flux[inV_flux <= 0] = np.nan  # set anything below 0 = 0
peakEnergy = 1000 # inverted-V peak enregy
for idx, ptchVal in enumerate(wPitchs_to_plot):
    engyIdx = np.abs(Energy - peakEnergy).argmin()
    ptchSlice = inV_flux[:, ptchVal, engyIdx:].T
    Energies = Energy[engyIdx:]
    # fluxSum = [sum(ptchSlice[engy]) / len(ptchSlice[engy]) for engy in range(len(Energies))]
    fluxSum = np.array([np.nanmean(ptchSlice[engy]) for engy in range(len(Energies))])
    fluxSum[np.where(Energies>400)] = 0
    subAxes[axID].plot(Energies, fluxSum, color=plot_Colors[idx], label=rf'$\alpha = {Pitch[ptchVal]}^\circ$', marker='.', ms=plot_MarkerSize-7)

props = dict(boxstyle='round', facecolor='white', alpha=1)
subAxes[axID].text(500, 2.4E7, s=f'Inverted-V (T={round(rktTime_counts[lowCut+1],1)} $\pm$ 0.05 s)', fontsize=text_FontSize-7, weight='bold', color='black', bbox=props, ha='center')
subAxes[axID].set_yscale('log')
subAxes[axID].set_ylim(diffNFlux_limit_Low, diffNFlux_limit_High)
subAxes[axID].set_xlim(-10, 1000)

subAxes[axID].tick_params(axis='y', which='major', labelsize=tick_SubplotLabelSize+4, width=tick_Width, length=tick_Length)
subAxes[axID].tick_params(axis='y', which='minor', labelsize=tick_SubplotLabelSize, width=tick_Width, length=tick_Length / 2)
subAxes[axID].tick_params(axis='x', which='major', labelsize=tick_SubplotLabelSize, width=tick_Width, length=tick_Length)
subAxes[axID].tick_params(axis='x', which='minor', labelsize=tick_SubplotLabelSize, width=tick_Width, length=tick_Length / 2)
subAxes[axID].grid(alpha=0.7)
subAxes[axID].set_yticklabels([])
subAxes[axID].set_xticklabels([])

# --- --- --- --- --- --
# The Inverted-V near S5
# --- --- --- --- --- --
axID = 5
lowCut, highCut = np.abs(data_dict_eepaa_high['Epoch'][0] - dt.datetime(2022,11,20,17,25,1,362000)).argmin(), np.abs(data_dict_eepaa_high['Epoch'][0] - dt.datetime(2022,11,20,17,25,1,462000)).argmin()
inV_flux = deepcopy(diffNFlux[lowCut:highCut + 1])
# inV_flux[inV_flux < 0] = np.nan  # set anything below 0 = 0
# peakEnergy = 800

inV_flux[inV_flux <= 0] = np.nan  # set anything below 0 = 0
peakEnergy = 1000
for idx, ptchVal in enumerate(wPitchs_to_plot):
    engyIdx = np.abs(Energy - peakEnergy).argmin()
    ptchSlice = inV_flux[:, ptchVal, engyIdx:].T
    Energies = Energy[engyIdx:]
    fluxSum = [np.nanmean(ptchSlice[engy]) for engy in range(len(Energies))]
    subAxes[axID].plot(Energies, fluxSum, color=plot_Colors[idx], label=rf'$\alpha = {Pitch[ptchVal]}^\circ$', marker='.', ms=plot_MarkerSize-7)

props = dict(boxstyle='round', facecolor='white', alpha=1)
subAxes[axID].text(500, 2.4E7, s=f'Inverted-V (T={round(rktTime_counts[lowCut+2],1)} $\pm$ 0.05 s)', fontsize=text_FontSize-7, weight='bold', color='black', bbox=props, ha='center')
subAxes[axID].set_yscale('log')
subAxes[axID].set_ylim(diffNFlux_limit_Low, diffNFlux_limit_High)
subAxes[axID].set_xlim(-10, 1000)
subAxes[axID].tick_params(axis='y', which='major', labelsize=tick_SubplotLabelSize+4, width=tick_Width, length=tick_Length)
subAxes[axID].tick_params(axis='y', which='minor', labelsize=tick_SubplotLabelSize, width=tick_Width, length=tick_Length / 2)
subAxes[axID].tick_params(axis='x', which='major', labelsize=tick_SubplotLabelSize, width=tick_Width, length=tick_Length)
subAxes[axID].tick_params(axis='x', which='minor', labelsize=tick_SubplotLabelSize, width=tick_Width, length=tick_Length / 2)
subAxes[axID].grid(alpha=0.6)
subAxes[axID].set_yticklabels([])
subAxes[axID].set_xlabel('Energy [eV]', fontsize=labels_FontSize, weight='bold')



# output the figure
plt.subplots_adjust(left=0.1, bottom=0.055, right=0.98, top=0.99, wspace=None, hspace=0.2)
fileOutName = rf'C:\Users\cfelt\Desktop\Research\ACESII\Feltman2025_ACESII_Alfven_Observations\PLOTS\Plot8\Plot8_keyObservations_base.png'
plt.savefig(fileOutName, dpi=dpi)
stl.Done(start_time)

