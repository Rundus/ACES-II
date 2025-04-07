# --- Plots6_pitchAnglePlots.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Recreation of DiffEFlux pitch angle Plots, focusing on
# a few particle signatures

# TODO: sort the white line on the spectrogram plots


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
plt.rcParams["font.family"] = "Arial"

print(stl.color.UNDERLINE + f'Plot7p5_DetermineSourceDistributionCurves' + stl.color.END)

#################
# --- TOGGLES ---
#################
wPitchsToFit = [2, 3, 4, 5]

######################
# --- PLOT TOGGLES ---
######################
figure_height = 6*(3)
figure_width = (16)

dpi = 200

Label_FontSize = 24
Label_Padding = 3
velSpace_Label_FontSize = 25

Legend_fontSize = 14
LimitsMargin = 0.2
Text_FontSize = 18
Title_FontSize = 18
Tick_FontSize = 20
Tick_Width = 2
Tick_Length = 7
Plot_Linewidth = 3
Errorbar_capsize = 5
Scatter_markersize = 18

vertical_lineStyles = ['dashdotdotted','dashdot','dashdotdotted','dashdot']
cbar_Fontsize = 24

my_cmap = stl.apl_rainbow_black0_cmap()
my_cmap.set_bad(color=(0, 0, 0))
my_cmap.set_under('black')
# cbarMin,cbarMax = 5E6, 1E9 # for diffEFlux
cbarMin,cbarMax = 1E5, 1E7


# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---

stl.prgMsg('Loading Data')

# Dispersive region fit data
inputFile = rf'C:\Users\cfelt\Desktop\Research\ACESII\Feltman2025_ACESII_Alfven_Observations\PLOTS\Plot7p5\primaryBeam_fitting_parameters_DispersiveRegion.cdf'
data_dict_dispersiveFitData = stl.loadDictFromFile(inputFilePath=inputFile)

# primary inverted-V data
inputFile = rf'C:\Users\cfelt\Desktop\Research\ACESII\Feltman2025_ACESII_Alfven_Observations\PLOTS\Plot7p5\primaryBeam_fitting_parameters_PrimaryInvertedV.cdf'
data_dict_invertedVFitData = stl.loadDictFromFile(inputFilePath=inputFile)

data_dicts = [data_dict_dispersiveFitData, data_dict_invertedVFitData]

# DiffNFlux data
inputFile = 'C:\Data\ACESII\L2\high\ACESII_36359_l2_eepaa_fullCal.cdf'
data_dict_diffNFlux = stl.loadDictFromFile(inputFilePath=inputFile)

# DistributionFunc Data
inputFile = 'C:\Data\ACESII\L3\DistFunc\high\ACESII_36359_distFunc_eepaa_fullCal.cdf'
data_dict_dist = stl.loadDictFromFile(inputFilePath=inputFile)


# Attitude Data
inputFile = glob(r'C:\Data\ACESII\attitude\high\*.cdf')[0]
data_dict_attitude = stl.loadDictFromFile(inputFilePath=inputFile)
stl.Done(start_time)


############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

stl.prgMsg('Beginning Plot')
fig = plt.figure()
fig.set_size_inches(figure_width, figure_height)

gs0 = gridspec.GridSpec(4, 1, figure=fig, height_ratios=[1.25/6, 3.25/6, 1/12, 1.5/6])
gs00 = gs0[0].subgridspec(1, 2)
gs01 = gs0[1].subgridspec(4, 2)
gs02 = gs0[3].subgridspec(1, 2)
gsFill = gs0[2].subgridspec(1, 1)

avg_fitParams = []

axesStored = []

for tmeRangeIdx in range(2):

    # get the data dict set
    data_dict = data_dicts[tmeRangeIdx]

    # pull out the individual data
    Epoch = data_dict['timestamp_fitData'][0].flatten()
    Density = data_dict['n'][0][:,2].flatten()
    Temperature = data_dict['Te'][0][:,2].flatten()
    Potential = data_dict['V0'][0][:,2].flatten()
    ChiSquare = data_dict['ChiSquare'][0][:,2].flatten()

    # --- Do some initial processing on the data ---

    ### find the average electrostatic potential for each time##
    avgPeakE = np.sum(data_dict['V0'][0],axis=1)/len(data_dict['V0'][0])
    avgPeakE_Epoch = Epoch[0]

    # --- [1] make the spectrogram from 0 to 40 deg ---
    axSpec = fig.add_subplot(gs00[0, tmeRangeIdx])
    titleVal = 'Dispersive Region Inverted-V' if tmeRangeIdx == 0 else 'Primary Inverted-V'
    axSpec.set_title(titleVal, fontsize=Title_FontSize+10,weight='bold')
    diffFlux_avg = np.array(data_dict_diffNFlux['Differential_Number_Flux'][0])
    diffFlux_avg[diffFlux_avg < 0] = 0
    diffFlux_avg = diffFlux_avg[:, wPitchsToFit, :]
    diffFlux_avg = np.sum(diffFlux_avg, axis=1) / (len(wPitchsToFit))

    # reduce the spectrogram to only the fit region
    lowIdx, highIdx = np.abs(data_dict_diffNFlux['Epoch'][0] - Epoch[0]).argmin(), np.abs(data_dict_diffNFlux['Epoch'][0] - Epoch[-1]).argmin()
    cmap = axSpec.pcolormesh(data_dict_diffNFlux['Epoch'][0][lowIdx:highIdx+1], data_dict_diffNFlux['Energy'][0], diffFlux_avg[lowIdx:highIdx+1].T, cmap=my_cmap, vmin=cbarMin, vmax=cbarMax, norm='log') # plot spectrogram
    # axSpec.scatter(Epoch, Potential, marker='o',color='white',s=Scatter_markersize)
    sortedLists = np.array(sorted(zip(Epoch, Potential))).T
    goodIndicies = np.where(sortedLists[1] > 0 )
    lp_Epoch = sortedLists[0][goodIndicies]
    lp_Potential = sortedLists[1][goodIndicies]
    # sort the data to make the line plot nicer
    if tmeRangeIdx == 0:
        axSpec.scatter(lp_Epoch, lp_Potential,s=Scatter_markersize, color='white')
    else:
        axSpec.plot(lp_Epoch, lp_Potential, color='white', linewidth=Plot_Linewidth)

    axSpec.set_ylabel(rf'Energy [eV]', fontsize=Label_FontSize,labelpad=Label_Padding) if tmeRangeIdx == 0 else None
    axSpec.set_yscale('log')
    axSpec.set_ylim(28.22, 2500)

    # --- [2] Make the Chi Square Plot ---
    axChi = fig.add_subplot(gs01[0, tmeRangeIdx],sharex=axSpec)
    axChi.scatter(Epoch, ChiSquare, marker='o',color='black',s=Scatter_markersize)
    axChi.set_ylim((1-LimitsMargin)*0.1, (1+LimitsMargin)*ChiSquare.max())
    axChi.set_ylabel(r'$\chi ^{2}_{\nu}$', fontsize=Label_FontSize,labelpad=Label_Padding) if tmeRangeIdx == 0 else None
    axChi.set_yscale('log')
    avg_Chi = round(sum([val for val in ChiSquare if val > 0]) / len(Potential), 1)
    axChi.axhline(avg_Chi, color='red', label=rf'$\chi^2_\nu$ (Avg) = {avg_Chi}')  # plot the average value
    axChi.legend(fontsize=Legend_fontSize,loc='upper right')

    # --- [3] Make the Density Plot ---
    axDensity = fig.add_subplot(gs01[1, tmeRangeIdx],sharex=axSpec)
    axDensity.scatter(Epoch, Density, marker='o',color='black',s=Scatter_markersize)
    axDensity.set_ylim((1-LimitsMargin)*0, (1+LimitsMargin)*Density.max())
    axDensity.set_ylabel('$n_{0}$  [cm$^{-3}$]', fontsize=Label_FontSize,labelpad=Label_Padding) if tmeRangeIdx == 0 else None
    avg_Density = round(sum([val for val in Density if val > 0]) / len(Potential), 1)
    axDensity.axhline(avg_Density, color='red', label=rf'$n_{0}$ (Avg) = {avg_Density}')  # plot the average value
    axDensity.legend(fontsize=Legend_fontSize,loc='upper right')

    # --- [4] Make the Parallel Potential Plot ---
    axV0 = fig.add_subplot(gs01[2, tmeRangeIdx],sharex=axSpec)
    axV0.scatter(Epoch, Potential, marker='o',color='black',s=Scatter_markersize)
    axV0.set_ylabel('$V_{0}$ [eV]', fontsize=Label_FontSize,labelpad=Label_Padding) if tmeRangeIdx == 0 else None
    axV0.set_yscale('log')
    avg_V0 = round(sum([val for val in Potential if val > 0]) / len(Potential), 1)
    axV0.axhline(avg_V0, color='red', label=rf'$V_0$ (Avg) = {avg_V0}')  # plot the average value
    axV0.legend(fontsize=Legend_fontSize,loc='upper right')
    axV0.set_ylim(90, 2000)

    # --- [4] Make the Temperature Plot ---
    axTemp = fig.add_subplot(gs01[3, tmeRangeIdx],sharex=axSpec)
    axTemp.scatter(Epoch, Temperature, marker='o',color='black',s=Scatter_markersize)
    axTemp.set_ylim((1-LimitsMargin)*10, (1+LimitsMargin)*Temperature.max())
    axTemp.set_ylabel('$T_{e}$ [eV]', fontsize=Label_FontSize,labelpad=Label_Padding) if tmeRangeIdx == 0 else None
    axTemp.set_yscale('log')
    avg_Temp = round(sum([val for val in Temperature if val > 0]) / len(Potential), 1)
    axTemp.axhline(avg_Temp, color='red', label=rf'$T_e$ (Avg) = {avg_Temp}')  # plot the average value
    axTemp.legend(fontsize=Legend_fontSize,loc='upper right')
    axTemp.set_xlabel('time [UTC]\nAlt [km]\nILat [deg]', fontsize=15, weight='bold')
    axTemp.xaxis.set_label_coords(-0.07, -0.09)

    avg_fitParams.append([avg_Density,avg_Temp,avg_V0])

    # Adjust the xticks
    xtickTimes = [[
                   dt.datetime(2022, 11, 20, 17, 25, 1, 300000),
                  dt.datetime(2022, 11, 20, 17, 25, 1, 700000),
                    dt.datetime(2022, 11, 20, 17, 25, 2, 100000),
                   dt.datetime(2022, 11, 20, 17, 25, 2, 500000)],

                  [dt.datetime(2022, 11, 20, 17, 25, 26,500000),
                   dt.datetime(2022, 11, 20, 17, 25, 31,500000),
                   dt.datetime(2022, 11, 20, 17, 25, 36,500000),
                   dt.datetime(2022, 11, 20, 17, 25, 41,500000),
                   dt.datetime(2022, 11, 20, 17, 25, 46,500000)],
                  ]

    xtick_indicies = np.array([np.abs(data_dict_attitude['Epoch'][0] - tick).argmin() for tick in xtickTimes[tmeRangeIdx]])
    ILat_ticks = [str(round(tick, 2)) for tick in data_dict_attitude['ILat'][0][xtick_indicies]]
    Alt_ticks = [str(round(tick/1000, 1)) for tick in data_dict_attitude['Alt'][0][xtick_indicies]]
    time_ticks = [tick.strftime("%M:%S.") + str(round(tick.microsecond/1000)) for tick in xtickTimes[tmeRangeIdx]]
    tickLabels = [f'{time_ticks[k]}\n{Alt_ticks[k]}\n{ILat_ticks[k]}' for k in range(len(xtick_indicies))]
    axTemp.set_xticks(xtickTimes[tmeRangeIdx])
    axTemp.set_xticklabels(tickLabels)

    # Make general adjustments
    axes = [axSpec, axChi, axDensity, axV0, axTemp]
    fig.align_ylabels(axes[:])
    for i, ax in enumerate(axes):
        ax.minorticks_on()
        ax.tick_params(axis='y', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
        ax.tick_params(axis='x', labelsize=Tick_FontSize-5, length=Tick_Length, width=Tick_Width)
        if i > 0:
            ax.grid(True,which='Major', alpha=0.5)
        if i < 4:
            ax.tick_params(axis='x', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width,labelbottom=False,bottom=False,top=False)
        if tmeRangeIdx:
            ax.tick_params(axis='y', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width, labelleft=False, bottom=False, top=False)

    # store the axes for later
    axesStored.append([axSpec,axChi,axDensity,axV0,axTemp])

cax = fig.add_axes([0.915, 0.814-0.015, 0.022, 0.171])
cbar = plt.colorbar(cmap, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbar_Fontsize)
cbar.set_label('[cm$^{-2}$s$^{-1}$sr$^{-1}$eV$^{-1}$]',fontsize=cbar_Fontsize, weight='bold')

##############################################
# --- --- --- --- --- --- --- --- --- --- ---
# --- MAKE THE DISTRIBUTIONS & EXAMPLE FIT ---
# --- --- --- --- --- --- --- --- --- --- ---
##############################################
wDataset = 3

exampleDataset = [
    [dt.datetime(2022, 11, 20, 17, 25, 1, 162210), 0], # 10deg
    [dt.datetime(2022, 11, 20, 17, 25, 1, 462208), 0], # 10deg
    [dt.datetime(2022, 11, 20, 17, 25, 2, 362214), 0],# 10deg
    [dt.datetime(2022, 11, 20, 17, 25, 2, 862215), 0],# 10deg
    [dt.datetime(2022, 11, 20, 17, 25, 1, 212210), 1]# 20deg
]

targetTime = exampleDataset[wDataset][0]
targetidx = np.abs(data_dict_dist['Epoch'][0]-targetTime).argmin()
wExamplePtch = exampleDataset[wDataset][1]

# plot text label on first plot
axesStored[0][0].axvline(targetTime,color='white',alpha=0.7,linewidth=Plot_Linewidth,linestyle="--")
axesStored[0][0].text(targetTime- dt.timedelta(0,0,30000),1500,s=f'Example Fit (T={targetTime.strftime("%H:%M:%S.%f")} UTC)',color='white',fontsize=Text_FontSize,ha='right',va='center',weight='bold')
axesStored[0][0].scatter(targetTime, 1500, s=2*Scatter_markersize, color='white')

##########################
# --- Example Fit Plot ---
##########################
axExampleFit = fig.add_subplot(gs02[0, 0])

# --- DEFINE THE NOISE LEVEL ---
countNoiseLevel = 4
def generateNoiseLevel(energyData, countNoiseLevel):
    count_interval = ACESII.ESA_count_interval*1E-6
    geo_factor = ACESII.ESA_geometric_factor_TRICEII[0][0]
    deadtime = ACESII.ESA_deadtime[0]


    # --- DEFINE THE NOISE LEVEL ---
    diffNFlux_NoiseCount = np.zeros(shape=(len(energyData)))

    for idx, engy in enumerate(energyData):
        deltaT = (count_interval) - (countNoiseLevel * deadtime)
        diffNFlux_NoiseCount[idx] = (countNoiseLevel) / (geo_factor * deltaT * engy)

    return diffNFlux_NoiseCount


stl.Done(start_time)

# --- Get Example fitted data ---
Energy = data_dict_diffNFlux['Energy'][0]
Pitch = data_dict_diffNFlux['Pitch_Angle'][0]
ptchChoiceIdx = np.abs(Pitch - 10).argmin()
exampleFitted_data = data_dict_diffNFlux['Differential_Number_Flux'][0][targetidx][ptchChoiceIdx]

# get the model parameters
modelParamsIdx = np.abs(data_dict_dispersiveFitData['timestamp_fitData'][0] - targetTime).argmin()


exampleFit_Density = data_dict_dispersiveFitData['n'][0][wExamplePtch][modelParamsIdx]
exampleFit_Temperature = data_dict_dispersiveFitData['Te'][0][wExamplePtch][modelParamsIdx]
exampleFit_Potential = data_dict_dispersiveFitData['V0'][0][wExamplePtch][modelParamsIdx]
exampleFit_ChiSquare = data_dict_dispersiveFitData['ChiSquare'][0][wExamplePtch][modelParamsIdx]
exampleFit_Energies = data_dict_dispersiveFitData['fitted_Energies'][0][wExamplePtch][modelParamsIdx]
exampleFit_diffNFlux = data_dict_dispersiveFitData['fitted_NFluxes'][0][wExamplePtch][modelParamsIdx]
exampleFit_StdDev = data_dict_diffNFlux['Differential_Number_Flux_stdDev'][0][targetidx][ptchChoiceIdx]

# plot the example fit
axExampleFit.set_title(rf'Example Fit ($\alpha$ = {data_dict_dispersiveFitData["pitchAngle"][0][wExamplePtch][0]}$^\circ$, {data_dict_dist["Epoch"][0][targetidx].strftime("%H:%M:%S.%f")} UTC)',weight='bold',fontsize=Title_FontSize)
axExampleFit.plot(exampleFit_Energies,exampleFit_diffNFlux,color='red',label=rf'n$_0$ =   {round(exampleFit_Density,2)} [cm$^{-3}$]' +'\n'+
                                                                             rf'$T_e=$ {round(exampleFit_Temperature,1)} [eV]' +'\n'+
                                                                             rf'$V_0=$ {round(exampleFit_Potential,0)} [eV]'+'\n'+
                                                                             rf'$\chi^2_\nu=$ {round(exampleFit_ChiSquare,2)}',zorder=2)
diffNFlux_NoiseCount = generateNoiseLevel(data_dict_dist['Energy'][0],countNoiseLevel)
axExampleFit.plot(data_dict_dist['Energy'][0],diffNFlux_NoiseCount,color='black',label=f'{countNoiseLevel}-count noise')
axExampleFit.vlines(exampleFit_Potential,ymin=0,ymax=1E13,color='red',linewidth=Plot_Linewidth,zorder=1)
axExampleFit.errorbar(data_dict_dist['Energy'][0], exampleFitted_data,yerr=exampleFit_StdDev,color='tab:blue',marker='o',capsize=Errorbar_capsize,zorder=0)
axExampleFit.set_yscale('log')
axExampleFit.set_ylabel('[cm$^{-2}$s$^{-1}$sr$^{-1}$eV$^{-1}$]',fontsize=Label_FontSize,weight='bold')
axExampleFit.set_ylim(3E4, 5E7)
axExampleFit.set_xscale('log')
axExampleFit.set_xlim(25, 1000)
axExampleFit.legend(fontsize=Legend_fontSize,loc='upper right')
axExampleFit.set_xlabel('Energy [eV]',fontsize=Label_FontSize,weight='bold')
axExampleFit.tick_params(axis='y', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
axExampleFit.tick_params(axis='x', labelsize=Tick_FontSize-5, length=Tick_Length, width=Tick_Width)
axExampleFit.grid(True, alpha=0.4,which='both')
axExampleFit.text(250, 5E4, s='Primaries',fontsize=Text_FontSize,weight='bold',va='center',ha='center')
axExampleFit.text(70, 5E4, s='Secondaries/Backscatter', fontsize=Text_FontSize, weight='bold', va='center', ha='center')

# ############################
# --- DISTRIBUTIONS PLOTS ---
# ############################
theoreticalMaximumPSengy = exampleFit_Energies[-1] - exampleFit_Energies[0] # Technically what it should be
theoreticalMaximumPSengy = 1+390.74 - exampleFit_Energies[0] # what it is if I'm above noise count 4
velSpace_engy = [-1*engy for engy in data_dict_dist['Energy'][0] if engy <= theoreticalMaximumPSengy]+[0] + [engy for engy in data_dict_dist['Energy'][0] if engy <= theoreticalMaximumPSengy] # in eV
velSpace_vel =  [-1*np.sqrt(2*stl.q0*engy/stl.m_e) for engy in data_dict_dist['Energy'][0] if engy <= theoreticalMaximumPSengy]+[0] + [np.sqrt(2*stl.q0*engy/stl.m_e) for engy in data_dict_dist['Energy'][0] if engy <= theoreticalMaximumPSengy] # in eV

def generate_Maxwellian_Espace(self, n, T, energy_Grid):
    '''
    :param n: density [cm^-3]
    :param T: Temperature [eV]
    :param energy_Grid: Energy Grid [eV]
    :return: distribution function in SI units [s^3 m^-6]
    '''
    return (1E6 * n) * np.power(stl.m_e / (2 * np.pi * stl.q0 * T), 3 / 2) * np.exp(-1 * energy_Grid / T)

dist_PrimaryV = generate_Maxwellian_Espace(avg_fitParams[0][0:2],velSpace_engy)
dist_DisReg = generate_Maxwellian_Espace(avg_fitParams[1][0:2],velSpace_engy)

# plot the distributions
axDist = fig.add_subplot(gs02[0, 1])
axDist.scatter(velSpace_engy,dist_PrimaryV,color='tab:blue',label=f'Dispersive Region Plasma Sheet Model (Avg)')
axDist.scatter(velSpace_engy,dist_DisReg,color='tab:red', label=f'Primary inverted-V Plasma Sheet Model (Avg)')
distDispersive = data_dict_dist['Distribution_Function'][0][targetidx][ptchChoiceIdx]
distDispersive_vel = [np.sqrt(2*stl.q0*engyVal/stl.m_e) for engyVal in data_dict_dist['Energy'][0]]
distDispersive_engy = [engyVal for engyVal in data_dict_dist['Energy'][0]]
axDist.set_title('Model Comparison',fontsize=Title_FontSize,weight='bold')
axDist.scatter(distDispersive_engy, distDispersive,color='black',label=f'Inverted-V (T={data_dict_dist["Epoch"][0][targetidx].strftime("%H:%M:%S.%f")} UTC)')
axDist.set_ylim(1E-17, 1E-12)
axDist.set_yscale('log')
axDist.set_xlim(25, 1E3)
axDist.set_xscale('log')
axDist.set_xlabel('Energy [eV]',fontsize=Label_FontSize,weight='bold')
axDist.set_ylabel('Dist. Function [$\mathbf{m^{-6}s^{3}}$]',fontsize=Label_FontSize-3,weight='bold')
axDist.legend(fontsize=Legend_fontSize,loc='upper right')
axDist.tick_params(axis='x', labelsize=Tick_FontSize-5, length=Tick_Length, width=Tick_Width)
axDist.tick_params(axis='y', labelsize=Tick_FontSize-5, length=Tick_Length, width=Tick_Width)
axDist.yaxis.set_label_position("right")
axDist.yaxis.tick_right()
axDist.grid(True,alpha=0.4,which='both')


# --- --- --- --- --- -
# --- SAVE THE PLOT ---
# --- --- --- --- --- -

fig.subplots_adjust(left=0.07, bottom=0.06, right=0.91, top=0.97, wspace=0.04,hspace=0.04)  # remove the space between plots
plt.savefig(rf'C:\Users\cfelt\Desktop\Research\ACESII\Feltman2025_ACESII_Alfven_Observations\PLOTS\Plot7p5\Plot7p5_base.png', dpi=dpi)
stl.Done(start_time)