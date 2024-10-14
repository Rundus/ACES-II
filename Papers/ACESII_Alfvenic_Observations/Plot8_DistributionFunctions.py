# --- Plot8_DistributionFunctions.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Plot the mapped distribution function at various altitudes


# TODO: Handle the betavalue in the backscatter mapping better because beta=1 should correspond to AT the inverted-V altitude
# TODO: NOT at the ionosphere (as it is currently programmed)

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
print(color.UNDERLINE + f'Plot8_DistributionFunctions' + color.END)

from myspaceToolsLib.time import EpochTo_T0_Rocket
import matplotlib.gridspec as gridspec
from functools import partial
from ACESII_code.Science.InvertedV.Evans_class_var_funcs import *
from scipy.interpolate import CubicSpline

#########################
# --- Physics Toggles ---
#########################

# --- Collect the Data ---
datasetReduction_TargetTime = [dt.datetime(2022, 11, 20, 17, 24, 50, 000000), dt.datetime(2022,11,20,17,25,15,000000)]
targetVar = [datasetReduction_TargetTime, 'Epoch']
mappingAltitudes = [6000, 5072, 3131]


# note: Beta = 5 0.817 --> 5210.826

# --- PRIMARY BEAM diffNFlux fit  ---
primaryBeam_TargetTimes_data = dt.datetime(2022, 11, 20, 17, 25, 1, 212210)
NSweeps = 2 # the WIDTH of How many Sweeps to Average Over. E.g. if NSweeps = 2 --> +2 sweeps on either side of the center point
wPitchToPlot = 3
fitParameters = [2.2, 104.8, 231.8]
diffNFlux_Fit_threshEngy = 100
countNoiseLevel = 4

# --- Model Primary Beam toggles ---
N = 401
PS_BeamThreshEnergy = 621.29*(1)-231.8 # in eV
# PS_BeamThreshEnergy = 621.29 # in eV
PS_BeamThreshPitch = 90
InV_BeamThreshEnergy = 621.29*(1)-231.8 # in eV
# InV_BeamThreshEnergy = 621.29*(1) # in eV
InV_BeamThreshPitch = 90

# backscatter toggles
BackScatter_TargetTimes_AverageData = [dt.datetime(2022, 11, 20, 17, 25, 1, 312207), dt.datetime(2022, 11, 20, 17, 25, 1, 512208)]
BackScatter_EnergyThreshold = [28, fitParameters[2]] # in eV
BackScatter_wPitchAngles = np.array([10, 11, 12, 13, 14, 15, 16, 17, 18, 19])

# Other Plot toggles
Plot_Individually = False
Plot_useVelGrid =False
Plot_useDistributionFunction = False # use distribution Function instead of diffNFlux

# Loss Cone
LossConeAlt = 100 # in km

# --- Alfven Wave Resonance Bands ---
from ACESII_code.Science.AlfvenSingatureAnalysis.Simulations.TestParticle.simToggles import EToggles
waveFreq = EToggles.waveFreq_Hz




################################
# --- Plot toggles - General ---
################################
figure_width = 16 # in inches
figure_height = 22 # in inches
Title_FontSize = 30
Label_FontSize = 25
Label_Padding = 8
Line_LineWidth = 2.5
Text_Fontsize = 25
Tick_FontSize = 25
Tick_FontSize_minor = 20
Tick_Length = 10
Tick_Width = 2
Tick_Length_minor = 5
Tick_Width_minor = 1
Plot_LineWidth = 0.5
Plot_MarkerSize = 14
Legend_fontSize = 20
dpi = 100

# --- Velocity Space ---
VelSpace_Norm = 1E7
VelSpace_Max = 2.8

# --- Cbar ---
mycmap = apl_rainbow_black0_cmap()
if Plot_useDistributionFunction:
    cbarMin, cbarMax = 1E-22, 1E-12
else:
    # cbarMin, cbarMax = 1E3, 5E7
    cbarMin, cbarMax = 5E4, 5E6

cbar_TickLabelSize = 20
cbar_LabelFontSize = 35


##########################
# --- --- --- --- --- ---
# --- LOADING THE DATA ---
# --- --- --- --- --- ---
##########################
prgMsg('Loading Data')

rocketAttrs, b, c = ACES_mission_dicts()

# EEPAA Distribution Data
inputEEPAA_dist_high = glob('C:\Data\ACESII\L3\DistFunc\high\*ACESII_36359_distFunc_eepaa_fullCal*')[0]
data_dict_dist = loadDictFromFile(inputFilePath=inputEEPAA_dist_high, targetVar=targetVar, wKeys_Reduce=['Distribution_Function', 'Epoch'])
inputEEPAA_counts_high = glob('C:\Data\ACESII\L1\high\*eepaa_fullCal*')[0]
data_dict_counts = loadDictFromFile(inputFilePath=inputEEPAA_counts_high, targetVar=targetVar, wKeys_Reduce=['eepaa', 'Epoch'])
inputEEPAA_diffFlux_high = glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')[0]
data_dict_diffFlux = loadDictFromFile(inputEEPAA_diffFlux_high, targetVar=targetVar, wKeys_Reduce=['Differential_Energy_Flux','Differential_Number_Flux', 'Epoch'])

# Define the data
counts = deepcopy(data_dict_counts['eepaa'][0])
distFunc = deepcopy(data_dict_dist['Distribution_Function'][0])
diffEFlux = deepcopy(data_dict_diffFlux['Differential_Energy_Flux'][0])
diffNFlux = deepcopy(data_dict_diffFlux['Differential_Number_Flux'][0])
Epoch = deepcopy(data_dict_dist['Epoch'][0])
Energy = deepcopy(data_dict_dist['Energy'][0])
Pitch = deepcopy(data_dict_dist['Pitch_Angle'][0])
Done(start_time)

###########################################
# --- --- --- --- --- --- --- --- --- --- -
# --- CALCULATED VARIABLES/COLLECT DATA ---
# --- --- --- --- --- --- --- --- --- --- -
###########################################
prgMsg('Preparing Data')

####################################
###### noise THRESHOLD LINE ######
####################################
diffFlux_NoiseCount = np.zeros(shape=(len(Energy)))
geo_factor = rocketAttrs.geometric_factor[0]
count_interval = 0.8992E-3

for engy in range(len(Energy)):
    deltaT = (count_interval - countNoiseLevel * rocketAttrs.deadtime[0])
    diffFlux_NoiseCount[engy] = countNoiseLevel / (geo_factor[0] * deltaT*Energy[engy])

# Time since launch
LaunchDateTime = pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 20, 00, 000000))
Epoch_timeSince = EpochTo_T0_Rocket(InputEpoch=Epoch, T0=LaunchDateTime)

#####################################
# --- TOP PLOT: diffNFlux fitting ---
#####################################
diffNFlux_ChoiceIdx = np.abs(Epoch - primaryBeam_TargetTimes_data).argmin()
EngyIdx = np.abs(Energy - diffNFlux_Fit_threshEngy).argmin()
EpochSliceTimeSinceLaunchValue = Epoch_timeSince[diffNFlux_ChoiceIdx]

# --- Perform the fit ---
fitData = diffNFlux[diffNFlux_ChoiceIdx][wPitchToPlot]
peakDiffNVal = fitData[:EngyIdx].max()
peakDiffNVal_index = np.argmax(fitData[:EngyIdx])
xData_fit = np.array(Energy[:peakDiffNVal_index+1])
yData_fit = np.array(fitData[:peakDiffNVal_index+1])
nonZeroIndicies = np.where(yData_fit!=0)[0]
xData_fit = xData_fit[nonZeroIndicies]
yData_fit = yData_fit[nonZeroIndicies]
fittedX = np.linspace(xData_fit.min(), xData_fit.max(), 100)
fitFuncAtPitch = partial(diffNFlux_for_mappedMaxwellian, alpha=Pitch[wPitchToPlot])
fittedY = fitFuncAtPitch(fittedX, *fitParameters)

# --- STD dev ---
dataToAverage = np.transpose(diffNFlux[(diffNFlux_ChoiceIdx-NSweeps):(diffNFlux_ChoiceIdx+NSweeps), wPitchToPlot, :]) # Collect the data. Note: Need ALL the energies if I want the error bars on all datapoints
data_stdDevs = np.array([np.std(engyArr[np.where(engyArr > 0)]) for engyArr in dataToAverage]) / np.sqrt(2*NSweeps+1)

# --- ChiSquare ---
stdDev_myEnergies = data_stdDevs[:peakDiffNVal_index+1]
stdDev_fit = stdDev_myEnergies[np.where(stdDev_myEnergies>=0)]# get the error bars only for the data in the fit range
ChiSquare = []
nu = 1/(3-1)
for i in range(len(xData_fit)):
    yVal = yData_fit[i]
    fxVal = fitFuncAtPitch(xData_fit[i],*fitParameters)
    error = stdDev_fit[i]
    if error == 0:
        error = yVal
    ChiSquare.append(((fxVal - yVal)/error)**2)

ChiSquare = nu*sum(ChiSquare)


Done(start_time)

##################################################
# --- VELOCITY SPACE GRID VALUES USED IN MODEL ---
##################################################
Emax_vel = np.sqrt(2*q0*2000/m_e)
Vperp_Vals = np.linspace(-1*Emax_vel, Emax_vel, N)
Vpara_Vals = np.linspace(0, Emax_vel, N)
Vperp_Grid, Vpara_Grid = np.meshgrid(Vperp_Vals, Vpara_Vals)
DistFunc_Grid = dist_Maxwellian(Vperp=Vperp_Grid,Vpara=Vpara_Grid,model_Params=fitParameters)

###########################
# --- BackScatter Data ---
###########################
prgMsg('Calculating BackScatter')

# Determine the Average distribution Function across all energies up to BackScatter_EnergyThreshold
# thresholdIdx = np.abs(Energy - Energy[peakDiffNVal_index]).argmin()
BackScatter_Energies = Energy[peakDiffNVal_index:]

# handle fillvals and Average the lower energies over all pitch angles
tempData = np.array(diffNFlux[diffNFlux_ChoiceIdx, BackScatter_wPitchAngles, peakDiffNVal_index:]).T
backScatterData = []
for engyRow in tempData:
    temp = [val for val in engyRow if np.abs(val) <= 1E11]
    backScatterData.append(sum(temp)/len(temp))

# Sort the data
BackScatter_Energies, backScatterData = zip(*sorted(zip(BackScatter_Energies,backScatterData)))


# interpolate the curve onto the Primary Beam's Energy range
BackScatter_diffNFlux_spline = CubicSpline(np.array(BackScatter_Energies), np.array(backScatterData))
BackScatter_diffNFlux_Grid = calc_BackScatter_onto_Velspace(VperpGrid=Vperp_Grid,
                                                            VparaGrid=Vpara_Grid,
                                                            BackScatterSpline=BackScatter_diffNFlux_spline,
                                                            EngyLimit=BackScatter_EnergyThreshold)
# print(BackScatter_diffNFlux_Grid)
BackScatter_distFunc_Grid = diffNFlux_to_distFunc(Vperp=Vperp_Grid,
                                                  Vpara=Vpara_Grid,
                                                  diffNFlux=BackScatter_diffNFlux_Grid)

# Map the backscatter FROM the ionosphere to the Magnetosphere at various beta values
Data_BackScatter_diffNFluxAtBeta = []
Data_BackScatter_VperpAtBeta = []
Data_BackScatter_VparaAtBeta = []
for mapAlt in mappingAltitudes:
    VperpGrid_mapped, VparaGrid_mapped, diffNFlux_mapped = mapping_VelSpace_magMirror(VperpGrid=Vperp_Grid,
                                                                                       VparaGrid=Vpara_Grid,
                                                                                       distFuncGrid=BackScatter_distFunc_Grid,
                                                                                       targetAlt=400,
                                                                                       startingAlt=mapAlt,
                                                                                       mapToMagSph=True)
    Data_BackScatter_VperpAtBeta.append(VperpGrid_mapped)
    Data_BackScatter_VparaAtBeta.append(VparaGrid_mapped)
    Data_BackScatter_diffNFluxAtBeta.append(diffNFlux_mapped)


Done(start_time)


##########################################
# --- Primary Beam + BackScatter Model ---
##########################################
prgMsg('Calculating Primary Beam')

# --- INVERTED V AT DIFFERENT BETA ---
Data_InvertedV_VperpAtBeta = []
Data_InvertedV_VparaAtBeta = []
Data_InvertedV_diffNFluxAtBeta = []

# Generate the new ACCELERATED distribution
VperpGrid_inV, VparaGrid_inV, distGrid_inV, diffNFluxGrid_inV = calc_velSpace_DistFuncDiffNFluxGrid(Vperp_gridVals=Vperp_Vals,
                                                                                    Vpara_gridVals=Vpara_Vals,
                                                                                    model_Params=fitParameters,
                                                                                    initalBeamParams=[InV_BeamThreshPitch,InV_BeamThreshEnergy])

for mapAlt in mappingAltitudes:

    # Map it to a different beta altitude

    VperpGrid_mapped, VparaGrid_mapped, diffNFlux_mapped = mapping_VelSpace_magMirror(
        VperpGrid=VperpGrid_inV,
        VparaGrid=VparaGrid_inV,
        distFuncGrid=distGrid_inV,
        targetAlt=mapAlt,
        startingAlt=max(mappingAltitudes),
        mapToMagSph=False)
    Data_InvertedV_VperpAtBeta.append(VperpGrid_mapped)
    Data_InvertedV_VparaAtBeta.append(VparaGrid_mapped)
    Data_InvertedV_diffNFluxAtBeta.append(diffNFlux_mapped)

# --- PLASMA SHEET AT DIFFERENT BETAs ---
Data_PS_VperpAtBeta = []
Data_PS_VparaAtBeta = []
Data_PS_diffNFluxAtBeta = []
fitParamsNoAccel = deepcopy(fitParameters)
fitParamsNoAccel[-1] = 0
# Generate the new ACCELERATED distribution
VperpGrid_PS, VparaGrid_PS, distGrid_PS, diffNFluxGrid_PS = calc_velSpace_DistFuncDiffNFluxGrid(
    Vperp_gridVals=Vperp_Vals,
    Vpara_gridVals=Vpara_Vals,
    model_Params=fitParamsNoAccel,
    initalBeamParams=[PS_BeamThreshPitch, PS_BeamThreshEnergy])

for mapAlt in mappingAltitudes:

    # Map it to a different beta altitude
    VperpGrid_mapped, VparaGrid_mapped, diffNFlux_mapped = mapping_VelSpace_magMirror(
        VperpGrid=VperpGrid_PS,
        VparaGrid=VparaGrid_PS,
        distFuncGrid=distGrid_PS,
        targetAlt=mapAlt,
        startingAlt=max(mappingAltitudes),
        mapToMagSph=False)

    for rowIdx in range(len(VperpGrid_mapped)):
        for colIdx in range(len(VperpGrid_mapped[0])):

            specificPitchAngle = np.arctan2(VperpGrid_mapped[rowIdx][colIdx],VparaGrid_mapped[rowIdx][colIdx])
            if VperpGrid_mapped[rowIdx][colIdx] < 0:
                specificPitchAngle = -1*specificPitchAngle

            lostPitchThresh = np.arcsin(np.sin(specificPitchAngle)*np.power((Re + mapAlt)/(Re + LossConeAlt),3/2))
            if -1*lostPitchThresh <= specificPitchAngle <= lostPitchThresh:
                diffNFlux_mapped[rowIdx][colIdx] = 0


    Data_PS_VperpAtBeta.append(VperpGrid_mapped)
    Data_PS_VparaAtBeta.append(VparaGrid_mapped)
    Data_PS_diffNFluxAtBeta.append(diffNFlux_mapped)

Done(start_time)


##################################################
# --- COLLECT INERTIAL ALFVEN SPEED MODEL DATA ---
##################################################
# Use my own model data and my chosen potential to determine some resonance Velocities
data_dict_alfModel = loadDictFromFile(r'C:\Data\ACESII\science\simulations\TestParticle\plasmaEnvironment\plasmaEnvironment.cdf')
targetAltAlfvenSim = 10000*1000 # in kilometers
targetIdx_Sim = np.abs(data_dict_alfModel['simAlt'][0] - targetAltAlfvenSim).argmin()
altRange_sim = data_dict_alfModel['simAlt'][0][:targetIdx_Sim]
alfSpeedMHD_sim = data_dict_alfModel['alfSpdMHD'][0][:targetIdx_Sim]
alfSpeedKinetic_sim = data_dict_alfModel['alfSpdInertial'][0][:targetIdx_Sim]

ResonanceLimits = []
AlfSpeedAtSpecificAlti = []
PhiPotential = []
for d, mapAlt in enumerate(mappingAltitudes):
    Altitude = mapAlt*1000
    altIdx = np.abs(data_dict_alfModel['simAlt'][0] - Altitude).argmin()
    localAlfSpdInertial = data_dict_alfModel['alfSpdInertial'][0][altIdx]
    # PhiVal = 2*EparallelMax[d]*localAlfSpdInertial/waveFreq
    AlfSpeedAtSpecificAlti.append(localAlfSpdInertial)
    if d == 2:
        PhiVal = round(0.5 * m_e * np.power(localAlfSpdInertial - np.sqrt(2 * (28) * q0 / m_e), 2) / q0)
        print(Altitude, localAlfSpdInertial, 0.5 * m_e * np.power(localAlfSpdInertial - np.sqrt(2 * (28) * q0 / m_e), 2) / q0)
    else:
        PhiVal = round(0.5 * m_e * np.power(localAlfSpdInertial - np.sqrt(2 * (InV_BeamThreshEnergy + fitParameters[-1]) * q0 / m_e), 2) / q0)
        print(Altitude, localAlfSpdInertial, 0.5 * m_e * np.power(localAlfSpdInertial - np.sqrt(2 * (InV_BeamThreshEnergy + fitParameters[-1]) * q0 / m_e), 2) / q0)

    PhiPotential.append(PhiVal)
    ResonanceLimits.append([localAlfSpdInertial - np.sqrt(2 * q0 * PhiVal / (m_e)), localAlfSpdInertial])





















################################
# --- --- --- --- --- --- --- --
# --- PLOT THE DISTRIBUTIONS ---
# --- --- --- --- --- --- --- --
################################
prgMsg('Plotting Fig')

# --- plot parameters ---

if Plot_Individually:
    fig, ax_diffNFluxFit = plt.subplots()
    fig.set_size_inches(13, 6)
else:
    fig = plt.figure()
    fig.set_size_inches(figure_width, figure_height)
    gs0 = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[3 / 16, 2 / 16, 11 / 16])
    ax_diffNFluxFit = fig.add_subplot(gs0[0, :])

# --- --- --- --- --- --- --
# --- diffNFlux fit plot ---
# --- --- --- --- --- --- --

ax_diffNFluxFit.errorbar(x=Energy, y=fitData, yerr=data_stdDevs,xerr=None, color='tab:blue', linewidth=Line_LineWidth, label= f'{Pitch[wPitchToPlot]}'+'$^{\circ}$ Pitch Data'+f' (T = {round(EpochSliceTimeSinceLaunchValue,1)} sec)', capsize=7,zorder=1)
ax_diffNFluxFit.plot(Energy, fitData, color='tab:blue', marker='o', markersize=4)
ax_diffNFluxFit.plot(BackScatter_Energies,backScatterData,marker='o',color='tab:purple',label='Avg. BackScatter (Up Going)', linewidth=Line_LineWidth)
ax_diffNFluxFit.plot(fittedX, fittedY, color='tab:red',label='Maxwellian Fit', linewidth=Line_LineWidth,zorder=2)
ax_diffNFluxFit.text(x=1000,y=5E5,s=f'      $n_0$={round(fitParameters[0],1)}' +' [cm$^{-3}$]   '+ f'T={round(fitParameters[1],1)} [eV]\n'
                                                          + f'$V_0$={round(fitParameters[2],1)} [eV]   '+r'$\chi^{2}_{\nu}$='+f'{round(ChiSquare,3)}', fontsize=Text_Fontsize-5, color='tab:red', va='center',ha='center')
ax_diffNFluxFit.text(x=100, y=7.7E4, s='Secondary/BackScatter Flux', fontsize=Text_Fontsize, ha='center')
ax_diffNFluxFit.text(x=1180, y=7.7E4, s='Primary Beam Flux', fontsize=Text_Fontsize,ha='center')
ax_diffNFluxFit.axvline(Energy[peakDiffNVal_index], color='black',linestyle='--')
ax_diffNFluxFit.set_yscale('log')
ax_diffNFluxFit.set_xscale('log')
ax_diffNFluxFit.set_xlabel('Energy [eV]', fontsize=Label_FontSize,labelpad=-15)
ax_diffNFluxFit.set_ylabel('[cm$^{-2}$s$^{-1}$str$^{-1}$ eV$^{-1}$]', fontsize=Label_FontSize)
ax_diffNFluxFit.set_xlim(28, 2E3)
ax_diffNFluxFit.set_ylim(6E4, 5E7)
ax_diffNFluxFit.plot(Energy, diffFlux_NoiseCount, color='black', label=f'{countNoiseLevel}-count level')
ax_diffNFluxFit.legend(fontsize=Legend_fontSize, loc='upper right')
ax_diffNFluxFit.tick_params(axis='y', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
ax_diffNFluxFit.tick_params(axis='y', which='minor', labelsize=Tick_FontSize_minor, width=Tick_Width_minor, length=Tick_Length_minor)
ax_diffNFluxFit.tick_params(axis='x', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
ax_diffNFluxFit.tick_params(axis='x', which='minor', labelsize=Tick_FontSize_minor, width=Tick_Width_minor, length=Tick_Length_minor)



# --- --- --- --- --- --- --
# --- diffNFlux fit plot ---
# --- --- --- --- --- --- --
if Plot_Individually:
    plt.tight_layout()
    fig.savefig(r'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\diffNFitting.png', dpi=dpi)
    fig, ax_alfvenSpeed = plt.subplots()
    fig.set_size_inches(12, 6)
else:
    ax_alfvenSpeed = fig.add_subplot(gs0[1, :])
ax_alfvenSpeed.plot(deepcopy(altRange_sim)/1000, deepcopy(alfSpeedMHD_sim)/VelSpace_Norm, color='blue',linewidth=Line_LineWidth+1, label='MHD',zorder=1)
ax_alfvenSpeed.plot(deepcopy(altRange_sim)/1000, deepcopy(alfSpeedKinetic_sim)/VelSpace_Norm, label='Inertial', color='red', linewidth=Line_LineWidth+1,zorder=1)
ax_alfvenSpeed.text(x=1250,y=3.5,s=r'$\lambda_{\perp 0}$ = '+f'{EToggles.lambdaPerp0/1000} km',fontsize=Text_Fontsize, ha='center', va='center',bbox=dict(facecolor='white', edgecolor='black',boxstyle='round'))
ax_alfvenSpeed.set_ylabel(r'Alfv$\'e$n Speed'+'\n  [10000 km/s]',fontsize=Label_FontSize)
ax_alfvenSpeed.set_xlabel(f'Altitude [km]', fontsize=Label_FontSize)
ax_alfvenSpeed.grid(True)
ax_alfvenSpeed.set_ylim(0,4.2)
ax_alfvenSpeed.set_xlim(100, 10000)
ax_alfvenSpeed.legend(loc='lower right',fontsize=Legend_fontSize)
keyPoints_alt = np.array(mappingAltitudes)
keyPoints_VA = np.array(AlfSpeedAtSpecificAlti)/VelSpace_Norm
ax_alfvenSpeed.scatter(x=keyPoints_alt,y=keyPoints_VA,color='red',s=200,zorder=2)
ax_alfvenSpeed.tick_params(axis='y', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
ax_alfvenSpeed.tick_params(axis='y', which='minor', labelsize=Tick_FontSize_minor, width=Tick_Width_minor, length=Tick_Length_minor)
ax_alfvenSpeed.tick_params(axis='x', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
ax_alfvenSpeed.tick_params(axis='x', which='minor', labelsize=Tick_FontSize_minor, width=Tick_Width_minor, length=Tick_Length_minor)

if Plot_Individually:
    plt.tight_layout()
    fig.savefig(r'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\AflvenSpeed.png', dpi=dpi)


# --- --- --- --- --- --- --- --- --- --- -
# --- distribution func plot - HIGH ALT ---
# --- --- --- --- --- --- --- --- --- --- -
if not Plot_Individually:
    gs00 = gridspec.GridSpecFromSubplotSpec(len(mappingAltitudes), 2, subplot_spec=gs0[2, :])

    VperpData = [Data_PS_VperpAtBeta, Data_InvertedV_VperpAtBeta]
    VparaData = [Data_PS_VparaAtBeta, Data_InvertedV_VparaAtBeta]
    ZData = [Data_PS_diffNFluxAtBeta, Data_InvertedV_diffNFluxAtBeta] if not Plot_useDistributionFunction else [[DistFunc_Grid for i in range(len(mappingAltitudes))], [DistFunc_Grid for i in range(len(mappingAltitudes))]]

    for i in range(len(mappingAltitudes)): # Which Beta I'm considering
        for j in range(2): # Plasma Sheet vs InvertedV

            ax = fig.add_subplot(gs00[i, j])

            ###############################################
            # --- VELOCITY SPACE GRID INSTRUMENTAL DATA ---
            ###############################################
            xData = VperpData[j][i]
            yData = VparaData[j][i]
            xData_Back = Data_BackScatter_VperpAtBeta[i]
            yData_Back = Data_BackScatter_VparaAtBeta[i]

            if Plot_useVelGrid:
                # --- Primary Beam ---
                PitchInterp = [-350 + i * 10 for i in range(72)]
                EnergyBins = np.logspace(0, 4, base=10, num=42)
                ZGrid_New_Prim, EnergyGrid_Instr_Prim, PitchGrid_Instr_Prim = velocitySpace_to_PitchEnergySpace(EnergyBins=EnergyBins, PitchBins=PitchInterp, VperpGrid=xData, VparaGrid=yData, ZGrid=ZData[j][i],method='average')

                if j == 1:
                    # --- Secondaries ---
                    ZGrid_New_BackScatter, EnergyGrid_Instr_BackScatter, PitchGrid_Instr_BackScatter = velocitySpace_to_PitchEnergySpace(EnergyBins=EnergyBins, PitchBins=PitchInterp, VperpGrid=xData_Back, VparaGrid=yData_Back, ZGrid=Data_BackScatter_diffNFluxAtBeta[i],method='average')
                    newZData = ZGrid_New_Prim + ZGrid_New_BackScatter
                else:
                    newZData = ZGrid_New_Prim

                # --- Create the Velocity Space Grid ---
                Vperp = deepcopy(ZGrid_New_Prim)
                Vpara = deepcopy(ZGrid_New_Prim)

                for ptch in range(len(ZGrid_New_Prim)):
                    for engy in range(len(EnergyBins)):
                        Vmag = np.sqrt(2 * q0 * EnergyBins[engy] / m_e)
                        Vperp[ptch][engy] = np.sin(np.radians(PitchInterp[ptch])) * Vmag
                        Vpara[ptch][engy] = np.cos(np.radians(PitchInterp[ptch])) * Vmag

                Vpara, Vperp = np.array(Vpara) / (VelSpace_Norm), np.array(Vperp) / (VelSpace_Norm)
                cmapObj = ax.pcolormesh(Vperp, Vpara, newZData, vmin=cbarMin, vmax=cbarMax, cmap=mycmap, norm='log')
            else:
                cmapObj = ax.pcolormesh(xData/VelSpace_Norm, yData/VelSpace_Norm, ZData[j][i], vmin=cbarMin, vmax=cbarMax, cmap=mycmap, norm='log')

            #########################
            # --- RESONANCE BANDS ---
            #########################
            ax.axhline(ResonanceLimits[i][0] / VelSpace_Norm, linestyle='--', color='black', linewidth=Line_LineWidth)
            # ax.text(x=0.65, y= sum(ResonanceLimits[i])/(2* VelSpace_Norm)+0.075, s=r'$|\frac{\omega}{k} - v_{\parallel} |< \sqrt{\frac{2e\phi}{m_{e}}}$', fontsize=Text_Fontsize-2)

            if i == 0:
                ax.text(x=0.7, y=sum(ResonanceLimits[i]) / (2 * VelSpace_Norm) + 0.075, s=r'$\Phi_{P} = $' + f'{PhiPotential[i]} eV', fontsize=Text_Fontsize - 2, ha='left')
            if i == 1:
                ax.text(x=0.7, y=sum(ResonanceLimits[i]) / (2 * VelSpace_Norm) + 0.075, s=r'$\Phi_{max} = $' + f'{PhiPotential[i]} eV', fontsize=Text_Fontsize - 2, ha='left')
            if i == 2:
                ax.text(x=0.7, y=sum(ResonanceLimits[i]) / (2 * VelSpace_Norm) + 0.075, s=r'$\Phi_{min} = $' + f'{PhiPotential[i]} eV', fontsize=Text_Fontsize - 2, ha='left')
            ax.axhline(ResonanceLimits[i][1] / VelSpace_Norm, color='black', linewidth=Line_LineWidth)
            ax.fill_between([-1 * VelSpace_Max, VelSpace_Max], ResonanceLimits[i][0] / VelSpace_Norm, ResonanceLimits[i][1] / VelSpace_Norm, color='grey', alpha=0.35)

            # # General Labels
            ax.set_ylim(0, VelSpace_Max)
            ax.set_xlim(-VelSpace_Max, VelSpace_Max)
            ax.invert_yaxis()
            ax.tick_params(axis='y', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
            ax.tick_params(axis='y', which='minor', labelsize=Tick_FontSize_minor, width=Tick_Width_minor, length=Tick_Length_minor)
            ax.tick_params(axis='x', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
            ax.tick_params(axis='x', which='minor', labelsize=Tick_FontSize_minor, width=Tick_Width_minor, length=Tick_Length_minor)
            ax.text(x=-VelSpace_Max+0.2, y=0.15, s=rf'{mappingAltitudes[i]} km', fontsize=Label_FontSize,bbox=dict(facecolor='white', edgecolor='black',boxstyle='round'),va='top',ha='left')
            # points = ['A','B','C']
            # ax.text(x=-VelSpace_Max+0.2, y=0.15, s=rf'{points[i]}', fontsize=Label_FontSize+5,bbox=dict(facecolor='white', edgecolor='black',boxstyle='round'),va='top',ha='left')

            # Specific things
            if i == 0 and j == 0:
                ax.set_title('Ambient Plasma Sheet', fontsize=Title_FontSize)
            if i == 0 and j == 1:
                ax.set_title('Inverted-V', fontsize=Title_FontSize)

            if j==1 and not Plot_useVelGrid: # plot the backscatter
                ax.pcolormesh(Data_BackScatter_VperpAtBeta[i] / VelSpace_Norm, Data_BackScatter_VparaAtBeta[i] / VelSpace_Norm, Data_BackScatter_diffNFluxAtBeta[i], vmin=cbarMin, vmax=cbarMax, cmap=mycmap, norm='log')
            if j == 1 and i == 0: # plot the accelerated potential line
                ax.axhline(np.sqrt(2 * fitParameters[2] * q0 / m_e) / VelSpace_Norm, color='tab:red',linewidth=Line_LineWidth)
                ax.text(x=0.5, y=(np.sqrt(2 * fitParameters[2] * q0 / m_e) / VelSpace_Norm)-0.075,s='$V_{0}$' + f'= {round(fitParameters[2],1)} eV',fontsize=Text_Fontsize, color='black'  )
            # if j == 1 and i == 2:
                # ax.axhline(np.sqrt(2 * 28 * q0 / m_e) / VelSpace_Norm, color='tab:red', linewidth=Line_LineWidth)
                # ax.text(x=1.72, y=(np.sqrt(2 * 28 * q0 / m_e) / VelSpace_Norm) - 0.075, s='28 eV', fontsize=Text_Fontsize, color='black')

            if j == 0: # plot Vparallel labels
                ax.set_ylabel('$V_{\parallel}$ ' + rf'[$10^4$ km/s]', fontsize=Label_FontSize)
            if i == 2:# plot Vperp labels
                ax.set_xlabel('$V_{\perp}$ ' + f'[$10^4$ km/s]', fontsize=Label_FontSize)



    # Cbar
    cax = fig.add_axes([0.90, 0.05, 0.025, 0.564])
    cbar = plt.colorbar(cmapObj,cax=cax)
    cbar.ax.minorticks_on()
    cbar.ax.tick_params(labelsize=cbar_TickLabelSize)
    cbar.ax.get_yaxis().labelpad = 40
    cbar.set_label(r'[cm$^{-2}$str$^{-1}$s$^{-1}$eV$^{-1}$]', fontsize=cbar_LabelFontSize, rotation=270)
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_weight("bold")
        l.set_fontsize(cbar_TickLabelSize)


    plt.subplots_adjust(left=0.08, bottom=0.05, right=0.89, top=0.98, wspace=None, hspace=None)
    try:
        os.remove(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\Plot8_DistributionFunc_Base.png')
    except:
        print('')
    plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot8\Plot8_DistributionFunc_Base.png',dpi=dpi)
    Done(start_time)