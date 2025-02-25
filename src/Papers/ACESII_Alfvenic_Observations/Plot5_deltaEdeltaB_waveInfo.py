# --- Plot3_detlaEdetlaB_waveInfo.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Stack Plot detailing field-aligned
# particle data along with electric and magnetic signatures. Also calculates the Local Alfven Speed
# Where the rocket is


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np
from src.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from numpy.fft import rfft, fftfreq
from scipy.integrate import simpson
from scipy.signal import spectrogram
import spaceToolsLib as stl

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
# figure_width = 7.5 # in inches
figure_width = 9.5 # in inches
figure_height = 4.3*(12/5) # in inches
dpi = 500
PlotLabelSize = 11
PlotTitleSize = 12
PlotLegendSize = 10
PlotTickLabelSize = 10
PlotLabelPad = 1
freqLimit = 15
plotLineWidth = 1
TitlePadding = 10
PoyntingScale = 1E3# convert from W/m^2 to ergs/cm^2
Escale = 1000 # convert V/m to mV/m
EB_ratio_limits = [1E5, 1E7]
plotColors = ['tab:blue', 'tab:red', 'tab:orange', 'tab:green']
RegionNames = ['Quiet Region', 'Temporal Align', 'ILat Align']
kSets = [0, 1] # determines the E/B pairs: (1) B_e/E_r (2) B_r/-E_e or [0,1] for both
NORMALIZE = False

# Phase Histogram
windowType, npersegN, scalingType = 'hann', 64*4, 'spectrum'  # spectrogram toggles
overlap = int(npersegN * (7 / 8))  # hanning filter overlap


# --- TARGET REDUCTION VARIABLES ---
targetTimes = [[dt.datetime(2022, 11, 20, 17, 24, 31, 500), dt.datetime(2022, 11, 20, 17, 24, 38, 00)], # Quiet Time
               [dt.datetime(2022, 11, 20, 17, 24, 54, 250000), dt.datetime(2022, 11, 20, 17, 25, 3, 250000)], # Dispersed Time
               [dt.datetime(2022, 11, 20, 17, 25, 18, 000000), dt.datetime(2022, 11, 20, 17, 25, 27, 000000)]  # Auroral Time ( i.e. the FULL aurora)
               ]

# targetILats = [[71.15, 71.25], # Quiet Time
#                [71.615, 71.67], # temporally Time
#                [71.86, 71.912]] # spatially alignedZ

# trying to make it fit
targetILats = [[71.17, 71.25], # Quiet Time
               [71.615, 71.67], # temporally Time
               [71.865, 71.905]] # spatially alignedZ

# OTHER regions of alfven activity
# extraILats = [[72.1, 72.136],
#               [72.5458, 72.6787],
#               [72.75, 72.8]]

# targetILats = [[72.1, 72.136],
#               [72.61, 72.68],
#               [72.75, 72.8]]

targetVar = targetILats
targetVarName = 'ILat'








##########################
# --- --- --- --- --- ---
# --- PLOTTING START ---
# --- --- --- --- --- ---
#########################


# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---

print(stl.color.UNDERLINE + f'Plot3_WaveAnalysis' + stl.color.END)
stl.prgMsg('Loading Data')

# inputFile_deltaB = glob('C:\Data\ACESII\L3\deltaB\low\*Field_Aligned*')[0] # get the deltaB data
# inputFile_deltaE = glob('C:\Data\ACESII\L3\deltaE\low\*Field_Aligned*')[0] # get the deltaE data
inputFile_deltaB = glob('C:\Data\ACESII\L2\low\ACESII_36364_l2_RingCore_Field_Aligned_HighPass_low0.7_high0.7.cdf')[0] # get the deltaB data
inputFile_deltaE = glob('C:\Data\ACESII\L2\low\ACESII_36364_l2_E_Field_Field_Aligned_HighPass_low0.7_high0.7.cdf')[0] # get the deltaE data

inputFile_B = 'C:\Data\ACESII\L2\low\ACESII_36364_l2_RingCore_Field_Aligned.cdf'  # get the B data
inputFile_E = 'C:\Data\ACESII\L2\low\ACESII_36364_l2_E_Field_Field_Aligned.cdf'  # get the E data
# inputFile_poynting = glob('C:\Data\ACESII\L3\deltaS\low\*Field_Aligned*')[0] # get the Poynting Flux data
inputFile_poynting = glob('C:\Data\ACESII\science\PoyntingFlux\low\*Field_Aligned*')[0] # get the Poynting Flux data

inputFile_ASpeed = glob('C:\Data\ACESII\science\AlfvenSpeed_rkt\low\ACESII_36364_AlfvenSpeed_flight.cdf')[0] # get the Alfven Speed data
inputFile_Langmuir = 'C:\Data\ACESII\L3\Langmuir\low\ACESII_36364_langmuir_fixed_LowPass_low0.3_high0.3.cdf'
inputFile_Bmag = 'C:\Data\ACESII\L1\low\ACESII_36364_l1_RingCore_rktFrm.cdf'
inputFile_EISCAT = 'C:\Data\ACESII\science\EISCAT_ACESII_Slice\low\ACESII_36364_EISCAT_Tromso_rktSlice.cdf'


# --- Break up DataDicts into targetTime sections ---

# FORMAT: [[deltaB, deltaE, poynting, AlfvenSpeedRatio], ...]
dict_sets = []
sectionTimeRange = []

for i in range(len(targetVar)):
    data_dict_E = stl.loadDictFromFile(inputFile_E,targetVar=[targetVar[i], targetVarName])
    data_dict_B = stl.loadDictFromFile(inputFile_B,targetVar=[targetVar[i], targetVarName])
    data_dict_deltaB = deepcopy(stl.loadDictFromFile(inputFile_deltaB,targetVar=[targetVar[i], targetVarName]))
    data_dict_deltaE = deepcopy(stl.loadDictFromFile(inputFile_deltaE,targetVar=[targetVar[i], targetVarName]))
    data_dict_poynting = deepcopy(stl.loadDictFromFile(inputFile_poynting,targetVar=[targetVar[i], targetVarName]))
    data_dict_langmuir = deepcopy(stl.loadDictFromFile(inputFile_Langmuir,targetVar=[targetVar[i], targetVarName],wKeys_Load=['ni', 'Epoch', 'ILat']))
    data_dict_EISCAT = deepcopy(stl.loadDictFromFile(inputFile_EISCAT, targetVar=[targetVar[i], targetVarName],wKeys_Load=['Ion_Comp', 'Op_Comp','ILat','Epoch']))

    # downsample the langmuir data
    indexVals = [np.abs(data_dict_langmuir['ILat'][0] - ILat).argmin() for ilt, ILat in enumerate(data_dict_B['ILat'][0])]
    data_dict_langmuir['ni'][0] = deepcopy(data_dict_langmuir['ni'][0][indexVals])
    data_dict_langmuir['Epoch'][0] = deepcopy(data_dict_langmuir['Epoch'][0][indexVals])
    data_dict_langmuir['ILat'][0] = deepcopy(data_dict_langmuir['ILat'][0][indexVals])

    # Up-sample the EISCAT data (it's fine since the EISCAT variables are very slowly varying)
    data_dict_EISCAT_interp = stl.InterpolateDataDict(InputDataDict=data_dict_EISCAT,InputEpochArray=data_dict_EISCAT['Epoch'][0],targetEpochArray=data_dict_deltaE['Epoch'][0],wKeys=[])
    data_dict_EISCAT = deepcopy(data_dict_EISCAT_interp)

    sectionTimeRange.append([data_dict_deltaB['Epoch'][0][0], data_dict_deltaB['Epoch'][0][-1]])

    # adjust E-field scale
    data_dict_deltaE['E_e'][0] = Escale*data_dict_deltaE['E_e'][0]
    data_dict_deltaE['E_p'][0] = Escale*data_dict_deltaE['E_p'][0]
    data_dict_deltaE['E_r'][0] = Escale*data_dict_deltaE['E_r'][0]

    dict_sets.append([data_dict_deltaB,
                      data_dict_deltaE,
                      data_dict_poynting,
                      data_dict_langmuir,
                      data_dict_B,
                      data_dict_E,
                      data_dict_EISCAT
                      ])


stl.Done(start_time)

#######################################
# --- PREPARE THE DATA FOR PLOTTING ---
#######################################
def Plot3_deltaEdeltaB_waveInfo(targetVar,dict_sets):

    # --- Calculate the FFT of deltaB, deltaE and their ratio ---
    # FORMAT: [[FFT_deltaBe, FFT_deltaBr,xf_deltaB, FFT_deltaEe, FFT_deltaEr,xf_DeltaE, FFT_ASpeed_ErBe, FFT_ASpeed_EeBr], ...]
    # Indicies:[          0,           1,        2,           3,           4,        5,               6,               7]
    data_dict_FFT = {
                'B0' : [],
                'Be_fft':[],
                'Br_fft':[],
                'Ee_fft':[],
                'Er_fft': [],
                'xf_deltaB':[],
                'xf_deltaE': [],
                'Be_norm_fft' : [],
                'Br_norm_fft' : [],
                'Ee_norm_fft' : [],
                'Er_norm_fft' : [],
                'ni' : [],
                'VA_avg' : [],
                'VA_t' : [],
                'xf_B':[],
                'xf_E':[]}

    for i in range(len(targetVar)):
        print(i)

        data_dicts = dict_sets[i]

        # Calculate the normalized Quantities
        B_e = 1E-9 * data_dicts[4]['B_e'][0] # in tesla
        B_p = 1E-9 * data_dicts[4]['B_p'][0]
        B_r = 1E-9 * data_dicts[4]['B_r'][0]
        dB_e = data_dicts[0]['B_e'][0] # in nanotesla
        dB_p = data_dicts[0]['B_p'][0]
        dB_r = data_dicts[0]['B_r'][0]

        B0 = 1E-9 * data_dicts[4]['Bmag'][0]
        # B0 = np.array([np.linalg.norm([dB_e[i], dB_p[i], dB_r[i]]) for i in range(len(data_dicts[0]['Epoch'][0]))])

        E_e = data_dicts[5]['E_e'][0] # in V/m
        E_p = data_dicts[5]['E_p'][0]
        E_r = data_dicts[5]['E_r'][0]

        dE_e = data_dicts[1]['E_e'][0] # in mV/m
        dE_p = data_dicts[1]['E_p'][0]
        dE_r = data_dicts[1]['E_r'][0]

        ni = (stl.cm_to_m ** 3) * data_dicts[3]['ni'][0]

        # calculate the Alfven speed using EISCAT mass concentrations
        # Note: At 180km the primary contributions to Ions are NO+, O2+, N2+, O+
        # With the exception of O+, everything has an atomic mass around 28-32 amu. EISCAT gives these ratios of concentration
        rhoCalc = (ni*data_dicts[6]['Ion_Comp'][0]*((stl.IonMasses[4] + stl.IonMasses[5] + stl.IonMasses[7])/3) + stl.IonMasses[1]*ni*data_dicts[6]['Op_Comp'][0])
        VA_t = B0 / np.sqrt(stl.u0 * rhoCalc)
        VA_avg = sum(VA_t)/len(VA_t)

        # calculate the normalization
        Er_norm = E_r / (np.abs(VA_t) * np.abs(B0))
        Be_norm = B_e / (np.abs(B0))
        Ee_norm = E_e / (np.abs(VA_t) * np.abs(B0))
        Br_norm = B_r / (np.abs(B0))

        # calculate the FFTs

        # B_norm
        N, T = len(data_dicts[4]['Epoch'][0]), 1 / 128
        xf_B = fftfreq(N, T)[:N // 2]
        yf_Be = rfft(Be_norm)
        FFT_Be = 2.0 / N * np.abs(yf_Be[0:N // 2])
        yf_Br = rfft(Br_norm)
        FFT_Br = 2.0 / N * np.abs(yf_Br[0:N // 2])

        # deltaB
        N, T = len(data_dicts[0]['Epoch'][0]), 1 / 128
        xf_deltaB = fftfreq(N, T)[:N // 2]
        yf_dBe = rfft(dB_e)
        FFT_dBe = 2.0 / N * np.abs(yf_dBe[0:N // 2])
        yf_dBr = rfft(dB_r)
        FFT_dBr = 2.0 / N * np.abs(yf_dBr[0:N // 2])

        # E_norm
        N, T = len(data_dicts[1]['Epoch'][0]), 1 / 128
        xf_E = fftfreq(N, T)[:N // 2]
        yf_Ee = rfft(Ee_norm)
        FFT_Ee = 2.0 / N * np.abs(yf_Ee[0:N // 2])
        yf_Er = rfft(Er_norm)
        FFT_Er = 2.0 / N * np.abs(yf_Er[0:N // 2])

        # deltaE
        N, T = len(data_dicts[5]['Epoch'][0]), 1 / 128
        xf_deltaE = fftfreq(N, T)[:N // 2]
        yf_dEe = rfft(dE_e)
        FFT_dEe = 2.0 / N * np.abs(yf_dEe[0:N // 2])
        yf_dEr = rfft(dE_r)
        FFT_dEr = 2.0 / N * np.abs(yf_dEr[0:N // 2])

        # store everything
        data_dict_FFT['B0'].append(B0)
        data_dict_FFT['Be_fft'].append(FFT_dBe)
        data_dict_FFT['Br_fft'].append(FFT_dBr)
        data_dict_FFT['Ee_fft'].append(FFT_dEe)
        data_dict_FFT['Er_fft'].append(FFT_dEr)
        data_dict_FFT['xf_deltaB'].append(xf_deltaB)
        data_dict_FFT['xf_deltaE'].append(xf_deltaE)
        data_dict_FFT['Be_norm_fft'].append(FFT_Be)
        data_dict_FFT['Br_norm_fft'].append(FFT_Br)
        data_dict_FFT['Ee_norm_fft'].append(FFT_Ee)
        data_dict_FFT['Er_norm_fft'].append(FFT_Er)
        data_dict_FFT['xf_B'].append(xf_B)
        data_dict_FFT['xf_E'].append(xf_E)
        data_dict_FFT['ni'].append(ni)
        data_dict_FFT['VA_avg'].append(VA_avg)
        data_dict_FFT['VA_t'].append(VA_t)

        print(np.average(data_dicts[6]['Ion_Comp'][0]), np.average(data_dicts[6]['Op_Comp'][0]))

    # ---------------------
    #######################
    # --- PLOT THE DATA ---
    #######################
    # ---------------------
    # prepare the data for plotting
    waveSetKeys = [
        ['B_e', 'E_r'],
        ['B_r', 'E_e']]

    waveSetLabels = [['$\delta B_{e}$', '$\delta E_{r}$'],
                     [' $\delta B_{r}$', '-$\delta E_{e}$']]

    if NORMALIZE:
        spectraFreqs = [data_dict_FFT['xf_B'], data_dict_FFT['xf_E']]

        spectraData = [
                       [data_dict_FFT['Be_norm_fft'], data_dict_FFT['Er_norm_fft']],
                       [data_dict_FFT['Br_norm_fft'], data_dict_FFT['Ee_norm_fft']]
                       ]
        spectraLabels = [['$B_{e}(t)/|B_{0}(t)|$', '$E_{r}(t)/|B_{0}(t)||V_{A}(t)|$'],
                         ['$B_{r}(t)/|B_{0}(t)|$', '$E_{e}(t)/|B_{0}(t)||V_{A}(t)|$']]

        alfvenData = [ [ np.array(spectraData[0][1][i])/np.array(spectraData[0][0][i]) for i in range(len(targetVar))],
                       [ np.array(spectraData[1][1][i])/np.array(spectraData[1][0][i]) for i in range(len(targetVar))],
                       ]
    else:
        spectraFreqs = [data_dict_FFT['xf_deltaB'], data_dict_FFT['xf_deltaE']]
        spectraData = [  [data_dict_FFT['Be_fft'], data_dict_FFT['Er_fft']],
                         [data_dict_FFT['Br_fft'], data_dict_FFT['Ee_fft']]
                         ]
        spectraLabels = [['$\delta B_{e}$ [nT]', '$\delta E_{r}$ [mV/m]'],
                         [' $\delta B_{r}$ [nT]', '-$\delta E_{e}$ [mV/m]']]
        alfvenData = [[(1E6)*np.array(spectraData[0][1][i]) / np.array(spectraData[0][0][i]) for i in range(len(targetVar))],
                      [(1E6)*np.array(spectraData[1][1][i]) / np.array(spectraData[1][0][i]) for i in range(len(targetVar))],
                      ]

    #############################
    ### DO THE PLOTTING LOOPS ###
    #############################

    for wKeySet in kSets:

        stl.prgMsg('Plotting Data')
        fig, ax = plt.subplots(nrows=4, ncols=len(targetVar))
        fig.set_size_inches(figure_width, figure_height)

        wWaveSetKeys = waveSetKeys[wKeySet]
        wSpectraLabel = spectraLabels[wKeySet]
        wSpectraData = spectraData[wKeySet]
        wAlfvenData = alfvenData[wKeySet]
        wWaveSetLabels = waveSetLabels[wKeySet]


        ######## LOOP THROUGH ALL INVARIANT LATTITUDE REGIONS ########
        for i in range(len(targetVar)):

            data_dicts = dict_sets[i] # raw data

            # Set the title
            lowTime = sectionTimeRange[i][0]
            highTime = sectionTimeRange[i][1]

            # ax[0][i].set_title(RegionNames[i] + f'\n {lowTime.hour}:{lowTime.minute}:{lowTime.second + round(lowTime.microsecond/1E6,2)} to {highTime.hour}:{highTime.minute}:0{highTime.second + round(highTime.microsecond/1E6,2)} UTC',
            #                    weight='bold',
            #                    fontsize=PlotTitleSize,
            #                    pad=TitlePadding)

            ax[0][i].set_title(RegionNames[i] + f'\n {lowTime.strftime("%H:%M:%S")} to {highTime.strftime("%H:%M:%S")} UTC',
                               weight='bold',
                               fontsize=PlotTitleSize,
                               pad=TitlePadding)

            ####################
            # --- Everything ---
            #####################

            # set the grid
            for j in range(4):
                ax[j][i].grid(which='both', alpha=0.5)

            #################################
            # --- delta Be, delta Er plot ---
            #################################

            # delta Be
            ax[0][i].set_xlabel('ILat [deg]',fontsize=PlotLabelSize, labelpad=PlotLabelPad+1)
            ax[0][i].set_ylim(-8, 8)
            ax[0][i].set_xmargin(0)
            ax[0][i].tick_params(axis='both', which='major', labelsize=PlotTickLabelSize)
            ax[0][i].tick_params(axis='y', which='minor', labelsize=PlotTickLabelSize,length=0,grid_alpha=0)
            ax[0][i].tick_params(axis='both', which='minor', labelsize=PlotTickLabelSize-2)
            ax[0][i].plot(data_dicts[0]['ILat'][0], data_dicts[0][wWaveSetKeys[0]][0], color=plotColors[0], linewidth=plotLineWidth, label=f'{wWaveSetLabels[0]} [nT]')

            # delta Er
            mod = -1 if wKeySet == 1 else 1
            ax[0][i].plot(data_dicts[1]['ILat'][0], mod*data_dicts[1][wWaveSetKeys[1]][0], color=plotColors[1], linewidth=plotLineWidth,label=f"{wWaveSetLabels[1]} [mV/m]")
            ax[0][i].set_xmargin(0)
            ax[0][i].minorticks_on()

            # fix the xticks
            # newTicks = data_dicts[1]['ILat'][0][::int(len(data_dicts[1]['ILat'][0])/2)]
            # newTicks = [round(tick, 2) for tick in newTicks]
            # newTicks.append(round(data_dicts[1]['ILat'][0][-1], 2))
            # newTickStr = [str(tick) for tick in newTicks]
            # ax[0][i].set_xticks(newTicks, newTickStr)

            ############################
            # --- Poynting Flux plot ---
            ############################
            ax[1][i].plot(data_dicts[2]['ILat'][0], PoyntingScale*np.array(data_dicts[2]['S_p'][0]), plotColors[2], label='$\delta S_{p}$ [erg/$cm^{2}$s]',linewidth=plotLineWidth)
            ax[1][i].set_xlabel('ILat [deg]',fontsize=PlotLabelSize, labelpad=PlotLabelPad+1)
            ax[1][i].set_ylim(-0.005, 0.03)
            ax[1][i].set_xmargin(0)
            # ax[1][i].invert_yaxis()
            ax[1][i].tick_params(axis='both', which='major', labelsize=PlotTickLabelSize)
            ax[1][i].tick_params(axis='both', which='minor', labelsize=PlotTickLabelSize-2)
            # newTicks = data_dicts[1]['ILat'][0][:-1:int(len(data_dicts[1]['ILat'][0])/2)]
            # newTicks = [round(tick, 2) for tick in newTicks]
            # newTicks.append(round(data_dicts[1]['ILat'][0][-1], 2))
            # newTickStr = [str(tick) for tick in newTicks]
            # ax[1][i].set_xticks(newTicks, newTickStr)
            # ax[1][i].minorticks_on()

            ######################
            # --- FFT E, FFT B ---
            ######################

            ln1 = ax[2][i].plot(spectraFreqs[0][i], wSpectraData[0][i], color=plotColors[0], label=f'{wSpectraLabel[0]}',linewidth=plotLineWidth)
            ax[2][i].set_xlabel('Frequency [Hz]', fontsize=PlotLabelSize, labelpad=PlotLabelPad)
            ax[2][i].set_xlim(0, freqLimit)
            ax[2][i].set_yscale('log')
            ax[2][i].set_ylim(1E-3, 1E1)

            # if NORMALIZE:
            #     print('potato')
            #     ax[2][i].set_ylim(1E-7, 1E-2)
            #     # freqS = [val**(-1.67) for val in spectraFreqs[0][i] if val !=0]
            #     # ax[2][i].plot(spectraFreqs[0][i][:-1],freqS,color='black',linestyle='--',label='$f^{-1.67}$')
            # else:
            #     # freqS = [val ** (-1.67) for val in spectraFreqs[0][i] if val != 0]
            #     # ax[2][i].plot(spectraFreqs[0][i][:-1], freqS, color='black', linestyle='--')
            #     # ax[2][i].axvline(0.8, color='green',linestyle='--',alpha=0.5)
            #
            #     # for k in range(3):
            #     #     ax[2][i].axvline(0.55*(k+1), color='green',linestyle='--',alpha=0.5)

            ax[2][i].tick_params(axis='both', which='major', labelsize=PlotTickLabelSize)
            ax[2][i].tick_params(axis='both', which='minor', labelsize=PlotTickLabelSize-2)

            axFFT_E = ax[2][i].twinx()
            ln2 = ax[2][i].plot(spectraFreqs[1][i], wSpectraData[1][i], color=plotColors[1], label=f'{wSpectraLabel[1]}',linewidth=plotLineWidth)

            axFFT_E.set_yscale('log')
            axFFT_E.set_ylim(1E-3, 1E1)
            axFFT_E.set_yticks([])
            ax[2][i].minorticks_on()

            ######################
            # --- ASpeed Ratio ---
            ######################
            ax[3][i].plot(spectraFreqs[0][i], wAlfvenData[i], color='black', label=f'', linewidth=0.75, linestyle='-', marker='o', markersize=1)
            # ax[3][i].plot(spectraFreqs[0][i], wAlfvenData[i], color='black', label=f'', linewidth=1.2, linestyle='-')

            if NORMALIZE:
                ax[3][i].axhline(1, linestyle='--', color='red', label='V$_{A}$ (Avg)')
                ax[3][i].set_ylim(0, 3)
            else:
                ax[3][i].axhline(data_dict_FFT['VA_avg'][0],linestyle='--', color='red', label='V$_{A}$ (Avg) [m/s]')
                ax[3][i].set_ylim(EB_ratio_limits[0], EB_ratio_limits[1])
                ax[3][i].set_yscale('log')
                # ax[3][i].axvline(0.8, color='green', linestyle='--', alpha=0.5)

                # for k in range(3):
                #     ax[3][i].axvline(0.55*(k+1), color='green',linestyle='--',alpha=0.5)

            ax[3][i].minorticks_on()
            ax[3][i].set_xlabel('Frequency [Hz]', fontsize=PlotLabelSize, labelpad=PlotLabelPad)
            ax[3][i].set_xlim(0, freqLimit)
            ax[3][i].tick_params(axis='both', which='major', labelsize=PlotTickLabelSize)
            ax[3][i].tick_params(axis='both', which='minor', labelsize=PlotTickLabelSize-2)


            #########################
            # --- PHASE HISTOGRAM ---
            #########################
            #
            # # B-Field
            # f_B, t_B, Sxx_B = spectrogram(data_dicts[0][wWaveSetKeys[0]][0],
            #                               fs=128,
            #                               window=windowType,
            #                               nperseg=npersegN,  # note: if ==None default size is 256
            #                               noverlap=overlap,
            #                               scaling=scalingType,
            #                               mode='angle')  # scaling = density or scaling = spectrum
            # # E-Field
            # f_E, t_E, Sxx_E = spectrogram(mod * data_dicts[1][wWaveSetKeys[1]][0],
            #                               fs=128,
            #                               window=windowType,
            #                               nperseg=npersegN,  # note: if ==None default size is 256
            #                               noverlap=overlap,
            #                               scaling=scalingType,
            #                               mode='angle')  # scaling = density or scaling = spectrum
            #
            # phaseDiff = np.degrees(Sxx_B) - np.degrees(Sxx_E)
            #
            # # --- Create the Histogram ---
            # phaseDiffBins = [-187.5 + 15 * i for i in range(26)]
            # plotBins = [-180 + 15 * i for i in range(25)]
            # phaseHistogram = np.zeros(shape=(len(f_B), len(phaseDiffBins) - 1))
            # medianPhase = []
            #
            # for frq in range(len(phaseDiff)):
            #     hist, bin_edges = np.histogram(phaseDiff[frq], bins=phaseDiffBins)
            #     phaseHistogram[frq] = hist
            #
            #     # Get the median value for each frequency
            #     # if i == 2:
            #     #     print(f_B[frq],np.median(phaseDiff[frq]))
            #
            #     medianPhase.append(np.median(phaseDiff[frq]))
            #
            #
            #
            # # Normalize the histogram
            # phaseHistogram = 100 * phaseHistogram / phaseHistogram.max()
            # cmap = ax[4][i].pcolormesh(f_B, plotBins, phaseHistogram.T, vmin=0, vmax=50, cmap=matlab_parula_cmap())
            # ax[4][i].plot(f_B,medianPhase, color='white', linewidth=plotLineWidth+0.5)
            #
            #
            # ax[4][i].set_xlabel('Frequency [Hz]', fontsize=PlotLabelSize)
            # ax[4][i].set_xlim(0, freqLimit)
            # ax[4][i].set_yticks([-180, -90, 0, 90, 180])
            # ax[4][i].minorticks_on()
            # ax[4][i].tick_params(axis='both', which='major', labelsize=PlotTickLabelSize)
            # ax[4][i].tick_params(axis='both', which='minor', labelsize=PlotTickLabelSize - 2)
            # ax[4][i].set_ylim(-180, 180)

            # for k in range(3):
            #     ax[4][i].axvline(0.55 * (k + 1), color='green', linestyle='--', alpha=0.5)


            # set the legend(s)
            lns = ln1 + ln2
            labs = [l.get_label() for l in lns]
            if i == 0:
                ax[0][i].set_ylabel(f'{wWaveSetLabels[0]} and {wWaveSetLabels[1]}' + "\n 0.7 Hz HighPass", fontsize=PlotLabelSize,labelpad=PlotLabelPad)
                ax[0][i].legend(loc='upper right', prop={'size': PlotLegendSize})
                ax[1][i].set_ylabel('Poynting Flux', fontsize=PlotLabelSize,labelpad=PlotLabelPad+2)
                ax[1][i].legend(loc='upper right', prop={'size': PlotLegendSize})
                ax[2][i].set_ylabel('Amplitude Spectra', fontsize=PlotLabelSize,labelpad=PlotLabelPad+2)
                axFFT_E.legend(lns, labs, loc='upper right', prop={'size': PlotLegendSize})
                ax[3][i].legend(loc='upper right', prop={'size': PlotLegendSize})
                # ax[4][i].set_ylabel('Phase Difference', fontsize=PlotLabelSize,labelpad=PlotLabelPad-1)

                if NORMALIZE:
                    ax[3][i].set_ylabel('$E_{\perp}(t)/B_{\perp}(t)$', fontsize=PlotLabelSize)
                else:
                    ax[3][i].set_ylabel(f'{wWaveSetLabels[1]}/{wWaveSetLabels[0]} [m/s]', fontsize=PlotLabelSize)
            else:
                ax[0][i].set_yticklabels([])
                ax[1][i].set_yticklabels([])
                ax[2][i].set_yticklabels([])
                ax[3][i].set_yticklabels([])



        ### OUTPUT THE DATA ###
        # plt.tight_layout(rect=[0, 0.0, 0.93, 0.99], pad=-0.2)
        # plt.tight_layout()
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.3)
        fileOutName = "Plot3_WaveAnalysis_BeEr_base.png" if wWaveSetKeys[0] == 'B_e' and wWaveSetKeys[1] =='E_r' else "Plot3_WaveAnalysis_mBrEe_base.png"
        outputPath = rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot3\{fileOutName}'

        # spectrogram colorbar
        # cax = fig.add_axes([0.932, 0.0325, 0.0155, 0.158])
        # cbar = plt.colorbar(cmap, cax=cax)
        # cbar.ax.minorticks_on()
        # cbar.set_label('Histogram Occurance (%)', labelpad=5)

        plt.savefig(outputPath, dpi=dpi,bbox_inches='tight')
        stl.Done(start_time)









#################
# --- EXECUTE ---
#################
Plot3_deltaEdeltaB_waveInfo(targetVar,dict_sets)
