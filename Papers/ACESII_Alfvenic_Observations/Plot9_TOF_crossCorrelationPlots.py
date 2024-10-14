# --- crossCorrelationTOFAnalysis.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Isolate Alfvenic Dispersions and perform cross-correlation analysis
# to determine values that are fitted linearly to get the height of the source region

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
from ACESII_code.myImports import *
start_time = time.time()

# --- --- --- --- ---
justPrintFileNames = False
wRocket = 4
modifier = ''
inputPath_modifier = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'science\AlfvenSignatureAnalysis' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder


################################
# --- Plot toggles - General ---
################################
figure_width = 12 # in inches
figure_height = 8.5 # in inches
Title_FontSize = 20
Label_FontSize = 15
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
Legend_fontSize = 15
dpi = 200


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
# plot all of the dispersion functions over a range of pitch angles (user input)
wDispersions = [2,3,4,5] # [] -> plot all dispersion traces, [#,#,#,...] plot specific ones. USE THE DISPERSION NUMBER NOT PYTHON -1 INDEX
wPitch = 2 # plots specific pitch angles by their index
# ---------------------------
justPlotKeyDispersions = False #IF ==TRUE no cross-correlation will occur
applyMaskVal = True
cbar_low, cbar_high = 0, 30
maskVal = 5 # apply a single mask to the dispersion feature
# ---------------------------
isolateAlfvenSignature = True # removes unwanted data from the alfven signature
# ---------------------------
DetectorTimeResolution = 0.05 # in seconds
DetectorEnergyResolution = 0.18
# ---------------------------
correlationAnalysis = True
showErrorBars = False
weightLinearFitByCounts = False
outputCorrelationPlot = True
# ---------------------------
LambdaPerpPlot = False
LambdaPerpFit = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---

from scipy.signal import correlate, correlation_lags
from itertools import combinations
from ACESII_code.Science.AlfvenSingatureAnalysis.Particles.dispersionAttributes import dispersionAttributes
from myspaceToolsLib.physicsVariables import Re

mycmap = apl_rainbow_black0_cmap()



# Linear Fitted Model
def fitFunc_linear(x, a, b):
    return a * x + b

# Polynomial Fitted Model
# def fitFunc_polynomial(x, a0, a1, a2):
#     return a0 + a1 * x + a2 * (x ** 2)

def fitFunc_polynomial(x,a0, a1, a2):
    return  a0+ a1 * x + a2 * (x ** 2)


def AlfvenSignatureCrossCorrelation(wRocket, rocketFolderPath, justPrintFileNames,wDis):

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}{modifier}\*eepaa_fullcal*')
    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}{modifier}\\', '') for ifile in inputFiles]
    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return

    # --- get the data from the ESA file ---
    data_dict = loadDictFromFile(inputFiles[0])

    # --- --- --- --- --- --- ----
    # --- Isolate a Dispersion ---
    # --- --- --- --- --- --- ----

    wDispersion_key = f's{wDis}'
    Energy = data_dict['Energy'][0]
    Pitch = data_dict['Pitch_Angle'][0]
    lowCut, highCut = np.abs(data_dict['Epoch'][0] - dispersionAttributes.keyDispersionDeltaT[wDis-1][0]).argmin(), np.abs(data_dict['Epoch'][0] - dispersionAttributes.keyDispersionDeltaT[wDis-1][1]).argmin()
    Epoch_dis = dateTimetoTT2000(data_dict['Epoch'][0][lowCut:highCut+1],inverse=False)
    Epoch_dis = (np.array(Epoch_dis) - Epoch_dis[0]) / 1E9 # converted data to TIME SINCE START OF DISPERSION (needed for proper plot fitting)

    # calculate the center point (in time) of the dispersion
    whenSTEBoccured_time = data_dict['Epoch'][0][int((highCut+lowCut)/2)]
    whenSTEBoccured_Alt = data_dict['Alt'][0][int((highCut + lowCut) / 2)]

    eepaa_dis = data_dict['eepaa'][0][lowCut:highCut+1]

    if applyMaskVal:
        # --- Apply the Masking Value to the dataset ---
        # Note: This removes all Fillvals
        eepaa_dis[eepaa_dis<=maskVal] = 0 # applies mask and zeros-out anything below 0
        eepaa_dis = eepaa_dis.clip(min=0)
        eepaa_dis[eepaa_dis>1E15] = 0 # zeros-out anything above 1E15

    if isolateAlfvenSignature:
        eepaa_dis = np.array(dispersionAttributes.isolationFunctions[wDispersion_key](eepaa_dis, Energy, Epoch_dis)) # pply the isolation functions found in dispersionAttributes.py





    # --- Reduce data to Pitch Slice of Interest ---
    # gets one slice in pitch of the data and flips
    # data from dimensions 0 == Time to --> dimension 0 == Energy so that the first index gets energy for all time
    eepaa_dis_onePitch = np.array(eepaa_dis[:,wPitch,:]).T

    #################################################
    # --- PLOT THE DISPERSION IN ORDER TO ISOLATE ---
    #################################################

    if justPlotKeyDispersions:
        fig, ax = plt.subplots(2)
        figure_width = 25
        figure_height = 15
        fig.set_figwidth(figure_width)
        fig.set_figheight(figure_height)
        fig.suptitle(f'STEB {wDis}\n Pitch = {Pitch[wPitch]}deg\n {data_dict["Epoch"][0][lowCut].strftime("%H:%M:%S.%f")} to {data_dict["Epoch"][0][highCut].strftime("%H:%M:%S.%f")}')
        cmap = ax[0].pcolormesh(Epoch_dis, Energy, np.array(eepaa_dis[:, wPitch, :]).T, cmap=mycmap, vmin=cbar_low, vmax=cbar_high)
        ax[0].set_yscale('log')
        ax[0].set_ylim(Energy[-1], Energy[np.abs(Energy - dispersionAttributes.keyDispersionEnergyLimits[wDis - 1][1]).argmin() - 1])
        ax[1].pcolormesh(Epoch_dis, Energy, np.array(data_dict['eepaa'][0][lowCut:highCut + 1, wPitch, :]).T,cmap=mycmap, vmin=cbar_low, vmax=cbar_high)
        ax[1].set_yscale('log')
        ax[1].set_ylim(Energy[-1], Energy[np.abs(Energy - dispersionAttributes.keyDispersionEnergyLimits[wDis - 1][1]).argmin() - 1])
        cbar = plt.colorbar(cmap, ax=ax.ravel().tolist())
        plt.show()
    else:
        ############################################
        # --- PERFORM CROSS CORRELATION ANALYSIS ---
        ############################################
        peakVal = 0

        for engy in range(len(Energy)):
            for tme in range(len(eepaa_dis_onePitch.T)):
                if eepaa_dis_onePitch[engy][tme] > peakVal:
                    peakVal = eepaa_dis_onePitch[engy][tme]
                    engyIndx = engy
                    tmeIndex = tme
                    # print(wDis, Energy[engyIndx], Epoch_dis[tmeIndex], peakVal)


        # for each energy in the eepaa data, check if there's at least 1 count in that energy "row". If so, consider this energy
        validEngyIndicies = []
        for engy,engyRow in enumerate(eepaa_dis_onePitch):
            aboveZero = np.where(engyRow > 0)[0]

            if aboveZero.size != 0:
                validEngyIndicies.append(engy)

        validEngyIndicies = np.array(validEngyIndicies)
        minEnergy = validEngyIndicies.max() # remember: LOWER index --> Higher energy
        # energyPairs = [comb for comb in combinations(validEngyIndicies,2) if minEnergy in comb]

        maxEnergy = validEngyIndicies.min()
        energyPairs = [comb for comb in combinations(validEngyIndicies, 2) if maxEnergy in comb]

        # for each pair of Energies, perform the cross-correlation analysis:
        # Note: the peak value in the lag-time determines the deltaT value between the velocities
        deltaTs = []
        deltaVs = []
        errorT = []
        errorV = []
        errorZ = []

        for engyPair in energyPairs:

            higherEnergyIndex = engyPair[0]
            lowerEnergyIndex = engyPair[1]
            EnergySlow = q0*Energy[lowerEnergyIndex]
            EnergyFast = q0*Energy[higherEnergyIndex]
            velSlow = np.cos(np.radians(Pitch[wPitch]))*np.sqrt(2*EnergySlow/m_e)
            errorVelSlow = np.sqrt(DetectorEnergyResolution) * velSlow
            velFast = np.cos(np.radians(Pitch[wPitch]))*np.sqrt(2*EnergyFast/m_e)
            errorVelFast = np.sqrt(DetectorEnergyResolution) * velFast

            # check to make sure both datasets have at least one non-zero value
            if not np.any(eepaa_dis_onePitch[higherEnergyIndex]) or not np.any(eepaa_dis_onePitch[lowerEnergyIndex]):
                continue
            else:
                corr = np.array(correlate(eepaa_dis_onePitch[lowerEnergyIndex],eepaa_dis_onePitch[higherEnergyIndex]))
                lags = correlation_lags(len(eepaa_dis_onePitch[higherEnergyIndex]), len(eepaa_dis_onePitch[lowerEnergyIndex]))

            try:
                # Find the x,y value of the peak in the correlation output
                indexMax = np.where(corr == np.max(corr))[0][0]
                delayTime = DetectorTimeResolution*lags[indexMax]
                velDiff = 1000*(np.sqrt(m_e) / (np.cos(np.radians(Pitch[wPitch]))*(np.sqrt(2)))) * (1/(np.sqrt(EnergySlow)) - 1/(np.sqrt(EnergyFast)))

                if weightLinearFitByCounts:
                    # find if counts_slow or counts_fast is more at the correlation peak
                    N_fast = np.array(eepaa_dis_onePitch[higherEnergyIndex]).max()
                    N_Slow = np.array(eepaa_dis_onePitch[lowerEnergyIndex]).max()

                    # whichever is less (N_Fast,N_slow) append that many values to the results in order to bias the linear fit
                    # towards values that are both higher in counts, instead of just one
                    for h in range(min(N_fast, N_Slow)):
                        deltaTs.append(delayTime)
                        deltaVs.append(velDiff)
                else:
                    deltaTs.append(delayTime)
                    deltaVs.append(velDiff)

                # calculate the error in the measurements
                errorT.append(DetectorTimeResolution)
                errorV.append(((DetectorEnergyResolution*np.sqrt(m_e)) / ((2**(3/2))*np.cos(np.radians(Pitch[wPitch])))) * np.sqrt(1/EnergySlow + 1/EnergyFast))

                errZ1 = (DetectorTimeResolution**2) * (1 / (1/velSlow - 1/velFast))
                errZ2 = (np.square(errorVelSlow))*(delayTime*np.square(velFast)/np.square(velSlow - velFast))**2
                errZ3 = (np.square(errorVelFast))*(delayTime*np.square(velSlow)/np.square(velSlow - velFast))**2
                errorZ.append(errZ1 + errZ2 + errZ3)

            except:
                print('indexMax = np.where(corr == np.max(corr))[0][0] Found no maximum. Correlation data is likely all zeros. Did you set your maskval too high?')

        deltaTs = np.array(deltaTs)
        deltaVs = np.array(deltaVs)

        # convert error in velocity to kilometers
        errorV = np.array(errorV)*1000

        ##########################################
        # --- FIT THE DATA TO DERIVE TOF Z_ACC ---
        ##########################################
        p0 = [-9.113E-1, 7.350E3, -1.108E7]
        paramsPoly, covPoly = curve_fit(fitFunc_polynomial, deltaVs, deltaTs, p0=p0, maxfev=int(1E6))
        paramsLin, covLin = curve_fit(fitFunc_linear, deltaVs, deltaTs)

        # --- CORRELATION ---
        # # linear

        fitData = fitFunc_linear(deltaVs, *paramsLin)
        r_corr_linear = np.corrcoef([deltaTs,fitData])[0][1]

        # Polynomial
        fitData_poly = fitFunc_polynomial(deltaVs, *paramsPoly)
        r_corr_poly = np.corrcoef([deltaTs,fitData_poly])[0][1]


        # calculate chisquare
        # varianceT = np.var(deltaTs)
        # varianceV = np.var(deltaVs)
        # variance = varianceT + varianceV
        # ExpectedLin = fitFunc_linear(deltaVs,*paramsLin)
        # ExpectedPoly = fitFunc_polynomial(deltaVs,*paramsPoly)
        # print('\n')
        # print('lin',ExpectedLin)
        # print('poly',ExpectedPoly)
        #
        # ChiLin =  sum( [ (ExpectedLin[i] - deltaTs[i])**2/(0.05**2) for i in range(len(deltaVs))])
        # ChiPoly = sum([ (ExpectedPoly[i] - deltaTs[i])**2/(0.05**2)   for i in range(len(deltaVs))])
        #
        # print('ChiLin', ChiLin)
        # print('ChiPoly', ChiPoly)
        #
        # fig, ax = plt.subplots()
        # ax.plot(deltaVs,deltaTs,color='blue')
        # ax.plot(deltaVs,ExpectedPoly,color='red')
        # plt.show()

        ########################################
        # --- Return results of the Analysis ---
        ########################################
        # return an array of: [STEB Number,Observation Time, Observation Altitude, Z_acc, Z_acc_error, correlationR]
        errorZ_avg = sum(errorZ) / len(errorZ)
        return [wDis, whenSTEBoccured_time, whenSTEBoccured_Alt, paramsLin, errorZ_avg, r_corr_linear, paramsPoly, deltaVs, deltaTs,r_corr_poly]


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
STEBfitResults = []
for wDis in tqdm(wDispersions):
    results = AlfvenSignatureCrossCorrelation(wRocket, rocketFolderPath, justPrintFileNames, wDis)
    STEBfitResults.append(results)




if not justPrintFileNames:
    ##################
    # --- PLOTTING ---
    ##################
    fig, ax = plt.subplots(nrows=2,ncols=2)
    fig.set_size_inches(figure_width, figure_height)
    counter = 0
    Pitch = [-10 + i*10 for i in range(21)]

    for rowIdx in range(2):
        for colIdx in range(2):
            wDis = STEBfitResults[counter][0]
            deltaVs = STEBfitResults[counter][7]
            deltaTs = STEBfitResults[counter][8]
            paramsLin = STEBfitResults[counter][3]
            paramsPoly = STEBfitResults[counter][6]
            r_corr_linear = STEBfitResults[counter][5]
            r_corr_poly = STEBfitResults[counter][-1]


            x_s = np.linspace(deltaVs.min(), deltaVs.max(), 20)
            fitData = fitFunc_linear(x_s, *paramsLin)
            fitData_poly = fitFunc_polynomial(x_s, *paramsPoly)

            ax[rowIdx, colIdx].set_title(f'S{wDis}', fontsize=Title_FontSize)


            # Raw Data
            ax[rowIdx, colIdx].scatter(deltaVs, deltaTs)

            if colIdx == 0:
                ax[rowIdx, colIdx].set_ylabel('Delay time [s]',fontsize=Label_FontSize)
            if rowIdx == 1:
                ax[rowIdx, colIdx].set_xlabel('1/v [s/km]',fontsize=Label_FontSize)
            ax[rowIdx, colIdx].plot(x_s, fitData, color="red", label=f'r-corr. (first order) = {round(r_corr_linear, 3)}')
            ax[rowIdx, colIdx].plot(x_s, fitData_poly, color='black', label=f'r-corr. (second order) = {round(r_corr_poly,3)}')
            # ax.plot(x_s, fitData_poly, color='black', label=f'a0 = {paramsPoly[0]}\n a1 = {paramsPoly[1]} \n a2 = {0}\n r_corr_poly = {r_corr_poly}')

            xticks = np.linspace(0, deltaVs.max(), 6)
            xtick_labels = [f'{x:.2e}' for x in xticks]
            xtick_labels[0] = '0'
            ax[rowIdx, colIdx].set_xticks(xticks, xtick_labels)
            ax[rowIdx, colIdx].legend(loc='lower right', fontsize=Legend_fontSize)
            counter += 1

    plt.tight_layout()
    plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot9\Plot9_TOFfunc_base.png', dpi=dpi)






if LambdaPerpFit:

    Tanaka_XDat = [1,3,9]
    Tanaka_YDat = [-2.77E7,-9.71E6,-8.29E6]


    def fitFunc(x, a0, a2):
        return a0*np.log(x+1) + a2

    guess = [1E6,-1.04E7]
    paramsFit,cov = curve_fit(fitFunc,Tanaka_XDat,Tanaka_YDat, p0=guess)


    x_s = np.linspace(min(Tanaka_XDat),max(Tanaka_XDat),200)
    fitDat = fitFunc(x_s,*paramsFit)
    fig,ax = plt.subplots()
    ax.scatter(Tanaka_XDat,Tanaka_YDat)
    ax.plot(x_s,fitDat)
    plt.show()

    def inverseFunc(x,a,c):
        return np.exp((x-c)/a)+1
    print(paramsFit)

    for data in STEBfitResults:
        if data[0] in [2, 3, 4, 5]:
            params = data[6]
            deltaV = data[7]
            deltaT = data[8]

            a2_param = params[2]

            print(data[0],a2_param,inverseFunc(a2_param,*paramsFit))

if LambdaPerpPlot:

    fig, ax = plt.subplots()

    for data in STEBfitResults:
        if data[0] in [2, 3, 4, 5]:
            params = data[6]
            deltaV = data[7]
            deltaT = data[8]

            x_s = np.linspace(deltaV.min(), deltaV.max(), 20)
            # x_s = np.linspace(deltaV.min(), 3E-4, 20)

            def fitFunc_polynomial(x, a0, a1, a2):
                return a0 + a1 * x + a2 * (x ** 2)
            fitData_poly = fitFunc_polynomial(x_s, *params)

            ax.plot(x_s, fitData_poly,label=f'STEB {data[0]}')

    ax.set_ylabel('Delay time [s]')
    ax.set_xlabel('1/v [s/km]')
    ax.set_xlim(0,3E-4)
    plt.legend()
    plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\TOFanalysis\STEB_LambdaPerp.png')




