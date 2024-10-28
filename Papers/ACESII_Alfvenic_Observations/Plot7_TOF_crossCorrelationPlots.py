# --- crossCorrelationTOFAnalysis.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Isolate Alfvenic Dispersions and perform cross-correlation analysis
# to determine values that are fitted linearly to get the height of the source region

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import spaceToolsLib as stl
from myImports import *
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
figure_height = 3*4.25 # in inches
Title_FontSize = 20
Label_FontSize = 20
Label_Padding = 8

Text_Fontsize = 25

Tick_FontSize = 15
Tick_FontSize_minor = 20
Tick_Length = 6
Tick_Width = 2

Plot_LineWidth = 0.5
Plot_MarkerSize = 14
Legend_fontSize = 15
dpi = 200


cbar_low, cbar_high = 0, 40
cbar_fontSize = 15

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
# plot all of the dispersion functions over a range of pitch angles (user input)
wDispersions = [2, 3, 4, 5] # [] -> plot all dispersion traces, [#,#,#,...] plot specific ones. USE THE DISPERSION NUMBER NOT PYTHON -1 INDEX
wPitch = 2 # plots specific pitch angles by their index
# ---------------------------
justPlotKeyDispersions = False #IF ==TRUE no cross-correlation will occur
applyMaskVal = True
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
from Science.AlfvenSingatureAnalysis.Particles.dispersionAttributes import dispersionAttributes

mycmap = stl.apl_rainbow_black0_cmap()



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
    Epoch_dis = stl.dateTimetoTT2000(data_dict['Epoch'][0][lowCut:highCut+1],inverse=False)
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
        cmap = ax[0].pcolormesh(Epoch_dis, Energy, eepaa_dis_onePitch, cmap=mycmap, vmin=cbar_low, vmax=cbar_high)
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

        ########################################
        # --- Return results of the Analysis ---
        ########################################
        # return an array of: [STEB Number,Observation Time, Observation Altitude, Z_acc, Z_acc_error, correlationR]
        errorZ_avg = sum(errorZ) / len(errorZ)
        return [wDis, whenSTEBoccured_time, whenSTEBoccured_Alt, paramsLin, errorZ_avg, r_corr_linear, paramsPoly, deltaVs, deltaTs, r_corr_poly, eepaa_dis_onePitch]


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
STEBfitResults = []
for wDis in tqdm(wDispersions):
    results = AlfvenSignatureCrossCorrelation(wRocket, rocketFolderPath, justPrintFileNames, wDis)
    STEBfitResults.append(results)

for thing in STEBfitResults:
    print('\n')
    print(f'Dispersion No: {thing[0]}')
    print(f'Params (Linear): {thing[3]}')
    print(f'R-corr (Linear): {thing[5]}')
    print(f'Params (Poly): {thing[6]}')
    print(f'R-corr (Poly): {thing[9]}')
    print('\n')



if not justPrintFileNames:
    ##################
    # --- PLOTTING ---
    ##################
    fig, ax = plt.subplots(nrows=4,ncols=2,height_ratios=[1,0.1,1,1])
    fig.set_size_inches(figure_width, figure_height)
    counter = 0
    Pitch = [-10 + i*10 for i in range(21)]

    # === Plot the example data for S2 and S5 ===
    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket - 4]}{modifier}\*eepaa_fullcal*')
    data_dict = loadDictFromFile(inputFiles[0])
    Energy = data_dict['Energy'][0]
    Pitch = data_dict['Pitch_Angle'][0]

    wDisKeys = ['s2', 's5']
    wDis_indicies = [2, 5]

    for k, wDispersion_key in enumerate(wDisKeys):

        lowCut, highCut = np.abs(data_dict['Epoch'][0] - dispersionAttributes.keyDispersionDeltaT[wDis_indicies[k] - 1][0]).argmin(), np.abs(data_dict['Epoch'][0] - dispersionAttributes.keyDispersionDeltaT[wDis_indicies[k] - 1][1]).argmin()
        Epoch_dis = deepcopy(stl.dateTimetoTT2000(data_dict['Epoch'][0][lowCut:highCut + 1], inverse=False))
        Epoch_dis = (np.array(Epoch_dis) - Epoch_dis[0]) / 1E9  # converted data to TIME SINCE START OF DISPERSION (needed for proper plot fitting)

        # eepaa_dis = deepcopy(data_dict['eepaa'][0][lowCut:highCut + 1])
        if wDis_indicies[k]==2:
            eepaa_dis = STEBfitResults[0][-1]
        else:
            eepaa_dis = STEBfitResults[-1][-1]

        # --- Plot everything ---
        cmap = ax[0, k].pcolormesh(Epoch_dis, Energy, eepaa_dis, cmap=mycmap, vmin=cbar_low, vmax=cbar_high)
        ax[0, k].set_yscale('log')
        props = dict(boxstyle='round', facecolor='white', alpha=1, lw=4)
        ax[0, k].text((Epoch_dis[-1] - Epoch_dis[0])/2, 900, wDispersion_key.upper(), fontsize=Text_Fontsize, weight='bold', color='black', bbox=props, ha='center')
        ax[0, k].set_xlabel('Seconds [s]',fontsize=Label_FontSize-8)
        ax[0, k].tick_params(axis='y', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
        ax[0, k].tick_params(axis='y', which='minor', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length * 0.6)
        ax[0, k].tick_params(axis='x', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
        ax[0, k].tick_params(axis='x', which='minor', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length * 0.6)
        ax[0, k].set_ylim(30, 1000)

        if k ==0:
            ax[0, k].set_ylabel(r'$\alpha = 10^{\circ}$' + '\nEnergy [eV]',fontsize=Label_FontSize,weight='bold')
        if k == 1:
            ax[0, k].tick_params(axis='y', which='both',labelleft=False)


    # turn off the intentially blank plot
    ax[1, 0].axis('off')
    ax[1, 1].axis('off')

    # === Plot fitted TOF data ===

    for rowIdx in [2,3]:
        for colIdx in range(2):
            wDis = STEBfitResults[counter][0]
            deltaVs = STEBfitResults[counter][7]
            deltaTs = STEBfitResults[counter][8]
            paramsLin = STEBfitResults[counter][3]
            paramsPoly = STEBfitResults[counter][6]
            r_corr_linear = STEBfitResults[counter][5]
            r_corr_poly = STEBfitResults[counter][9]

            x_s = np.linspace(deltaVs.min(), deltaVs.max(), 20)
            fitData = fitFunc_linear(x_s, *paramsLin)
            fitData_poly = fitFunc_polynomial(x_s, *paramsPoly)

            props = dict(boxstyle='round', facecolor='white', alpha=1, lw=4)
            ax[rowIdx, colIdx].text(0.25E-4, 0.62, f'S{wDis}', fontsize=Title_FontSize, weight='bold', color='black', bbox=props, ha='center')


            # Raw Data
            ax[rowIdx, colIdx].scatter(deltaVs, deltaTs)
            ax[rowIdx, colIdx].tick_params(axis='both', which='both', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)

            if colIdx == 0:
                ax[rowIdx, colIdx].set_ylabel('Delay time [s]',fontsize=Label_FontSize,weight='bold')

            if colIdx == 1:
                ax[rowIdx, colIdx].tick_params(axis='y', which='both', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length, labelleft=False)

            if rowIdx == 2:
                ax[rowIdx, colIdx].tick_params(axis='x', which='both', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length,labelbottom=False)

            if rowIdx == 3:
                ax[rowIdx, colIdx].set_xlabel('1/v [s/km]',fontsize=Label_FontSize,weight='bold')

            ax[rowIdx, colIdx].plot(x_s, fitData, color="red", label=f'r-corr. (first order) = {round(r_corr_linear, 3)}')
            ax[rowIdx, colIdx].plot(x_s, fitData_poly, color='black', label=f'r-corr. (second order) = {round(r_corr_poly,3)}')
            ax[rowIdx, colIdx].set_xlim(0, 2.8E-4)
            ax[rowIdx, colIdx].set_ylim(-0.1, 0.71)
            ax[rowIdx, colIdx].legend(loc='lower right', fontsize=Legend_fontSize)
            ax[rowIdx, colIdx].ticklabel_format(axis='x',style='sci',scilimits=(0,0))
            counter += 1


    ##################
    # --- Colorbar ---
    ##################
    cax = fig.add_axes([0.89, 0.691, 0.025, 0.284])
    cbar = plt.colorbar(cmap, cax=cax)
    cbar.ax.minorticks_on()
    cbar.ax.tick_params(labelsize=Tick_FontSize, length=Tick_Length)

    cbar.ax.get_yaxis().labelpad = 17
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_weight("bold")
        l.set_fontsize(cbar_fontSize)

    cbar.set_label('Counts',rotation=90,fontsize=cbar_fontSize+10,weight='bold')

    plt.subplots_adjust(left=0.1, bottom=0.06, right=0.87, top=0.975, wspace=0.05, hspace=0.05)
    fig.align_ylabels(ax[:])
    plt.savefig(rf'C:\Users\cfelt\Desktop\Research\Feltman2024_ACESII_Alfven_Observations\PLOTS\Plot7\Plot7_TOFfunc_base.png', dpi=dpi)




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
    plt.savefig(rf'C:\Users\cfelt\Desktop\Research\Feltman2024_ACESII_Alfven_Observations\PLOTS\Plot7\STEB_LambdaPerp.png')




