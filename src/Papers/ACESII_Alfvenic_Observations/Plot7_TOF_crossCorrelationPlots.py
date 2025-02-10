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
from src.myImports import *
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
figure_height = 15 # in inches
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
mycmap = stl.apl_rainbow_black0_cmap()
mycmap.set_bad(color=(0, 0, 0))

cbar_low, cbar_high = 0, 45
cbar_fontSize = 20

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
# plot all of the dispersion functions over a range of pitch angles (user input)
wDispersions = [2,3,4,5] # [] -> plot all dispersion traces, [#,#,#,...] plot specific ones. USE THE DISPERSION NUMBER NOT PYTHON -1 INDEX
wPitch = 2 # plots specific pitch angles by their index
# ---------------------------
justPlotKeyDispersions = False #IF ==TRUE no cross-correlation will occur
applyMaskVal = True
maskVal = 3 # apply a single mask to the dispersion feature
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
generateLatexTable = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---

from scipy.signal import correlate, correlation_lags
from itertools import combinations
from src.Science.AlfvenSingatureAnalysis.Particles.dispersionAttributes import dispersionAttributes
import matplotlib.gridspec as gridspec




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
        # print(paramsLin)
        # print(paramsPoly)

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

if not justPrintFileNames:

    ##################
    # --- PLOTTING ---
    ##################

    fig = plt.figure()
    fig.set_size_inches(figure_width, figure_height)
    gs0 = gridspec.GridSpec(3,1, figure=fig, height_ratios=[1, 0.05, 1])

    gs00 = gridspec.GridSpecFromSubplotSpec(2,2,subplot_spec=gs0[0],hspace=0.16)
    gsNull = gridspec.GridSpecFromSubplotSpec(1,1,subplot_spec=gs0[1],hspace=0.0)
    gs01 = gridspec.GridSpecFromSubplotSpec(2,2,subplot_spec=gs0[2],hspace=0.05)

    counter = 0

    # === Plot the example data for S2 and S5 ===
    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket - 4]}{modifier}\*eepaa_fullcal*')
    data_dict = loadDictFromFile(inputFiles[0])
    Energy = data_dict['Energy'][0]
    Pitch = data_dict['Pitch_Angle'][0]

    wDisKeys = ['s2', 's3','s4','s5']
    wDis_indicies = [2, 3, 4, 5]

    for k, wDispersion_key in enumerate(wDisKeys):

        lowCut, highCut = np.abs(data_dict['Epoch'][0] - dispersionAttributes.keyDispersionDeltaT[wDis_indicies[k] - 1][0]).argmin(), np.abs(data_dict['Epoch'][0] - dispersionAttributes.keyDispersionDeltaT[wDis_indicies[k] - 1][1]).argmin()
        Epoch_dis = deepcopy(stl.dateTimetoTT2000(data_dict['Epoch'][0][lowCut:highCut + 1], inverse=False))
        Epoch_dis = (np.array(Epoch_dis) - Epoch_dis[0]) / 1E9  # converted data to TIME SINCE START OF DISPERSION (needed for proper plot fitting)

        # --- Plot everything ---
        if wDis_indicies[k] == 2:
            eepaa_dis = STEBfitResults[0][-1]
            isoIdx_row = 0
            isoIdx_col = 0
            ax = fig.add_subplot(gs00[isoIdx_row,isoIdx_col])
            ax.set_ylabel(r'$\mathbf{\alpha = 10}^{\circ}$' + '\nEnergy [eV]', fontsize=Label_FontSize, weight='bold')
            # ax.tick_params(axis='x', which='both', labelbottom=False)
        elif wDis_indicies[k] == 3:
            eepaa_dis = STEBfitResults[1][-1]
            isoIdx_row = 0
            isoIdx_col = 1
            ax = fig.add_subplot(gs00[isoIdx_row, isoIdx_col])
            ax.tick_params(axis='y', which='both', labelleft=False)
            # ax.tick_params(axis='x', which='both', labelbottom=False)
        elif wDis_indicies[k] == 4:
            eepaa_dis = STEBfitResults[2][-1]
            isoIdx_row = 1
            isoIdx_col = 0
            ax = fig.add_subplot(gs00[isoIdx_row, isoIdx_col])
            ax.set_ylabel(r'$\mathbf{\alpha = 10}^{\circ}$' + '\nEnergy [eV]', fontsize=Label_FontSize, weight='bold')
            ax.set_xlabel('Seconds [s]', fontsize=Label_FontSize, weight='bold')
        elif wDis_indicies[k] == 5:
            eepaa_dis = STEBfitResults[-1][-1]
            isoIdx_row = 1
            isoIdx_col = 1
            ax = fig.add_subplot(gs00[isoIdx_row, isoIdx_col])
            ax.tick_params(axis='y', which='both', labelleft=False)
            ax.set_xlabel('Seconds [s]', fontsize=Label_FontSize, weight='bold')

        cmap = ax.pcolormesh(Epoch_dis, Energy, eepaa_dis, cmap=mycmap, vmin=cbar_low, vmax=cbar_high,norm=None)
        ax.set_yscale('log')
        props = dict(boxstyle='round', facecolor='white', alpha=1, lw=4)
        ax.text(0.91*Epoch_dis[-1], 600, wDispersion_key.upper(), fontsize=Text_Fontsize, weight='bold', color='black', bbox=props, ha='center')
        ax.tick_params(axis='y', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
        ax.tick_params(axis='y', which='minor', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length * 0.6)
        ax.tick_params(axis='x', which='major', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)
        ax.tick_params(axis='x', which='minor', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length * 0.6)
        ax.set_ylim(27, 1000)



    # turn off the intentially blank plot
    # ax = fig.add_subplot(gsNull[0])
    # ax.axis('off')

    # === Plot fitted TOF data ===

    for rowIdx in range(2):
        for colIdx in range(2):
            ax = fig.add_subplot(gs01[rowIdx, colIdx])
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
            ax.text(0.25E-4, 0.6, f'S{wDis}', fontsize=Title_FontSize, weight='bold', color='black', bbox=props, ha='center',va='center')

            # Raw Data
            ax.scatter(deltaVs, deltaTs)
            ax.tick_params(axis='both', which='both', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length)

            if colIdx == 0:
                ax.set_ylabel('Delay time [s]',fontsize=Label_FontSize,weight='bold')

            if colIdx == 1:
                ax.tick_params(axis='y', which='both', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length, labelleft=False)

            if rowIdx == 0:
                ax.tick_params(axis='x', which='both', labelsize=Tick_FontSize, width=Tick_Width, length=Tick_Length,labelbottom=False)

            if rowIdx == 1:
                ax.set_xlabel('$\mathbf{\Delta (v^{-1})}$ [s/km]',fontsize=Label_FontSize,weight='bold')

            ax.plot(x_s, fitData, color="red", label=f'r-corr. (first order) = {round(r_corr_linear, 3)}')
            ax.plot(x_s, fitData_poly, color='black', label=f'r-corr. (second order) = {round(r_corr_poly,3)}')
            ax.set_xlim(0, 2.8E-4)
            ax.set_ylim(-0.15, 0.7)
            ax.legend(loc='lower right', fontsize=Legend_fontSize)
            ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
            ax.xaxis.offsetText.set_fontsize(24)
            counter += 1


    ##################
    # --- Colorbar ---
    ##################
    cax = fig.add_axes([0.87, 0.555, 0.025, 0.42])
    cbar = plt.colorbar(cmap, cax=cax)
    cbar.ax.minorticks_on()
    cbar.ax.tick_params(labelsize=Tick_FontSize, length=Tick_Length)

    cbar.ax.get_yaxis().labelpad = 17
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_weight("bold")
        l.set_fontsize(cbar_fontSize)

    cbar.set_label('Counts',rotation=90,fontsize=cbar_fontSize+10,weight='bold')

    plt.subplots_adjust(left=0.1, bottom=0.06, right=0.85, top=0.975, wspace=0.05, hspace=0.1)
    plt.savefig(rf'C:\Users\cfelt\Desktop\Research\Feltman2024_ACESII_Alfven_Observations\PLOTS\Plot7\Plot7_TOFfunc_base.png', dpi=dpi)


if generateLatexTable:
    from decimal import Decimal

    # Calculate the z_src
    m = stl.m_e
    q = stl.q0
    Emin = [28.22,28.22,28.22,28.22]
    Emax = [245.74,245.74,725.16,987.89]
    zsrc_linear = [round(val[3][0]) for val in STEBfitResults]
    a0 = [val[6][0] for val in STEBfitResults]
    a1 = [val[6][1]for val in STEBfitResults]
    a2 = [val[6][2]for val in STEBfitResults]
    r_corr_linear = [val[5] for val in STEBfitResults]
    r_corr_poly = [val[9] for val in STEBfitResults]
    def TOF_polyfit(x,Emax,a1,a2):
        return a2*(1E3 / np.cos(np.radians(10))) * np.sqrt(m/(2*q)) * ( x**(-0.5)  -Emax**(-0.5)) + a1

    z_src_poly = [ f'{ round(TOF_polyfit(Emin[i],Emax[i],a1[i],a2[i]) )} to { round(TOF_polyfit(Emax[i],Emax[i],a1[i],a2[i]) )}' for i in range(len(STEBfitResults))]


    from texttable import Texttable
    import latextable

    a1_temp = np.array(["%.2E" % Decimal(f"{val[6][1]}") for val in STEBfitResults])
    a2_temp = np.array(["%.2E" % Decimal(f"{val[6][2]}") for val in STEBfitResults])
    a1 = [rf'{val[0:4]} $\times 10^{val[-1]}$' for val in a1_temp]
    a2 = [rf'{val[0:5]} $\times 10^{val[-1]}$' for val in a2_temp]
    print(a1)
    print(a2)
    table_1 = Texttable()
    table_1.set_cols_align(["c", "c", "c","c","c","c"])
    table_1.add_rows([
        ['Parameters','S2','S3','S4','S5','unit'],
        ['E$_{min}$'] + Emin + ['eV'],
        ['E$_{max}$'] + Emax+ ['eV'],
        ['z$_{src}$ (linear)'] + zsrc_linear+['km'],
        ['r-corr (linear)'] +r_corr_linear+ ['-'],
        ['$a_{0}$'] + a0 + ['s'],
        ['$a_{1}$'] +a1 +['km'],
        ['$a_{2}$'] +a2 +['km$^{2}$/s'],
        ['r-corr (second order)'] + r_corr_poly+ ['-'],
        ['$z_{src}$ (second order)'] + z_src_poly + ['km'],
    ]
    )
    print(table_1.draw())
    print(latextable.draw_latex(table_1))




