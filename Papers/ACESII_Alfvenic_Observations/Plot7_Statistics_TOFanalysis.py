# --- crossCorrelationTOFAnalysis.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Plot the TOF results. DATA IS GOTTEN FROM  Sciece>AlfvenSignatureAnalysis>Particles>CrossCorrelationTOFAnaylsis

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---

wRocket = 4

# --- IMPORTS ---
from ACESII_code.Science.AlfvenSingatureAnalysis.Particles.dispersionAttributes import dispersionAttributes
from ACESII_code.class_var_func import Re,u0,ep0,q0,m_e
from matplotlib.dates import date2num, num2date
from matplotlib.lines import Line2D


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
figure_width = 20
figure_height = 10
textFontSize = 25
labelSize = 27
labelPadding = 13
tickLabelSize = 25
tickWidth = 3
tickLength = 16
titleSize = 10
legendSize = 24
errorBarMarkerSize = 10
legendMarkerSize = 18
lineWidth = 10
scatterSize = 650
dpi = 400

showErrorBars = False
remove_STEB18_and_STEB12 = True # steb18 is really far away from the rest and it makes the plot hard to understand. STEB12 is barely perceptible
dT_colors = ['tab:blue', 'tab:green', 'tab:red']

# Numbered Labels Alignment
verticalAligns = 0.08*np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
horizontalAligns = 0.05*np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])


def Plot5_TOFanalysis(rocketFolderPath):

    # load in the data
    inputFile_TOFresults = glob('C:\Data\ACESII\science\TOFanalysis\high\ACESII_36359_crossCorrelationAnalysis.cdf')[0]  # get the deltaB data
    inputFile_EEPAA = glob(f'C:\Data\ACESII\L1\high\*eepaa_fullcal*')[0]
    inputFile_attitude = r'C:\Data\ACESII\attitude\high\ACESII_36359_Attitude_Solution.cdf'

    # Load the data ESA... again
    data_dict_TOF = loadDictFromFile(inputFile_TOFresults)
    data_dict_eepaa = loadDictFromFile(inputFile_EEPAA)
    data_dict_attitude = loadDictFromFile(inputFile_attitude)

    # define some variables
    Energy = np.array(data_dict_eepaa['Energy'][0])
    Pitch = data_dict_eepaa['Pitch_Angle'][0]

    ###########################
    # --- Organize the data ---
    ###########################

    # FORMAT of STEB fit results:
    # [[STEB Number,Observation Time, Observation Altitude, Z_acc, Z_acc_error_avg, correlationR],
    #  [...],
    #  ]
    STEBtimes = data_dict_TOF['whenSTEBoccured_time'][0]
    STEB_ILats = [data_dict_attitude['ILat'][0][np.abs(data_dict_attitude['Epoch'][0] - tme).argmin()] for tme in STEBtimes]
    STEB_Zacc = data_dict_TOF['whenSTEBoccured_Alt'][0]/(1000*Re)+data_dict_TOF['Z_acc'][0]
    STEB_Zacc_avg = sum(STEB_Zacc) / len(STEB_Zacc)
    STEB_errorZacc_avg = data_dict_TOF['Z_acc_error'][0]
    STEB_deltaE = []
    STEB_engyBars = []

    for engyPair in dispersionAttributes.keyDispersionDeltaE:
        ELow = Energy[np.abs(Energy - engyPair[0]).argmin()]
        EHigh = Energy[np.abs(Energy - engyPair[1]).argmin()]
        STEB_deltaE.append(EHigh-ELow)
        STEB_engyBars.append([EHigh - 2*ELow, ELow])

    # Find deltaT
    STEB_deltaT = np.array([pycdf.lib.datetime_to_tt2000(tmePair[1]) - pycdf.lib.datetime_to_tt2000(tmePair[0]) for tmePair in dispersionAttributes.keyDispersionDeltaT]) / 1E9

    # Find the labels to use for each STEB
    STEB_text_labels = []
    colorChoice = ['tab:orange', 'tab:purple', 'tab:cyan']
    for type in dispersionAttributes.keyDispersionType:
        if 'isolat' in type.lower():
            STEB_text_labels.append(['I', colorChoice[0]])
        elif 'aurora' in type.lower():
            STEB_text_labels.append(['A', colorChoice[1]])
        elif 'edge' in type.lower():
            STEB_text_labels.append(['E', colorChoice[2]])

    STEB_text_labels = np.array(STEB_text_labels)

    # determine the colors
    dT_color_labels = []
    for delT in STEB_deltaT:
        if delT <= 0.5:
            dT_color_labels.append(dT_colors[0])
        elif 0.5 < delT <= 0.9:
            dT_color_labels.append(dT_colors[1])
        elif delT > 0.95:
            dT_color_labels.append(dT_colors[2])

    # remove STEB 18 from the data
    if remove_STEB18_and_STEB12:
        badIndicies = [12, 17]
        STEBtimes = np.delete(STEBtimes, badIndicies)
        STEB_ILats = np.delete(STEB_ILats, badIndicies)
        STEB_Zacc = np.delete(STEB_Zacc, badIndicies)
        STEB_deltaE = np.delete(STEB_deltaE, badIndicies)
        STEB_engyBars = np.delete(STEB_engyBars, badIndicies, axis=0)
        STEB_text_labels = np.delete(STEB_text_labels, badIndicies, axis=0)
        STEB_deltaT = np.delete(STEB_deltaT, badIndicies)

    ##########################
    # --- PLOT THE RESULTS ---
    ##########################
    prgMsg('Making Plot5')
    fig, ax = plt.subplots(nrows=2)
    fig.set_figwidth(figure_width)
    fig.set_figheight(figure_height)

    #--- Z_ACC PLOT ---
    cmap = ax[0].scatter(x=STEBtimes, y=STEB_Zacc, s=scatterSize, c=STEB_text_labels[:, 1])
    ax[0].axhline(STEB_Zacc_avg, linestyle='--', color='black', linewidth=3, zorder=0)
    ax[0].text(dt.datetime(2022, 11, 20, 17, 24, 40), 1.05*STEB_Zacc_avg, f"Average: {round(STEB_Zacc_avg, 2)} $R_E$ ", fontsize=textFontSize, color='black', va='bottom')
    ax[0].set_ylabel('Source Altitude [$R_{E}$]', color='black', fontsize=labelSize, labelpad=labelPadding+8,weight='bold')
    ax[0].set_ylim(0, 1)

    # Scatter Labels
    for l in range(len(STEBtimes)):
        ax[0].text(x=STEBtimes[l], y=STEB_Zacc[l]+verticalAligns[l],
                   s=f'{l + 1}', va='center', ha='center',
                   weight='bold', fontsize=textFontSize-5, zorder=1)

    # Ticks
    # startPoint, endPoint = date2num(dt.datetime(2022, 11, 20, 17, 24, 18, 000)), date2num(dt.datetime(2022, 11, 20, 17, 25, 25))
    startPoint, endPoint = date2num(min(STEBtimes)), date2num(max(STEBtimes))
    EpochTicks_dt = num2date(np.linspace(startPoint, endPoint, 8))
    EpochTicks = np.array([ round((pycdf.lib.datetime_to_tt2000(tick) - pycdf.lib.datetime_to_tt2000(dt.datetime(2022,11,20,17,20,000,000)))/1E9,1) for tick in EpochTicks_dt])
    offSetAwareAttitudeEpoch = date2num(data_dict_attitude['Epoch'][0])
    iLatTicks = [round(data_dict_attitude['ILat'][0][np.abs(offSetAwareAttitudeEpoch - date2num(tme)).argmin()], 2) for tme in EpochTicks_dt]

    ax[0].set_xticks(ticks=EpochTicks_dt, labels=iLatTicks)
    ax[0].minorticks_on()
    ax[0].tick_params(axis='both', which='major', labelsize=tickLabelSize, width=tickWidth, length=tickLength)
    ax[0].tick_params(axis='both', which='minor', labelsize=tickLabelSize, width=tickWidth, length=tickLength/2)
    ax[0].grid(True)
    ax[0].set_xticks(ticks=EpochTicks_dt, labels=['' for thing in EpochTicks])

    # create Epoch axis
    axTime = ax[0].twiny()
    axTime.scatter(STEBtimes, STEB_Zacc, alpha=0,s=scatterSize+50,zorder=1)
    axTime.set_xlabel('Time from launch [s]', fontsize=labelSize, labelpad=labelPadding+5,weight='bold')

    # Ticks
    # axTime.set_xticks(ticks=EpochTicks, labels=[ep.strftime('%H:%M.%S') for ep in EpochTicks])
    axTime.set_xticks(ticks=EpochTicks_dt, labels=[str(tick) for tick in EpochTicks])
    axTime.tick_params(axis='x', which='major', labelsize=tickLabelSize, width=tickWidth, length=tickLength)
    axTime.tick_params(axis='x', which='minor', labelsize=tickLabelSize, width=tickWidth, length=tickLength/2)
    axTime.minorticks_on()

    # custom legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=colorChoice[0], label='Isolated', markersize=legendMarkerSize),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=colorChoice[1], label='Auroral', markersize=legendMarkerSize),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=colorChoice[2], label='Edge-Type', markersize=legendMarkerSize)
    ]

    ax[0].legend(handles=legend_elements, loc='upper right', prop={'size': legendSize})

    # --- DELTA E PLOT ---
    for k in range(len(STEBtimes)):
        ax[1].errorbar(x=STEBtimes[k], y=[STEB_deltaE[k]], yerr=np.array([[STEB_engyBars[k][0]],[STEB_engyBars[k][1]]]),
                       fmt='o',markersize=errorBarMarkerSize,color=dT_color_labels[k],linewidth=4, capsize=legendMarkerSize-8)

    ax[1].grid(True)
    ax[1].set_ylim(0,1000)
    ax[1].set_xticks(ticks=EpochTicks_dt, labels=iLatTicks)
    ax[1].minorticks_on()
    ax[1].tick_params(axis='both', which='major', labelsize=tickLabelSize, width=tickWidth, length=tickLength)
    ax[1].tick_params(axis='both', which='minor', labelsize=tickLabelSize, width=tickWidth, length=tickLength/2)
    ax[1].set_ylabel('STEB $\Delta E$ [eV]', fontsize=labelSize, labelpad=labelPadding-3,weight='bold')
    ax[1].set_xlabel('ILat [deg]', color='black', fontsize=labelSize, labelpad=labelPadding-3,weight='bold')
    legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=dT_colors[0], label='0.242s $\leq \Delta T \leq$ 0.5s', markersize=legendMarkerSize),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=dT_colors[1], label='0.5s < $\Delta T \leq$ 1s', markersize=legendMarkerSize),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=dT_colors[2], label='1s $\leq \Delta T \leq$ 1.405s', markersize=legendMarkerSize)
    ]
    ax[1].legend(handles=legend_elements, loc='upper right', prop={'size': legendSize})

    # for l in range(len(STEBtimes)):
    #     ax[1].text(x=STEBtimes[l], y=STEB_deltaE[l], s=f'{l + 1}', va=verticalAligns[l], ha=horizontalAligns[l], weight='bold', fontsize=textFontSize)

    # Show the figure
    plt.tight_layout()
    plt.savefig(rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot5\STEB_CorrelationAnalysis_Results.png', dpi=dpi)
    Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
Plot5_TOFanalysis(rocketFolderPath)
