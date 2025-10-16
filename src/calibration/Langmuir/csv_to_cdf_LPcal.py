# --- csv_to_cdf_LPcal.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: convert the LPcalibration files that were performed at iowa into .cdf files
# with all the data into one file.
import matplotlib.pyplot as plt
import numpy as np

# NOTE: The instrument gain was set too high, so for both instruments, we must do a semi-log plot of
# ln(current) vs ADC voltage and fit the linear regions as best we can. The swept probe will have
# two linear fits since the two ampliers at the end probably had different responses. Diode saturation
# will be considered to be 1nA. The "saturation" points in the cal curve shouldn't even be there unless we are at the extremes of the ADC +/- 4095
# Recall that our calibration setup ONLY had current coming from the power supply, so there is not natural
# atmospheric saturation occuring, thus the ADC should've seen current applied grow linearly in Ln(current), but it didn't
# This is because the instrument gain was set too high and we reached a point were the instrument could no
# long measure current accurately, thus we fit semi-log plots and find the linear regime and only fit those regions.


# --- --- --- --- ---
from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

justPrintFileNames = False
wRocket = 4 # 4 or 5 for high and low, respectively
inputPath_modifier = "calibration\LP_calibration\sweptCals"
outputPath_modifier = 'calibration\LP_calibration'
modifier = ''
# information to read the .csv files
wRow_names = 0 ; wRow_startData = 1  # row index corresponding to the name of the variable # which row does the data begin
getVars = [1, 3, 5, 7] # column index corresponding to the desired variables to collect from each file. These should be 2 - DAC, 4 - P1, 6 - ni_stepped, 8- ne_stepped

### ANALYSIS TOGGLES ###

useNanoAmps = False # all currents in nA, else it's Amps
unitConv = 1E9

plotIndividualCurveData = False # plot each individual curve
plotCalCurvesOverlay = False # plot the overlay of all the individual curves

# Size of analog window when determining the error in current derived from V_applied/R_resistor
analogWindowSize = 39 # 39 is the smallest window that has >=2 datapoints in each window bin
plotErrorsInCurrent = False

# plot the final set of data for ni and ne which will be fitted for calibration
plotFlattenedANDCleanedOverlay = False

# Average the calibration curves to get a smoothed result
AverageCalCurves = True # The point of averaging is to remove the assumed 60Hz noise that was present on the output of the LP input point. We then fit the averaged dataset with a linear line
plotAverageCalCurves = False
AveragingWindow = 100
polyOrder = 2

# plot the final calibration curve with the fits
n_e_fitRange = [0, 1820] # ADC range to fit calibration curve
n_i_fitRange = [0, 4250] # ADC range to fit calibration curve
plotCalibrationCurveFitted = True

# output the data
outputData = False



# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from csv import reader
from spacepy import coordinates as coord
import warnings
warnings.filterwarnings('ignore')
coord.DEFAULTS.set_values(use_irbem=False, itol=5)  # maximum separation, in seconds, for which the coordinate transformations will not be recalculated. To force all transformations to use an exact transform for the time, set ``itol`` to zero.
print(color.BOLD + color.CYAN + 'csv_to_cdf_LPcal.py' + color.END + color.END)

def csv_to_cdf_LPcal(wRocket, rocketFolderPath, justPrintFileNames):

    if wRocket == 4:
        wflyer = 0
    elif wRocket == 5:
        wflyer = 1

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, missionDict, data_dict_temp = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'science'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}{modifier}\*.csv')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}{modifier}\\', '') for ifile in inputFiles]

    fileoutName = 'ACESII_LP_Cals_swept_Iowa.cdf'

    if justPrintFileNames:
        for i, file in enumerate(input_names):
            print('[{:.0f}] {:20s}{:5.1f} MB'.format(i, input_names[i], round(getsize(inputFiles[i]) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data' + color.END)

        # --- Create Data Dict ---
        data_dict = {}

        #############################
        # --- STORE THE .CSV DATA ---
        #############################
        Rnames = []
        calCurveData_knownI = []
        calCurveData_ni = []
        calCurveData_ne = []
        resistanceVals = []
        for i, inputfile in enumerate(inputFiles):

            # get the filename identifier
            filenameOnly = inputfile.replace(f'{rocketFolderPath}{inputPath_modifier}{modifier}\\', '').replace('.csv','').replace('LP_','')
            Rnames.append(filenameOnly)
            prgMsg(f'Converting {filenameOnly}')

            # determine the resistance value of the file
            resitanceNum = filenameOnly[:filenameOnly.index("p")]
            resitanceDeci =filenameOnly[filenameOnly.index("p") + 1 : filenameOnly.index("k")]
            resistanceVal = (float(resitanceNum) + float(resitanceDeci)/1000)*1000
            resistanceVals.append(resistanceVal)

            # --- collect the csv data ---
            with open(inputfile) as csvFile:
                csvAllData = [row for row in reader(csvFile)]

            csvData = np.array(csvAllData[wRow_startData:], dtype='float64').transpose()

            # --- get v_stepped ---
            v_stepped = np.array(csvData[getVars[0]][:], dtype='float64')
            R = 4.997/(1.21 + 4.997)
            v_stepped = np.array([4*(R*vStep + (-5)*(1 - R)) for vStep in v_stepped]) # apply function to convert DAC output to v_stepped
            v_stepped_current = np.array(v_stepped / resistanceVal)


            # --- determine what the ADC would've reported. The ADC was a 12-bit, 5V converter ---
            ni_stepped = np.array(csvData[getVars[2]][:], dtype='float64')
            ne_stepped = np.array(csvData[getVars[3]][:], dtype='float64')

            # firstly, there was a diode that killed ANY voltage that was <0, so:
            ni_stepped = np.array([val if val > 0 else 0 for val in ni_stepped])
            ne_stepped = np.array([val if val > 0 else 0 for val in ne_stepped])

            # convert the measured ni/ne into the ideal ADC values
            ni_stepped_digital = np.array([(4096/5) * (ni_stepped[i]) for i in range(len(ni_stepped))])
            ne_stepped_digital = np.array([(4096/5) * (ne_stepped[i]) for i in range(len(ne_stepped))])

            # use the individual circuit data curves
            calCurveData_ni.append(ni_stepped_digital)
            calCurveData_ne.append(ne_stepped_digital)


            if plotIndividualCurveData:
                if i in range(len(inputFiles)):
                    fig, ax = plt.subplots(2)
                    ax[0].plot([j for j in range(len(ni_stepped_digital))], ni_stepped_digital, [j for j in range(len(ne_stepped_digital))], ne_stepped_digital)
                    ax[1].plot([j for j in range(len(v_stepped))],v_stepped)
                    plt.suptitle(filenameOnly)
                    ax[0].set_xlabel('Datapoint No.')
                    ax[0].set_ylabel('ADC Value')
                    ax[1].set_xlabel('Datapoint No.')
                    ax[1].set_ylabel('Probe Voltage [V]')
                    plt.show()

            # --- store data ---
            calCurveData_knownI.append(v_stepped_current)

            Done(start_time)

        # convert all the numpy array
        calCurveData_ni, calCurveData_ne, calCurveData_knownI = map(np.array, (calCurveData_ni, calCurveData_ne, calCurveData_knownI))

        # apply nanoAmp condition
        if useNanoAmps:
            currentUnits = 'nA'
            calCurveData_knownI = np.array([calCurveData_knownI[i]*unitConv for i in range(len(calCurveData_knownI))])
        else:
            currentUnits = 'A'

        # store everything in a dictonary
        cali_dict = {
            'ni': calCurveData_ni,
            'ni_knownI_avg' : [],
            'ni_knownI': calCurveData_knownI,
            'ni_errors':[],
            'ni_error_bins':[],
            'ne': calCurveData_ne,
            'ne_knownI_avg': [],
            'ne_knownI': calCurveData_knownI,
            'ne_errors': [],
            'ne_error_bins': []
        }

        # --- --- --- --- --- --- --- --- --- --- -
        # --- FILTER THE CALIBRATION CURVE DATA ---
        # --- --- --- --- --- --- --- --- --- --- -


        if plotCalCurvesOverlay:
            legend = []
            fig, ax = plt.subplots(2)
            currents = [cali_dict['ne'], cali_dict['ni']]
            labels = ['ne','ni']
            for j in range(2):
                for i in range(len(calCurveData_knownI)): #loop through all the inidividual curves
                    stepped_current = np.array(calCurveData_knownI[i])
                    sweep = currents[j][i]
                    ax[j].scatter(sweep, stepped_current)
                    ax[j].set_xlabel('ADC value it wouldve Recorded [Analog]')
                    ax[j].set_ylabel(f'Current derived from V_applied/R [{currentUnits}]')
                    # ax[j].set_yscale('symlog',linthresh=0.1,base=np.exp(1))
                    ax[j].set_title(f'{labels[j]}')
                    legend.append(Rnames[i])

                plt.legend(legend)
            plt.show()

        # Flatten all the data and convert it to nano-amps
        cali_dict['ne_knownI'] = np.array([item for sublist in deepcopy(cali_dict['ne_knownI']) for item in sublist])
        cali_dict['ni_knownI'] = np.array([item for sublist in deepcopy(cali_dict['ni_knownI']) for item in sublist])
        cali_dict['ne'] = np.array([item for sublist in deepcopy(cali_dict['ne']) for item in sublist])
        cali_dict['ni'] = np.array([item for sublist in deepcopy(cali_dict['ni']) for item in sublist])

        # --- --- --- --- --- --- --- --- --- -
        # --- REMOVE THE ==0 ADC ==0 points ---
        # --- --- --- --- --- --- --- --- --- -
        # we don't care about when n_i/n_e == 0 since the v_applied was pos/neg, respectively, lets remove those points
        # ni
        badIndicies = np.where(cali_dict['ni'] == 0)[0]
        cali_dict['ni_knownI'] = np.delete(cali_dict['ni_knownI'], badIndicies)
        cali_dict['ni'] = np.delete(cali_dict['ni'], badIndicies)
        # ne
        badIndicies = np.where(cali_dict['ne'] == 0)[0]
        cali_dict['ne_knownI'] = np.delete(cali_dict['ne_knownI'], badIndicies)
        cali_dict['ne'] = np.delete(cali_dict['ne'], badIndicies)

        # --- Sort the data ---
        # assign the data and also sort it based on the analog values
        cali_dict['ni'], cali_dict['ni_knownI'] = list(zip(*sorted(zip(cali_dict['ni'], cali_dict['ni_knownI']))))
        cali_dict['ne'], cali_dict['ne_knownI'] = list(zip(*sorted(zip(cali_dict['ne'], cali_dict['ne_knownI']))))

        # convert everything to numpy
        for key, val in cali_dict.items():
            if len(val) > 0:
                cali_dict[key] = np.array(cali_dict[key])

        # --- CORRECT INPUT VOLTAGE NOISE PROBLEM ---
        # NOTE: The input voltage to the circuit was noisey sometimes. This means that around 0V, the current swung positive and
        # negative, which means n_i can sometimes see positive current and n_e can see negative current, which shouldn't be possible
        # normally. Here I just remove these points for fitting purposes.
        labels = ['ne', 'ni']
        for i in range(2):
            if labels[i] == 'ne':
                badIndicies = np.where(cali_dict[f'{labels[i]}_knownI'] < 0)[0]
            else:
                badIndicies = np.where(cali_dict[f'{labels[i]}_knownI'] > 0)[0]
            cali_dict[f'{labels[i]}'] = np.delete(deepcopy(cali_dict[f'{labels[i]}']), badIndicies)
            cali_dict[f'{labels[i]}_knownI'] = np.delete(deepcopy(cali_dict[f'{labels[i]}_knownI']),badIndicies)

        if plotFlattenedANDCleanedOverlay:
            fig, ax = plt.subplots(2)
            labels = ['ne', 'ni']
            for j in range(2):
                ax[j].scatter(cali_dict[f'{labels[j]}'], cali_dict[f'{labels[j]}_knownI'])
                ax[j].set_xlabel(f'{labels[j]} ADC [Analog]')
                ax[j].set_ylabel(f'Known I from Cal. Resistor [{currentUnits}]')
                ax[j].set_yscale('symlog',linthresh=0.001,base=np.exp(1))
                ax[j].set_title(f'{labels[j]}')
                # ax[j].set_xscale('log')
            plt.show()

        # --- COLLECT LP ERROR DATA ---
        # Both ne/ni have error in their determination of the known current, gather that info now

        prgMsg('Collecting Error Data')
        analogBins_ne = np.array([[num, num+analogWindowSize] for num in range(int(cali_dict['ne'].min()) - 1, int(cali_dict['ne'].max())+1, analogWindowSize)])
        analogBins_ni = np.array([[num, num+analogWindowSize] for num in range(int(cali_dict['ni'].min()) - 1, int(cali_dict['ni'].max())+1, analogWindowSize)])
        stdDevErrors_ne = []
        stdDevErrors_ni = []

        for bin in analogBins_ne:
            errors = [ np.log(cali_dict['ne_knownI'][i]) for i in range(len(cali_dict['ne'])) if bin[0] <= cali_dict['ne'][i] < bin[1]]
            stdDevErrors_ne.append(np.std(errors))

        for bin in analogBins_ni:
            errors = [np.log(-1*cali_dict['ni_knownI'][i]) for i in range(len(cali_dict['ni'])) if bin[0] <= cali_dict['ni'][i] < bin[1]]
            stdDevErrors_ni.append(np.std(errors))

        # store all the information
        cali_dict['ne_error_bins'] = np.array(analogBins_ne)
        cali_dict['ni_error_bins'] = np.array(analogBins_ni)
        cali_dict['ne_errors'] = np.array(stdDevErrors_ne)
        cali_dict['ni_errors'] = np.array(stdDevErrors_ni)

        if plotErrorsInCurrent:
            xData_ne = [bin[0] for bin in analogBins_ne]
            yData_ne = stdDevErrors_ne
            xData_ni = [bin[0] for bin in analogBins_ni]
            yData_ni = stdDevErrors_ni

            fig, ax = plt.subplots(2)
            fig.suptitle(f'$\Delta I$ vs Analog Value\n Analog Bin Size: {analogWindowSize}')
            ax[0].scatter(xData_ne, yData_ne)
            ax[0].set_xlabel('Lower Bound of Analog Bin [Analog]')
            ax[0].set_ylabel(f'ne Ln(I_known) Std Dev [Ln({currentUnits})]')
            ax[1].scatter(xData_ni, yData_ni)
            ax[1].set_xlabel('Lower Bound of Analog Bin [Analog]')
            ax[1].set_ylabel(f'ne Ln(I_known) Std Dev [Ln({currentUnits})]')
            plt.show()
        Done(start_time)


        # --- AVERAGE LP CAL CURVES ---
        # The point of averaging is to remove the assumed 60Hz noise that was present on the output of the LP input point. We then
        # fit the averaged dataset with a linear line

        if AverageCalCurves:
            prgMsg('Averaging calibration Curves')

            cali_dict['ni_knownI'] = np.array([np.log(-1*cur) for cur in cali_dict['ni_knownI']])
            cali_dict['ne_knownI'] = np.array([np.log(cur) for cur in cali_dict['ne_knownI']])

            from scipy.signal import savgol_filter
            cali_dict['ni_knownI_avg'] = savgol_filter(cali_dict['ni_knownI'], window_length=AveragingWindow, polyorder=polyOrder)
            cali_dict['ne_knownI_avg'] = savgol_filter(cali_dict['ne_knownI'], window_length=AveragingWindow, polyorder=polyOrder)

            if plotAverageCalCurves:

                fig, ax = plt.subplots(2)
                labels = ['ne', 'ni']
                fig.suptitle('SG filtering Cal data')
                for j in range(2):
                    ax[j].scatter(cali_dict[f'{labels[j]}'], cali_dict[f'{labels[j]}_knownI'])
                    ax[j].plot(cali_dict[f'{labels[j]}'], cali_dict[f'{labels[j]}_knownI_avg'],color='red')
                    ax[j].set_xlabel(f'{labels[j]} ADC [Analog]')
                    ax[j].set_ylabel(f'Ln(I from Cal) [{currentUnits}]')
                    # ax[j].set_yscale('symlog', linthresh=0.001, base=np.exp(1))
                    ax[j].set_title(f'{labels[j]}')
                    # ax[j].set_xscale('log')
                plt.show()

            Done(start_time)
        else:
            R1 = 100
            R2 = 1000
            testDat =  cali_dict['ni_knownI'][R1:R2]
            cali_dict['ni_knownI'] = np.array([np.log(-1 * cur) for cur in cali_dict['ni_knownI']])
            testDat2 = cali_dict['ni_knownI'][R1:R2]
            cali_dict['ne_knownI'] = np.array([np.log(cur) for cur in cali_dict['ne_knownI']])

            for i in range(len(testDat)):
                print(testDat[i],testDat2[i])

        # --- --- --- --- --- --- -
        # --- LINEARLY FIT DATA ---
        # --- --- --- --- --- --- -
        from scipy.optimize import curve_fit
        def linear(x, A, B):
            y = A * x + B
            return y

        # only get the data to fit within my specified range

        dataModifier = '_avg' if AverageCalCurves else ''
        xData_ne = []
        yData_ne = []
        xData_ni = []
        yData_ni = []

        for i in range(len(cali_dict['ne'])):
            if n_e_fitRange[0] <= cali_dict['ne'][i] <= n_e_fitRange[1]:
                xData_ne.append(cali_dict['ne'][i])
                yData_ne.append(cali_dict[f'ne_knownI{dataModifier}'][i])

        for i in range(len(cali_dict['ni'])):
            if n_i_fitRange[0] <= cali_dict['ni'][i] <= n_i_fitRange[1]:
                xData_ni.append(cali_dict['ni'][i])
                yData_ni.append(cali_dict[f'ni_knownI{dataModifier}'][i])

        params_ne, cov = curve_fit(linear, xData_ne, yData_ne)
        params_ni, cov = curve_fit(linear, xData_ni, yData_ni)

        if plotCalibrationCurveFitted:

            fig, ax = plt.subplots(2)

            # ne
            ax[0].set_title('ne')
            ax[0].scatter(cali_dict['ne'], cali_dict['ne_knownI'])
            xData = np.linspace(min(cali_dict['ne']), max(cali_dict['ne']), 1000)
            yData = [linear(x, *params_ne) for x in xData]
            ax[0].plot(xData, yData, color='purple')
            ax[0].set_xlabel(r'ADC Value')
            ax[0].set_ylabel(r'Ln($I_{cal}$)')

            #ni
            ax[1].set_title('ni')
            ax[1].scatter(cali_dict['ni'], cali_dict['ni_knownI'])
            xData = np.linspace(min(cali_dict['ni']), max(cali_dict['ni']), 1000)
            yData = [linear(x, *params_ni) for x in xData]
            ax[1].plot(xData, yData, color='purple')
            ax[1].set_xlabel(r'ADC Value')
            ax[1].set_ylabel(r'Ln($I_{cal}$)')

            if AverageCalCurves:
                ax[0].plot(cali_dict['ne'], cali_dict['ne_knownI_avg'], color='red')
                ax[1].plot(cali_dict['ni'], cali_dict['ni_knownI_avg'], color='red')

            plt.legend([f'Func n_e: A*x + B\n'
                        f'A: {params_ne[0]}\n'
                        f'B: {params_ne[1]}\n',
                        f'Func n_i: A*x + B\n'
                        f'A: {params_ni[0]}\n'
                        f'B: {params_ni[1]}\n'
                        ])
            plt.suptitle('V/R resistor Current vs Analog Circuit Response\n'
                         f'ACESII - {rocketID}')
            plt.show()


        ###################################
        ### CALULCATE ERRORS IN CURRENT ###
        ###################################

        # The analog bins for the error in the current need to be converted into current themselves
        # in order to be useful

        analogBins_current_ne = []  # bins for the error in the current
        analogBins_current_ni = []

        for i in range(len(cali_dict['ne_error_bins'])):
                analogBins_current_ne.append(
                    [linear(cali_dict['ne_error_bins'][i][0], *params_ne),
                     linear(cali_dict['ne_error_bins'][i][1], *params_ne)])

        for i in range(len(cali_dict['ni_error_bins'])):
                analogBins_current_ni.append(
                    [linear(cali_dict['ni_error_bins'][i][0], *params_ni),
                     linear(cali_dict['ni_error_bins'][i][1], *params_ni)])

        cali_dict['ne_error_bins'], cali_dict['ni_error_bins'] = np.array(analogBins_current_ne), np.array(analogBins_current_ni)

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        if outputData:
            prgMsg('Creating output file')

            # store the fit parameters
            fit_params = np.array([params_ne,params_ni])
            data_dict = {**data_dict,
                         **{'fit_params': [fit_params, {'LABLAXIS': 'fit_params',
                                                              'DEPEND_0': None,
                                                              'DEPEND_1': None,
                                                              'DEPEND_2': None,
                                                              'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                              'UNITS': 'Units',
                                                              'VALIDMIN': fit_params.min(),
                                                              'VALIDMAX': fit_params.max(),
                                                              'VAR_TYPE': 'data',
                                                              'SCALETYP': 'linear'}]}}

            # store the error in the current and the respective bins
            stdDevErrors = np.array(cali_dict['ne_errors'])
            data_dict = {**data_dict,
                         **{'ne_errorInCurrent': [stdDevErrors, {'LABLAXIS': 'ne_errorInCaldCurrent',
                                                                        'DEPEND_0': None,
                                                                        'DEPEND_1': None,
                                                                        'DEPEND_2': None,
                                                                        'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                                        'UNITS': f'{currentUnits}',
                                                                        'VALIDMIN': stdDevErrors.min(),
                                                                        'VALIDMAX': stdDevErrors.max(),
                                                                        'VAR_TYPE': 'data',
                                                                        'SCALETYP': 'linear'}]}}
            stdDevErrors = np.array(cali_dict['ni_errors'])
            data_dict = {**data_dict,
                         **{'ni_errorInCurrent': [stdDevErrors, {'LABLAXIS': 'ni_errorInCaldCurrent',
                                                                 'DEPEND_0': None,
                                                                 'DEPEND_1': None,
                                                                 'DEPEND_2': None,
                                                                 'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                 'FORMAT': 'E12.2',
                                                                 'UNITS': f'{currentUnits}',
                                                                 'VALIDMIN': stdDevErrors.min(),
                                                                 'VALIDMAX': stdDevErrors.max(),
                                                                 'VAR_TYPE': 'data',
                                                                 'SCALETYP': 'linear'}]}}

            analogBins = np.array(cali_dict['ne_error_bins'])
            data_dict = {**data_dict,
                         **{'ne_errorBins': [analogBins, {'LABLAXIS':'ne_errorBins',
                                                              'DEPEND_0': None,
                                                              'DEPEND_1': None,
                                                              'DEPEND_2': None,
                                                              'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                              'UNITS': 'analogValue',
                                                              'VALIDMIN': analogBins.min(),
                                                              'VALIDMAX': analogBins.max(),
                                                              'VAR_TYPE': 'support_data',
                                                              'SCALETYP': 'linear'}]}}
            analogBins = np.array(cali_dict['ni_error_bins'])
            data_dict = {**data_dict,
                         **{'ni_errorBins': [analogBins, {'LABLAXIS': 'ni_errorBins',
                                                       'DEPEND_0': None,
                                                       'DEPEND_1': None,
                                                       'DEPEND_2': None,
                                                       'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                       'UNITS': 'analogValue',
                                                       'VALIDMIN': analogBins.min(),
                                                       'VALIDMAX': analogBins.max(),
                                                       'VAR_TYPE': 'support_data',
                                                       'SCALETYP': 'linear'}]}}

            # store the regions where the different fits apply
            fitRegions = np.array([n_e_fitRange, n_i_fitRange])
            data_dict = {**data_dict,
                         **{'fitRegions': [fitRegions, {'LABLAXIS': 'fitRegions',
                                                                       'DEPEND_0': None,
                                                                       'DEPEND_1': None,
                                                                       'DEPEND_2': None,
                                                                       'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                       'FORMAT': 'E12.2',
                                                                       'UNITS': 'ADC',
                                                                       'VALIDMIN': 0,
                                                                       'VALIDMAX': 1,
                                                                       'VAR_TYPE': 'support_data',
                                                                       'SCALETYP': 'linear'}]}}

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'
            globalAttrsMod['Descriptor'] = "Langmuir_Probe"
            outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod, 'LangmuirProbe')

            Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 4:  # ACES II Integration High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Integration Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\*.csv')) == 0:
    print(color.RED + 'There are no .csv files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        csv_to_cdf_LPcal(wRocket, rocketFolderPath, justPrintFileNames)
    else:
        csv_to_cdf_LPcal(wRocket, rocketFolderPath, justPrintFileNames)


