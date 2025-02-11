#--- L1_to_L2_langmuir_fixed.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert the engineering Langmuir data to scientifically useful units. Also renames
# "Boom_Monitor_1/2" --> "Fixed_Boom_Monitor_1/2" etc

# it was discovered that ni_swept for the low AND high Flyer start with an extra value that should not be there.
# This is fixed in L0_to_L1.py

# It was discovered that the gain for the Swept LPs is probably too high. This means the instrument enters saturation
# very quickly. So both the real data and calibration data saturate too quickly


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

# --- --- --- --- ---
from src.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

#################
# --- TOGGLES ---
#################

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4

#################
# --- TOGGLES ---
#################
SECTION_fixedProbeCal = True
applyFixedCalCurve = False
plotFixedCalCurve = True

# --- DATA OUTPUT ---
outputData = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import scipy

#######################
# --- MAIN FUNCTION ---
#######################

def L1_to_L2_langmuir_fixed(wflyer, wFile, justPrintFileNames):

    # --- Get ACESII rocket Attributes ---
    rocket_folder_path = DataPaths.ACES_data_folder
    rocketID = ACESII.payload_IDs[wflyer]
    globalAttrsMod = deepcopy(ACESII.globalAttributes[wflyer])
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'Langmuir'

    # Set the paths for the file names
    L1Files = glob(f'{rocketFolderPath}L1\{fliers[wflyer]}\*_lp_*')
    LangmuirFiles = glob(f'{outputFolderPath}\{fliers[wflyer]}\*_langmuir_*')
    LangmuirSweptCalFiles = glob(f'{rocketFolderPath}\calibration\LP_calibration\{fliers[wflyer]}\*_345deg_*')
    L1_names = [ifile.replace(f'{rocketFolderPath}L1\{fliers[wflyer]}\\', '') for ifile in L1Files]
    L1_names_searchable = [ifile.replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace('l1_', '').replace('_v00', '') for ifile in L1_names]
    dataFile_name = L1_names_searchable[0].replace(f'{rocketFolderPath}L1\{fliers[wflyer]}\\', '').replace('lp_','').replace('ni_','').replace('ne_swept_','').replace('step_','').replace('ni_swept_','').replace('deltaNdivN_','')
    fileoutName = rf'ACESII_{rocketAttrs.rocketID[wRocket-4]}_l2_langmuir' + '.cdf'
    wInstr = [0, 'LangmuirProbe']

    if justPrintFileNames:
        for i, file in enumerate(L1Files):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, L1_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return

    print('\n')
    print(color.UNDERLINE + f'Creating langmuir data for {fliers[wRocket-4]} flyer' + color.END)

    # --- get the data from the tmCDF file ---
    prgMsg('Loading data from L1Files')

    # Collect the LP data except deltaNdivN into a data dict
    data_dict = {}
    for file in L1Files:
        if 'deltaNdivN' not in file:
            with pycdf.CDF(file) as L1DataFile:
                for key, val in L1DataFile.items():
                    if key not in data_dict:
                        data_dict = {**data_dict, **{key : [L1DataFile[key][...] , {key:val for key,val in L1DataFile[key].attrs.items()  }  ]  }  }

    # Collect the LP data except deltaNdivN into a data dict
    data_dict_cal = {}
    for file in LangmuirSweptCalFiles:
        if 'deltaNdivN' not in file:
            with pycdf.CDF(file) as LangmuirSweptCalFiles:
                for key, val in LangmuirSweptCalFiles.items():
                    if key not in data_dict_cal:
                        data_dict_cal = {**data_dict_cal, **{key: [LangmuirSweptCalFiles[key][...], {key: val for key, val in LangmuirSweptCalFiles[key].attrs.items()}]}}
    Done(start_time)

    #####################
    # --- Fixed Probe ---
    #####################
    # description: Loads in the fixed calibrated data (provided by scott, found in missionAttributes.py).

    prgMsg('Converting fixed probe to Voltage')
    def calFunction_fixed(x, A, B):
        y = A*x + B
        return y

    fixedCalResistances = rocketAttrs.LPFixed_calResistances[wRocket - 4]
    probeBias = rocketAttrs.LPFixedProbeBias[wRocket - 4]
    calibrationCurrents = []
    analog_vals = []

    # convert calibration data to current (in nanoamps)
    # NOTE: WE don't care about the "open case"
    for key, val in fixedCalResistances.items():
        if key != 'Open':
            analog_vals.append(val)

            if useNanoAmps:
                calibrationCurrents.append(unit_conversion * probeBias / key)
            else:
                calibrationCurrents.append(probeBias / key)

    analog_vals, calibrationCurrents = np.array(analog_vals), np.array(calibrationCurrents)

    # apply a natural log scale to the data in order to prepare it for fitting
    calibrationCurrents = np.array([np.log(-1*cur) for cur in calibrationCurrents])

    # Fit a linear line to the log'd data
    parameters, covariance = scipy.optimize.curve_fit(calFunction_fixed, analog_vals, calibrationCurrents, maxfev=10000)

    if plotFixedCalCurve:
        figWidth =10
        figHeight =10
        Label_Fontsize = 15
        Title_Fontsize = 15

        import matplotlib
        matplotlib.rc('figure', figsize=(5, 5))
        plt.rcParams["figure.figsize"] = (figHeight, figWidth)
        xDataFit = np.array([i for i in range(1, 4096)])
        yDataFit = [calFunction_fixed(val, *parameters) for val in xDataFit]
        plt.plot(xDataFit, yDataFit, color='red')
        plt.scatter(analog_vals, calibrationCurrents)
        plt.xlabel('ADC Value',fontsize=Label_Fontsize)
        plt.ylabel(r'Ln($I_{cal}$) [nA]',fontsize=Label_Fontsize)
        plt.suptitle(f'FIXED LP - {rocketAttrs.rocketID[wRocket-4]}\n'
                     'Calculated Calibration Current vs Analog Value',fontsize=Title_Fontsize)
        plt.legend(['ln(y) = mx + b\n'
                   f'm: {parameters[0]}\n'
                   f'b: {parameters[1]}'])
        plt.savefig(rf'C:\Users\cfelt\Desktop\rockets\ACES-II\Instruments\fixedLangmuirProbe\Calibration\ACESII_{fliers[wflyer]}_CalCurve.png')

    # Apply the calibration function curve
    index = np.abs(data_dict['Epoch_ni'][0] - dt.datetime(2022,11,20,17,25,44,500000)).argmin()

    if applyFixedCalCurve:
        caldCurrent = np.array([
            np.exp(calFunction_fixed(data_dict['ni'][0][i], parameters[0],parameters[1]))
            for i in range(len(data_dict['ni'][0]))])
        fixedLPunits = 'nA'
    else:
        caldCurrent = np.array(data_dict['ni'][0])
        fixedLPunits = 'ADC'

    # apply a quality assurance step:
    if applyFixedCalCurve:
        if wRocket == 4:
            for i in range(len(caldCurrent)):
                if np.abs(caldCurrent[i]) > 2000:
                    caldCurrent[i] = rocketAttrs.epoch_fillVal
        elif wRocket == 5:
            for i in range(len(caldCurrent)):
                if np.abs(caldCurrent[i]) > 2000:
                    caldCurrent[i] = rocketAttrs.epoch_fillVal

    # --- --- --- --- --- --- ----
    # --- FIXED ERROR ANALYSIS ---
    # --- --- --- --- --- --- ----

    # (1) The error in the plasma density comes from the error in the fit coefficents for the conversion between I_probe and ADC
    # % error n_i = sqrt(a^2 * delta a^2 + delta b^2)

    # (2) to get delta a, delta b we can do a reduced chi-square analysis on the fit:
    # chi^2 = 1/nu SUM (f(x_i) - y_i)^2 / (sigmaX_i^2 + sigmaY_i^2)
    # here:
    # (i) f: a*ADC+ b
    # (ii) y_i: ith value of ln(I_probe)
    # (iii): sigmaX_i: error in ADC value
    # (iv): SigmaY_i: error in ln(I_probe)
    #
    # we can work out what delta a, delta b are from optimization calculus to get:
    # (i): delta a = Sxx/DELTA
    # (ii): delta b = S/DELTA
    # for S = SUM(1/sqrt(sigmaX_i^2 +sigmaY_i^2)), Sxx = SUM(x_i^2/sqrt(sigmaX_i^2 +sigmaY_i^2))

    # (3) we need delta ADC and delta Ln(I_probe).
    # (a) The best we can do is delta ADC = +/- 0.5 ADC
    # (b) for Ln(I_probe), we assume voltage error of +/- 0.01V and known resistance error or 1% we have"
    IprobePercentError = np.sqrt((0.01/5.05)**2 + 0.01**2) # = about 0.01
    # (c) we convert ^^^ to delta Ln(I_probe) by delta ln(I_probe)= \delta I_probe /Iprobe BUT delta I_probe ~ 0.01 * Iprobe / Iprobe ===> 0.01
    deltaLnProbe = IprobePercentError


    Sxx = sum([ analog_vals[i]**2 / np.sqrt((0.5**2) + deltaLnProbe**2) for i in range(len(analog_vals))])
    Sx = sum([analog_vals[i] / np.sqrt((0.5 ** 2) + deltaLnProbe**2) for i in range(len(analog_vals))])
    S = sum([1 / np.sqrt((0.5 ** 2) + deltaLnProbe**2) for i in range(len(analog_vals))])
    DELTA = S*Sxx - (Sx)**2
    delta_a = Sxx/DELTA
    delta_b = S/DELTA
    delta_n = np.array([np.sqrt(  (parameters[0]*delta_a)**2 + (delta_b)**2  )])

    data_dict = {**data_dict, **{'fixed_current': [caldCurrent, {'LABLAXIS': 'current',
                                                    'DEPEND_0': 'fixed_Epoch',
                                                    'DEPEND_1': None,
                                                    'DEPEND_2': None,
                                                    'FILLVAL': rocketAttrs.epoch_fillVal,
                                                    'FORMAT': 'E12.2',
                                                    'UNITS': fixedLPunits,
                                                    'VALIDMIN': caldCurrent.min(),
                                                    'VALIDMAX': caldCurrent.max(),
                                                    'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
    data_dict = {**data_dict, **{'ni_percent_error': [delta_n, {'LABLAXIS': 'ni_percent_error',
                                                                 'DEPEND_0': None,
                                                                 'DEPEND_1': None,
                                                                 'DEPEND_2': None,
                                                                 'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                 'FORMAT': 'E12.2',
                                                                 'UNITS': 'cm^-3',
                                                                 'VALIDMIN': 0,
                                                                 'VALIDMAX': 1,
                                                                 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]}}
    data_dict['fixed_Epoch'] = data_dict.pop('Epoch_ni') # rename to fixed
    Done(start_time)



    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        # remove uneeded data.
        removeThese = ['ni', 'ni_swept', 'ne_swept', "Epoch_ni_swept", 'Epoch_ne_swept', 'step', 'EXP_Current', '28V_Monitor', 'Epoch_monitors']
        for thing in removeThese:
            data_dict.pop(thing)

        # rename Boom_Monitors
        data_dict['Swept_Boom_Monitor'] = data_dict.pop('Boom_Monitor_1')
        data_dict['Swept_Boom_Monitor'][1]['DEPEND_0'] = "Epoch_Swept_Boom_Monitor"
        data_dict['Epoch_Swept_Boom_Monitor'] = data_dict.pop('Epoch_monitor_1')

        data_dict['Fixed_Boom_Monitor'] = data_dict.pop('Boom_Monitor_2')
        data_dict['Fixed_Boom_Monitor'][1]['DEPEND_0'] = "Epoch_Fixed_Boom_Monitor"
        data_dict['Epoch_Fixed_Boom_Monitor'] = data_dict.pop('Epoch_monitor_2')

        prgMsg('Creating output file')

        outputPath = f'{outputFolderPath}\\' + f'{fliers[wRocket - 4]}\{fileoutName}'

        globalAttrsMod['Descriptor'] = wInstr[1]
        stl.outputCDFdata(outputPath, data_dict,instrNam='Langmuir',globalAttrsMod=globalAttrsMod)

        Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
data_directory = f'{DataPaths.ACES_data_folder}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}{modifier}\*.cdf'

if len(glob(data_directory)) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    for file_idx in wFiles:
        L2_to_L3_langmuir_fixed(wRocket-4, file_idx, justPrintFileNames)