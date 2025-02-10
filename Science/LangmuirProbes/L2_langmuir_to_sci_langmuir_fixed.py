# --- L2_langmuir_to_sci_langmuir_fixed.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Takes in the L2 Langmuir Data and
# (1) Converts fixed Fixed LP current to ni
# (2) performs the Chi Square analysis to get the Temperature and Density from the
# characteristic curves of the swept probe


# NOTE: to convert to plasma density, the average ion mass and temperature was assumed to be
# Ti_assumed = 1 eV OR uses the smoothed IRI model to estimate
# m_i = average mass from the smoothed IRI model


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

from myImports import *

start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

modifier = ''
inputPath_modifier = 'L2'  # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'L3\Langmuir'  # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
errorPath_modifier = 'calibration\LP_calibration'

# --- Fitting Toggles ---
unitConv = 1E9  # converts from A to nA

# --- OutputData ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from Science.LangmuirProbes.toggles import FixedLPToggles

##############################
# --- FITTED CHI FUNCTIONS ---
##############################
rocketAttrs, b, c = ACES_mission_dicts()
e_coefficient = (((q0 ** 3) * (rocketAttrs.LP_probe_areas[0][0] ** (2))) / (2 * np.pi * m_e)) ** (1 / 2)

def transitionFunc(x, a0, a1, a2):
    y = (e_coefficient) * a0 * np.exp((x - a1) / a2)
    return y
def saturationFunc(x, a0, a1, a2):
    y = a0*(e_coefficient)*(x/a2 + (1- a1/a2))
    return y


#######################
# --- MAIN FUNCTION ---
#######################
def L2_Langmuir_to_SciLangmuir(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer,wSweeps):
    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'LangmuirData'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*langmuir*')
    outputFiles = glob(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*Temp&Density*')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    output_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace(inputPath_modifier.lower() + '_', '').replace('_v00', '') for ifile in input_names]
    output_names_searchable = [ofile.replace(outputPath_modifier.lower() + '_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\','')
    fileoutName_fixed = f'ACESII_{rocketID}_langmuir_fixed.cdf'
    fileoutName_swept = f'ACESII_{rocketID}_langmuir_swept.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB   Made L2: {:3s} '.format(i, input_names_searchable[i], round(os.path.getsize(file) / (10 ** 6), 1), anws[0]))
        return

    print('\n')
    print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + color.END)
    print('[' + str(wFile) + ']   ' + str(round(os.path.getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

    # --- get the data from the L2 file ---
    prgMsg(f'Loading data from {inputPath_modifier} Files')
    data_dict = loadDictFromFile(inputFiles[wFile])
    if SECTION_SweptProbeniTe:
        data_dict['Epoch_swept_Current'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_swept_Current'][0][i]) for i in (range(len(data_dict['Epoch_swept_Current'][0])))])
    Done(start_time)

    # --- get IRI data ---
    prgMsg('Loading IRI data')
    IRI_path = glob(f'C:\Data\ACESII\science\IRI_ACESII_Slice\{fliers[wflyer]}\*smoothed*')
    data_dict_IRI = loadDictFromFile(IRI_path[0])
    Done(start_time)

    # --- get error data from calibration ---
    prgMsg(f'Loading ERROR data from {inputPath_modifier} Files')
    errorsPath = glob(f'{ACES_data_folder}{errorPath_modifier}\{fliers[wflyer]}\*Iowa*')
    data_dict_errors = loadDictFromFile(errorsPath[0])
    Done(start_time)

    if SECTION_calculateFixedni:

        prgMsg('Calculating Fixed ni')
        ##################################################################
        # --- Calculate the plasma density from Ion Saturation Current ---
        ##################################################################

        ni = data_dict['fixed_current'][0]

        # interpolate IRI model onto LP data
        data_dict_IRI_interp = InterpolateDataDict(InputDataDict=data_dict_IRI,
                                                   InputEpochArray=data_dict_IRI['Epoch'][0],
                                                   wKeys=[],
                                                   targetEpochArray=data_dict['fixed_Epoch'][0])

        # determining n_i from Ion saturation
        # using the fixed LP data (now calibrated), determine n_i from the basic ion saturation current equation
        if fixedTi_assumed: # use fixed Ti and m_i
            vth_i = [np.sqrt( 2*(q0*Ti_assumed)/(data_dict_IRI_interp['m_i_avg'][0][k])) for k in range(len(ni))]
        else: # use IRI model
            vth_i = [np.sqrt((kB*data_dict_IRI_interp['Ti'][0][k])/(data_dict_IRI_interp['m_i_avg'][0][k])) for k in range(len(ni))]

        ni_density = (1/(1E6))*np.array([(4 * (current*(1E-9))) / (q0 * vth_i[h] * rocketAttrs.LP_probe_areas[0][0]) for h, current in enumerate(ni)])
        ni_density[np.abs(ni_density) > 1E10] = rocketAttrs.epoch_fillVal # remove outliers

        # create an output data_dict for the fixed data

        varAttrs = {'LABLAXIS': 'plasma density', 'DEPEND_0': 'Epoch',
                                                                       'DEPEND_1': None,
                                                                       'DEPEND_2': None,
                                                                       'FILLVAL': rocketAttrs.epoch_fillVal,
                                                                       'FORMAT': 'E12.2',
                                                                       'UNITS': '!Ncm!A-3!N',
                                                                       'VALIDMIN': 0,
                                                                       'VALIDMAX': ni_density.max(),
                                                                       'VAR_TYPE': 'data', 'SCALETYP': 'linear'}

        data_dict = {**data_dict, **{'ni': [ni_density, ]}}
        data_dict_fixed = {'ni': [ni_density, varAttrs],
                           'Epoch': data_dict['fixed_Epoch'],
                           'Fixed_Boom_Monitor': data_dict['Fixed_Boom_Monitor'],
                           'Epoch_Fixed_Boom_Monitor': data_dict['Epoch_Fixed_Boom_Monitor']}


        # --- ---- QUALTIY ASSURANCE --- ----
        # some of the Epoch values are fillvals. Lets interpolate these datapoints

        # find fillval indicies
        fillVals_time = list(np.where(dateTimetoTT2000(data_dict_fixed['Epoch'][0], inverse=False) <= 0)[0])
        negative_ni = list(np.where(data_dict_fixed['ni'][0] <= 0)[0])
        fillVals_ni = list(np.where(data_dict_fixed['ni'][0] == rocketAttrs.epoch_fillVal)[0])
        badIndicies = list(set(fillVals_time + negative_ni + fillVals_ni))


        # find enormous jump gap indicies
        data_dict_fixed['ni'][0] = np.delete(data_dict_fixed['ni'][0], obj=badIndicies)
        data_dict_fixed['Epoch'][0] = np.delete(data_dict_fixed['Epoch'][0], obj=badIndicies)

        if tromsoCal:
            data_dict_fixed['ni'][0] = np.array(data_dict_fixed['ni'][0]) * tromsoScales[wRocket-4]

        # --- ---- ERROR ANALYSIS --- ----
        deltaNi = data_dict['ni_percent_error'][0][0]*np.array(data_dict_fixed['ni'][0])
        data_dict_fixed = {**data_dict_fixed, **{'ni_error':[deltaNi,deepcopy(data_dict_fixed['ni'][1])]}}
        Done(start_time)

    if outputData:

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        prgMsg('Creating output file')


        if SECTION_calculateFixedni:

            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName_fixed}'
            outputCDFdata(outputPath, data_dict_fixed, instrNam= 'Langmuir Probe')

        if SECTION_SweptProbeniTe:
            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName_swept}'
            outputCDFdata(outputPath,data_dict_swept,instrNam='Langmuir Probe')

        Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5:  # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        L2_Langmuir_to_SciLangmuir(wRocket, 0, rocketFolderPath, justPrintFileNames, wflyer,wSweeps)
    else:
        L2_Langmuir_to_SciLangmuir(wRocket, 0, rocketFolderPath, justPrintFileNames, wflyer,wSweeps)
