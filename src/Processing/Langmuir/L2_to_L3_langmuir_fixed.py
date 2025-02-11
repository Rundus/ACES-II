# --- L2_to_L3_langmuir_fixed.py ---
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
from src.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False # Just print the names of files

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

# --- OutputData ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.Processing.Langmuir.toggles import L2toL3_FixedLPToggles

#######################
# --- MAIN FUNCTION ---
#######################
def L2_to_L3_langmuir_fixed(wflyer, wFile, justPrintFileNames):

    # --- FILE I/O ---
    rocket_folder_path = DataPaths.ACES_data_folder
    rocketID = ACESII.payload_IDs[wflyer]
    globalAttrsMod = deepcopy(ACESII.global_attributes[wflyer])
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'LangmuirData'

    inputFiles = glob(f'{rocket_folder_path}{L2toL3_FixedLPToggles.inputPath_modifier}\{ACESII.fliers[wflyer]}{L2toL3_FixedLPToggles.modifier}\*langmuir*')
    outputFiles = glob(f'{rocket_folder_path}{L2toL3_FixedLPToggles.outputPath_modifier}\{ACESII.fliers[wflyer]}\*Temp&Density*')

    input_names = [ifile.replace(f'{rocket_folder_path}{L2toL3_FixedLPToggles.inputPath_modifier}\{ACESII.fliers[wflyer]}{L2toL3_FixedLPToggles.modifier}\\', '') for ifile in inputFiles]
    output_names = [ofile.replace(f'{rocket_folder_path}{L2toL3_FixedLPToggles.outputPath_modifier}\{ACESII.fliers[wflyer]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace(L2toL3_FixedLPToggles.inputPath_modifier.lower() + '_', '').replace('_v00', '') for ifile in input_names]
    output_names_searchable = [ofile.replace(L2toL3_FixedLPToggles.outputPath_modifier.lower() + '_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = inputFiles[wFile].replace(f'{rocket_folder_path}{L2toL3_FixedLPToggles.inputPath_modifier}\{ACESII.fliers[wflyer]}{L2toL3_FixedLPToggles.modifier}\\','')

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB   Made L2: {:3s} '.format(i, input_names_searchable[i], round(os.path.getsize(file) / (10 ** 6), 1), anws[0]))
        return

    print('\n')
    print(stl.color.UNDERLINE + f'Converting to {L2toL3_FixedLPToggles.outputPath_modifier} data for {dataFile_name}' + stl.color.END)
    print('[' + str(wFile) + ']   ' + str(round(os.path.getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

    # --- get the data from the L2 file ---
    stl.prgMsg(f'Loading data from {L2toL3_FixedLPToggles.inputPath_modifier} Files')
    data_dict = stl.loadDictFromFile(inputFiles[wFile])
    data_dict['Epoch_swept_Current'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict['Epoch_swept_Current'][0][i]) for i in (range(len(data_dict['Epoch_swept_Current'][0])))])
    stl.Done(start_time)

    # --- get IRI data ---
    stl.prgMsg('Loading IRI data')
    IRI_path = glob(f'C:\Data\ACESII\science\IRI_ACESII_Slice\{ACESII.fliers[wflyer]}\*smoothed*')
    data_dict_IRI = stl.loadDictFromFile(IRI_path[0])
    stl.Done(start_time)

    ##################################################################
    # --- Calculate the plasma density from Ion Saturation Current ---
    ##################################################################
    stl.prgMsg('Calculating Fixed ni')

    # get the data
    fixed_current = data_dict['fixed_current'][0]

    # interpolate IRI model onto LP data
    data_dict_IRI_interp = stl.InterpolateDataDict(InputDataDict=data_dict_IRI,
                                               InputEpochArray=data_dict_IRI['Epoch'][0],
                                               wKeys=[],
                                               targetEpochArray=data_dict['fixed_Epoch'][0])

    # determining n_i from Ion saturation
    # using the fixed LP data (now calibrated), determine n_i from the basic ion saturation current equation
    if L2toL3_FixedLPToggles.fixed_Ti_assumed: # use fixed Ti and m_i
        vth_i = np.sqrt(2*(stl.q0*L2toL3_FixedLPToggles.Ti_assumed)/ (data_dict_IRI_interp['m_i_avg'][0]*np.pi))
        I_th = (1 / 2) * stl.q0 * ACESII.LP_probe_areas[wflyer][0] * vth_i
        ni_density_m3 = (fixed_current * (1E-9)) / (I_th * (1 + (ACESII.LP_fixed_probe_bias[wflyer] - L2toL3_FixedLPToggles.V_plas_assumed) / L2toL3_FixedLPToggles.fixed_Ti_assumed))

    else: # use IRI model
        Ti_eV = stl.kB* data_dict_IRI_interp['Ti'][0]/stl.q0
        vth_i = np.sqrt(2*(stl.q0*Ti_eV)/ (data_dict_IRI_interp['m_i_avg'][0]*np.pi))
        I_th = (1/2)*stl.q0*ACESII.LP_probe_areas[wflyer][0]*vth_i
        ni_density_m3 = (fixed_current*(1E-9)) / (I_th*(1 + (ACESII.LP_fixed_probe_bias[wflyer] - L2toL3_FixedLPToggles.V_plas_assumed)/Ti_eV))

    ni_density_cm3 = (1/(1E6))*ni_density_m3
    # ni_density_cm3[np.abs(ni_density_cm3) > 1E10] = ACESII.epoch_fillVal # remove outliers

    # create an output data_dict for the fixed data

    varAttrs = {'LABLAXIS': 'plasma density', 'DEPEND_0': 'Epoch',
                                                                   'DEPEND_1': None,
                                                                   'DEPEND_2': None,
                                                                   'FILLVAL': ACESII.epoch_fillVal,
                                                                   'FORMAT': 'E12.2',
                                                                   'UNITS': '!Ncm!A-3!N',
                                                                   'VALIDMIN': 0,
                                                                   'VALIDMAX': ni_density_cm3.max(),
                                                                   'VAR_TYPE': 'data', 'SCALETYP': 'linear'}

    data_dict_fixed = {'ni': [ni_density_cm3, varAttrs],
                       'Epoch': deepcopy(data_dict['fixed_Epoch']),
                       'Fixed_Boom_Monitor': deepcopy(data_dict['Fixed_Boom_Monitor']),
                       'Epoch_Fixed_Boom_Monitor': deepcopy(data_dict['Epoch_Fixed_Boom_Monitor'])
                       }


    # --- ---- QUALTIY ASSURANCE --- ----
    # some of the Epoch values are fillvals. Lets interpolate these datapoints

    # # find fillval indicies
    # fillVals_time = list(np.where(stl.dateTimetoTT2000(data_dict_fixed['Epoch'][0], inverse=False) <= 0)[0])
    # negative_ni = list(np.where(data_dict_fixed['ni'][0] <= 0)[0])
    # fillVals_ni = list(np.where(data_dict_fixed['ni'][0] == ACESII.epoch_fillVal)[0])
    # badIndicies = list(set(fillVals_time + negative_ni + fillVals_ni))
    #
    #
    # # find enormous jump gap indicies
    # data_dict_fixed['ni'][0] = np.delete(data_dict_fixed['ni'][0], obj=badIndicies)
    # data_dict_fixed['Epoch'][0] = np.delete(data_dict_fixed['Epoch'][0], obj=badIndicies)
    #
    # # --- ---- ERROR ANALYSIS --- ----
    # deltaNi = data_dict['ni_percent_error'][0][0]*np.array(data_dict_fixed['ni'][0])
    # data_dict_fixed = {**data_dict_fixed, **{'ni_error':[deltaNi,deepcopy(data_dict_fixed['ni'][1])]}}
    # stl.Done(start_time)


    if outputData:

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        stl.prgMsg('Creating output file')
        fileoutName_fixed = f'ACESII_{rocketID}_langmuir_fixed.cdf'
        outputPath = f'{rocket_folder_path}{L2toL3_FixedLPToggles.outputPath_modifier}\{ACESII.fliers[wflyer]}\\{fileoutName_fixed}'
        stl.outputCDFdata(outputPath, data_dict_fixed, instrNam= 'Langmuir Probe',globalAttrsMod=globalAttrsMod)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
data_directory = f'{DataPaths.ACES_data_folder}{L2toL3_FixedLPToggles.inputPath_modifier}\{ACESII.fliers[wRocket-4]}{L2toL3_FixedLPToggles.modifier}\*.cdf'

if len(glob(data_directory)) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    for file_idx in wFiles:
        L2_to_L3_langmuir_fixed(wRocket-4, file_idx, justPrintFileNames)

