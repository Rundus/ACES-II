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

# --- OutputData ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.Processing.Langmuir.toggles import L2toL3_FixedLPToggles as fToggles

#######################
# --- MAIN FUNCTION ---
#######################
def L2_to_L3_langmuir_fixed(wflyer, justPrintFileNames):

    # --- FILE I/O ---
    rocket_folder_path = DataPaths.ACES_data_folder
    rocketID = ACESII.payload_IDs[wflyer]
    globalAttrsMod = deepcopy(ACESII.global_attributes[wflyer])
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'LangmuirData'

    inputFiles = glob(f'{rocket_folder_path}{fToggles.inputPath_modifier}\{ACESII.fliers[wflyer]}{fToggles.modifier}\*langmuir_fixed*')
    input_names = [ifile.replace(f'{rocket_folder_path}{fToggles.inputPath_modifier}\{ACESII.fliers[wflyer]}{fToggles.modifier}\\', '') for ifile in inputFiles]
    input_names_searchable = [ifile.replace(fToggles.inputPath_modifier.lower() + '_', '').replace('_v00', '') for ifile in input_names]
    dataFile_name = inputFiles[0].replace(f'{rocket_folder_path}{fToggles.inputPath_modifier}\{ACESII.fliers[wflyer]}{fToggles.modifier}\\','')

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(os.path.getsize(file) / (10 ** 6), 1)))
        return

    print('\n')
    print(stl.color.UNDERLINE + f'Converting to {fToggles.outputPath_modifier} data for {dataFile_name}' + stl.color.END)
    print('[' + str(0) + ']   ' + str(round(os.path.getsize(inputFiles[0]) / (10 ** 6), 1)) + 'MiB')

    # --- get the data from the L2 file ---
    stl.prgMsg(f'Loading data from {fToggles.inputPath_modifier} Files')
    data_dict = stl.loadDictFromFile(inputFiles[0])
    stl.Done(start_time)

    # --- get IRI data ---
    stl.prgMsg('Loading IRI data')
    IRI_path = glob(f'C:\Data\ACESII\science\IRI_ACESII_Slice\{ACESII.fliers[wflyer]}\*smoothed*')
    data_dict_IRI = stl.loadDictFromFile(IRI_path[0])
    stl.Done(start_time)

    ##################################################################
    # --- Calculate the plasma density from Ion Saturation Current ---
    ##################################################################

    stl.prgMsg('Interpolating IRI model')

    # interpolate IRI model onto LP data
    data_dict_IRI_interp = stl.InterpolateDataDict(InputDataDict=data_dict_IRI,
                                               InputEpochArray=data_dict_IRI['Epoch'][0],
                                               wKeys=['Ti', 'm_i_avg'],
                                               targetEpochArray=data_dict['Epoch'][0])
    stl.Done(start_time)

    stl.prgMsg('Calculating Fixed ni')
    # get the data
    fixed_current = data_dict['fixed_current'][0]  # note: current is in Nano Amps

    # determining n_i from Ion saturation
    # using the fixed LP data (now calibrated), determine n_i from the basic ion saturation current equation
    if fToggles.fixed_Ti_assumed: # use fixed Ti and m_i
        vth_i = np.array(np.sqrt(2*(stl.q0*fToggles.Ti_assumed)/ (data_dict_IRI_interp['m_i_avg'][0]*np.pi)))
        I_th = (1 / 2) * stl.q0 * ACESII.LP_probe_areas[wflyer][0] * vth_i
        ni_density_m3 = (-1*fixed_current * (1E-9)) / (I_th * (1 + (ACESII.LP_fixed_probe_bias[wflyer] - fToggles.V_plas_assumed) / fToggles.fixed_Ti_assumed))

    else: # use IRI model
        Ti_eV = stl.kB* data_dict_IRI_interp['Ti'][0]/stl.q0
        vth_i = np.array(np.sqrt(2*(stl.q0*Ti_eV)/ (data_dict_IRI_interp['m_i_avg'][0]*np.pi)))
        I_th = (1/2)*stl.q0*ACESII.LP_probe_areas[wflyer][0]*vth_i
        ni_density_m3 = (-1*fixed_current*(1E-9)) / (I_th*(1 + (ACESII.LP_fixed_probe_bias[wflyer] - fToggles.V_plas_assumed[wflyer])/Ti_eV))

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
                       'Epoch': deepcopy(data_dict['Epoch']),
                       'fixed_boom_monitor': deepcopy(data_dict['fixed_boom_monitor']),
                       'Epoch_boom_monitor': deepcopy(data_dict['Epoch_boom_monitor'])
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
    stl.Done(start_time)


    if outputData:

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        stl.prgMsg('Creating output file')
        fileoutName_fixed = f'ACESII_{rocketID}_l3_langmuir_fixed.cdf'
        outputPath = f'{rocket_folder_path}{fToggles.outputPath_modifier}\{ACESII.fliers[wflyer]}\\{fileoutName_fixed}'
        stl.outputCDFdata(outputPath, data_dict_fixed, instrNam= 'Langmuir Probe',globalAttrsMod=globalAttrsMod)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
data_directory = f'{DataPaths.ACES_data_folder}{fToggles.inputPath_modifier}\{ACESII.fliers[wRocket-4]}{fToggles.modifier}\*.cdf'

if len(glob(data_directory)) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    L2_to_L3_langmuir_fixed(wRocket-4, justPrintFileNames)

