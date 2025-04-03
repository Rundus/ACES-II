# --- L0_to_L1.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert electrostatic analyzer data from sector counts to raw counts
# Truncates all data to only that past 17:20:00.
# Special case given to LP's s.t. its dataset starts on sfid = 0 when reducing dataset


# NOTE: I have aligned the sector counts to start on an sfid 0 and end with an sfid 39 BUT
#  that doesn't mean I've started my energy sweeps at the beginning. It takes 4 minorframes
#  to complete 1 energy record and there are 49 + 1 retrace energy records per instrument sweep.
#  It takes 360, 280 and 360 words for the eepaa, iepaa, and leesa to complete a sweep,
#  which means for each major frame eepaa, iepaa, leesa complete 7 sweeps with 10 words left over

# I have chosen to timestamp an energy sweep at the beginning of the Sweep instead of the end

# There are several energy steps that are likely unusable, they are:
# 21.05,   18.01,   15.40,   13.17,   11.27,    9.64,    8.24,  7.05
# These channels exhibit a "mirroring" effect that suggests these channels
# were not at their assigned voltage during a sweep. Best seen by looking
# at the alfven signature in the HF data

# remember that we don't want the retrace value in our subsets


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

from src.my_imports import *

start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files in
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
wFiles = [[1, 2, 3], [1, 2]]

# EEPAA: how many energy values not to keep, starting from the lowest values e.g. adjusts = 8 --> remove the bottom 8 values
energy_adjusts = [8, 0, 0]  # [EEPAA,IEPAA,LEESA]
countNoiseThresh = 2 # anything LESS than this value is zero'd out

# Truncates all data to everything past 17:20:00 or whatever you wish
truncate_target_time = dt.datetime(2022, 11, 20, 17, 20)

# output the data
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from gc import collect
from warnings import \
    filterwarnings  # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.

filterwarnings("ignore")


def L0_to_L1_ESA(wRocket, wFile, rocketFolderPath, justPrintFileNames):

    # --- ACES II Flight/Integration Data ---
    globalAttrsMod = ACESII.global_attributes[wRocket - 4]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L1'

    # Set the paths for the file names
    L0Files = glob(f'{rocketFolderPath}L0\{ACESII.fliers[wRocket - 4]}\*.cdf')
    L0_names = [ifile.replace(f'{rocketFolderPath}L0\{ACESII.fliers[wRocket - 4]}\\', '') for ifile in L0Files]
    L0_names_searchable = [ ifile.replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace('l0_', '').replace('_v00', '') for ifile in L0_names]
    dataFile_name = L0Files[wFile].replace(f'{rocketFolderPath}L0\{ACESII.fliers[wRocket - 4]}\\', '')
    fileoutName = dataFile_name.replace('l0', 'l1')

    if justPrintFileNames:
        for i, file in enumerate(L0Files):
            print('[{:.0f}] {:70s}{:5.1f} MB'.format(i, L0_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return

    print('\n')
    print(stl.color.UNDERLINE + f'Converting to L1 data for {dataFile_name}' + stl.color.END)
    print('[' + str(wFile) + ']   ' + str(round(getsize(L0Files[wFile]) / (10 ** 6), 1)) + 'MiB')

    #############################
    # --- LOAD THE tmCDF FILE ---
    #############################
    stl.prgMsg('Loading data from L0Files')
    data_dict_tmCDF = stl.loadDictFromFile(L0Files[wFile])
    stl.Done(start_time)

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in dataFile_name:
            wInstr = [index, instr]

    # --- prepare the output ---
    data_dict_output = {'Energy': [[], {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None,
                                                  'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2', 'UNITS': 'eV',
                                                  'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data',
                                                  'SCALETYP': 'log', 'LABLAXIS': 'Energy'}],
                        'geometric_factor': [[], {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None,
                                                   'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2',
                                                   'UNITS': 'cm!U2!N str ev ev!U-1!N', 'VALIDMIN': None,
                                                   'VALIDMAX': None, 'VAR_TYPE': 'support_data',
                                                   'SCALETYP': 'linear'}],
                        'Epoch': [[], {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None,
                                           'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2',
                                           'UNITS': 'ns', 'VALIDMIN': None, 'VALIDMAX': None,
                                           'VAR_TYPE': 'support_data', 'MONOTON': 'INCREASE',
                                           'TIME_BASE': 'J2000', 'TIME_SCALE': 'Terrestrial Time',
                                           'REFERENCE_POSITION': 'Rotating Earth Geoid',
                                           'SCALETYP': 'linear'}],
                        'Pitch_Angle': [[], {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None,
                                              'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2', 'UNITS': 'deg',
                                              'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data',
                                              'SCALETYP': 'linear', 'LABLAXIS': 'Pitch_Angle'}],
                        'counts':
                            [[], {'LABLAXIS': ACESII.ESA_names_lower_case[wInstr[0]],
                                  'DEPEND_0': 'Epoch_esa', 'DEPEND_1': 'Pitch_Angle',
                                  'DEPEND_2': 'Energy',
                                  'FILLVAL': 65535, 'FORMAT': 'E12.2',
                                  'UNITS': 'counts', 'VALIDMIN': 0,
                                  'VALIDMAX': ACESII.ESA_max_counts, 'VAR_TYPE': 'data',
                                  'SCALETYP': 'linear'}]
                         }


    # --- --- --- --- --- --- --- --- ---
    # --- Calculate Instrument Data ---
    # --- --- --- --- --- --- --- --- ---
    stl.prgMsg('\nCreating L1 instrument data')
    if wInstr[1] == 'eepaa':  # EEPAA
        data_dict_output['Pitch_Angle'][0] = np.array(ACESII.ESA_instr_sector_to_pitch[0])
        data_dict_output['Energy'][0] = deepcopy(np.array(ACESII.ESA_instr_Energy[0]))
        data_dict_output['geometric_factor'][0] = np.array(ACESII.ESA_geometric_factor_TRACERS_ACE[0])

    elif wInstr[1] == 'leesa':  # LEESA
        data_dict_output['Pitch_Angle'][0] = np.array(ACESII.ESA_instr_sector_to_pitch[1])
        data_dict_output['Energy'][0] = deepcopy(np.array(ACESII.ESA_instr_Energy[1]))
        data_dict_output['geometric_factor'][0] = np.array(ACESII.ESA_geometric_factor_TRACERS_ACE[1])

    elif wInstr[1] == 'iepaa':  # IEPAA
        data_dict_output['Pitch_Angle'][0] = np.array(ACESII.ESA_instr_sector_to_pitch[2][::-1])
        data_dict_output['Energy'][0] = deepcopy(np.array(ACESII.ESA_instr_Energy[2]))
        data_dict_output['geometric_factor'][0] = np.array(ACESII.ESA_geometric_factor_TRACERS_ACE[2])

    # --- --- --- --- --- ----
    # --- PROCESS ESA DATA ---
    # --- --- --- --- --- ----
    sweepLength = len(data_dict_output['Energy'][0]) + 1  # must be a total of 50 values in a sweep

    # --- find index locations of the start and end point of the sweeps ---
    sectorCounts = data_dict_tmCDF['Sector_Counts'][0]
    start_of_sweeps = np.where(data_dict_tmCDF['sweep_step'][0] == 49)[0]  # before its turned on, sweep_steps are just zeros
    sweeps_start = start_of_sweeps[0] + 1  # index for the start of the first sweep, beginning on a retrace
    sweeps_end = start_of_sweeps[-1] + 1
    sectorCounts_trimmed = sectorCounts[sweeps_start:sweeps_end]
    epoch_trimmed = data_dict_tmCDF['Epoch'][0][sweeps_start:sweeps_end]
    no_of_sweeps = int(len(sectorCounts_trimmed) / (sweepLength))

    ###############################
    # --- STORE THE OUTPUT DATA ---
    ###############################

    # --- get the ESA Epoch ---
    # downsample the ESA's epoch so it starts at the beginning of each sweep
    # NOTE: the sweeps for all instruments are low Energy to High Energy
    data_dict_output['Epoch'][0] = np.array([epoch_trimmed[i * sweepLength] for i in range(no_of_sweeps) ])

    # --- get the ESA Counts ---
    # define the counts variable
    data_dict_output['counts'][0] = np.zeros(shape=(no_of_sweeps, len(data_dict_output['Pitch_Angle'][0]), len(data_dict_output['Energy'][0])))

    # loop through all the sweeps to cover the entire sector variable
    for i in tqdm(range(no_of_sweeps)):

        # take subset of sector counts and epoch data. This will include all pitch angles
        # NOTE: subSet is [[E_1 for all pitches], [E_2 for all pitches], ... [E_49 for all pitches]  ]
        subSet = sectorCounts_trimmed[i * sweepLength: (1 + i) * sweepLength]

        # place the sector count data
        for j in (range(len(data_dict_output['Pitch_Angle'][0]))):
            # We don't want the retrace value, so we can ignore the first value in the subSet
            sweep = [int(subSet[k][j]) for k in (range(sweepLength)) if k not in [1]]
            data_dict_output['counts'][0][i][j] = sweep[::-1]  # Invert the order so that highest energies are first. This is for convention

    # --- Clean up the Counts Variable ---
    # Apply a Counts treshold to the data so that only counts > threshold
    data_dict_output['counts'][0][data_dict_output['counts'][0] < countNoiseThresh] = 0

    # --- --- --- --- --- --- ---
    # --- ADJUST PITCH ANGLES ---
    # --- --- --- --- --- --- ---
    # the iepaa had the sector number inverse to the pitch angle i.e. sector1 --> 180deg, sector 2 --> 150 deg, etc.
    # This is opposite to the eepaa and leesa. Here we correct it by reversing the order of the pitch.
    if wInstr[1] == 'iepaa':
        data_dict_output['counts'][0] = data_dict_output['counts'][0][:,::-1,:]

    # --- --- --- --- --- ---
    # --- ADJUST ENERGIES ---
    # --- --- --- --- --- ---
    # The lowest energies of our ESA experienced unwanted RC decay. Here we remove them and update the counts variable
    # Adjust data if energy_adjust != 0:
    if energy_adjusts[wInstr[0]]:
        stl.prgMsg('Adjusting Energy')
        energy_adjust = energy_adjusts[wInstr[0]]
        data_dict_output['Energy'][0] = np.array(data_dict_output['Energy'][0][0:-energy_adjust])  # reduce the energy
        data_dict_output['counts'][0] = data_dict_output['counts'][0][:, :, 0:-energy_adjust]
        stl.Done(start_time)

    #########################################
    # --- ADD HOUSEKEEPING DATA TO OUTPUT ---
    #########################################
    stl.prgMsg('Downsampling HousingKeeping Data')
    Epoch_output_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_output['Epoch'][0]])
    for key in ['STEPPER_Voltage_ADC', 'MCP_Voltage_ADC', 'MCP_Current_ADC', 'STACK_Voltage_ADC', 'STACK_Current_ADC',
                'Count_Interval', '625kHz_Clock_Input', '3p3V_Voltage_Monitor_ADC', '5V_Voltage_Monitor_ADC',
                'Temperature_Monitor_ADC', 'EXP_Current','HV_div16','HV_enable','sweep_step','TP5_enable']:

        # first trim the data
        trimmed_var = data_dict_tmCDF[key][0][sweeps_start:sweeps_end]
        data_dict_output = {
            **data_dict_output,
            **{key:[np.array([trimmed_var[i * sweepLength] for i in range(no_of_sweeps) ]),  data_dict_tmCDF[key][1]]}
        }

    Epoch_monitors_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_tmCDF['Epoch_monitors'][0]])
    downsampled_indicies = np.array([np.abs(Epoch_monitors_tt2000 - val).argmin() for val in Epoch_output_tt2000])
    for key in ['28V_Monitor','Boom_Monitor']:
        data_dict_output = {
            **data_dict_output,
            **{key:[ data_dict_tmCDF[key][0][downsampled_indicies], data_dict_tmCDF[key][1]  ]}
        }

    # Update Dependencies
    for key, val in data_dict_output.items():
        if key not in ['geometric_factor', 'Energy', 'Pitch_Angle', 'Sector_Number']:
            data_dict_output[key][1]['DEPEND_0'] = 'Epoch'
    stl.Done(start_time)

    # --- --- --- --- --- ---
    # --- TRUNCATE DATASET ---
    # --- --- --- --- --- ---
    stl.prgMsg('Truncating Data')
    # Nothing interesting happens before 17:20 in the data, find this point and only keep data after it
    truncate_idx = np.abs(data_dict_output['Epoch'][0] - truncate_target_time).argmin()

    for key, val in data_dict_output.items():
        if key not in ['geometric_factor', 'Energy', 'Pitch_Angle', 'Sector_Number']:
            data_dict_output[key][0] = data_dict_output[key][0][truncate_idx:]
    stl.Done(start_time)

    # --- --- --- --- --- --- --- --- -
    # --- ADD ATTITUDE DATA TO FILE ---
    # --- --- --- --- --- --- --- --- -
    stl.prgMsg('Interpolating Attitude data into ESA files')
    inputFile_attitude = glob(f'C:\Data\ACESII\\attitude\{ACESII.fliers[wRocket - 4]}\\*Attitude_Solution*')[0]
    data_dict_attitude = stl.loadDictFromFile(inputFile_attitude, wKeys=['Lat', 'Long', 'Alt', 'Lat_geom', 'Long_geom', 'Alt_geom', 'Epoch'])

    # the low flyer attitude data needs to be extended to 17:20:00
    if wRocket == 5:
        # find the average step size in the attitude data
        epochAttitude = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_attitude['Epoch'][0]])
        deltaT = np.average([epochAttitude[i + 1] - epochAttitude[i] for i in range(len(epochAttitude) - 1)])
        targetStart = pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 20, 00, 14453))
        NumOfPoints = int((epochAttitude[0] - targetStart) / deltaT) + 1

        for key, val in data_dict_attitude.items():
            if key != 'Epoch':
                newVals = [val[0][0] for i in range(NumOfPoints)]
                data_dict_attitude[key][0] = np.array(newVals + list(data_dict_attitude[key][0]))
            else:
                newTimes = [pycdf.lib.tt2000_to_datetime(int(targetStart + deltaT * i)) for i in range(NumOfPoints)]
                data_dict_attitude[key][0] = np.array(newTimes + list(data_dict_attitude[key][0]))

    # interpolate attitude data onto ESA dataset
    data_dict_attitudeInterp = stl.InterpolateDataDict(InputDataDict=data_dict_attitude,
                                                       InputEpochArray=data_dict_attitude['Epoch'][0],
                                                       targetEpochArray=data_dict_output['Epoch'][0], wKeys=[])

    for key, val in data_dict_attitudeInterp.items():
        val[1]['VAR_TYPE'] = 'support_data'
        data_dict_output = {**data_dict_output, **{key: val}}
    stl.Done(start_time)


    # --- --- --- --- --- --- --- ---
    # --- WRITE OUT THE ESA DATA ---
    # --- --- --- --- --- --- --- ---
    if outputData:
        stl.prgMsg('Creating output ESA file')
        outputPath = f'{rocketFolderPath}L1\{ACESII.fliers[wRocket - 4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output, globalAttrsMod=globalAttrsMod, instrNam=wInstr[1])
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 0:  # ACES II Integration High
    rocketFolderPath = DataPaths.Integration_data_folder
elif wRocket == 1:  # ACES II Integration Low
    rocketFolderPath = DataPaths.Integration_data_folder
elif wRocket == 4:  # ACES II High
    rocketFolderPath = DataPaths.ACES_data_folder
elif wRocket == 5:  # ACES II Low
    rocketFolderPath = DataPaths.ACES_data_folder

if len(glob(f'{rocketFolderPath}L0\{ACESII.fliers[wRocket - 4]}\*.cdf')) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    if justPrintFileNames:
        L0_to_L1_ESA(wRocket, 0, rocketFolderPath, justPrintFileNames)
    elif not wFiles[wRocket - 4]:
        for fileNo in (range(len(glob(f'{rocketFolderPath}L0\{ACESII.fliers[wRocket - 4]}\*.cdf')))):
            L0_to_L1_ESA(wRocket, fileNo, rocketFolderPath, justPrintFileNames)
    else:
        for filesNo in wFiles[wRocket - 4]:
            L0_to_L1_ESA(wRocket, filesNo, rocketFolderPath, justPrintFileNames)