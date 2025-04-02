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
energy_adjusts = [8, 0, 0] #  [EEPAA,IEPAA,LEESA]
countNoiseThresh = 2

# Truncates all data to everything past 17:20:00 or whatever you wish
truncateData = True
targetTruncDate, targetTruncTime = [2022, 11, 20], [17, 20, 0, 0]

# output the data
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from gc import collect
from warnings import filterwarnings # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.
filterwarnings("ignore")


def L0_to_L1(wRocket, wFile, rocketFolderPath, justPrintFileNames):

    # --- ACES II Flight/Integration Data ---
    globalAttrsMod = ACESII.global_attributes[wRocket-4]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L1'
    # L1ModelData = ACESII.L1_TRICE_Quick(wRocket-4)

    # Set the paths for the file names
    L0Files = glob(f'{rocketFolderPath}L0\{ACESII.fliers[wRocket-4]}\*.cdf')
    L1Files = glob(f'{rocketFolderPath}L1\{ACESII.fliers[wRocket-4]}\*.cdf')

    L0_names = [ifile.replace(f'{rocketFolderPath}L0\{ACESII.fliers[wRocket-4]}\\', '') for ifile in L0Files]
    L1_names = [ofile.replace(f'{rocketFolderPath}L1\{ACESII.fliers[wRocket-4]}\\', '') for ofile in L1Files]

    L0_names_searchable = [ifile.replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace('l0_', '').replace('_v00', '') for ifile in L0_names]
    L1_names_searchable = [ofile.replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace('l1_', '').replace('_v00','').replace('__', '_') for ofile in L1_names]

    dataFile_name = L0Files[wFile].replace(f'{rocketFolderPath}L0\{ACESII.fliers[wRocket-4]}\\', '')
    fileoutName = dataFile_name.replace('l0', 'l1')

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in dataFile_name:
            wInstr = [index, instr]


    if justPrintFileNames:
        for i, file in enumerate(L0Files):
            anws = ["yes" if L0_names_searchable[i] in L1_names_searchable else "no"]
            print('[{:.0f}] {:70s}{:5.1f} MB'.format(i, L0_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return

    print('\n')
    print(stl.color.UNDERLINE + f'Converting to L1 data for {dataFile_name}' + stl.color.END)
    print('[' + str(wFile) + ']   ' + str(round(getsize(L0Files[wFile]) / (10 ** 6), 1)) + 'MiB')

    # --- get the data from the tmCDF file ---
    stl.prgMsg('Loading data from L0Files')

    data_dict = stl.loadDictFromFile(L0Files[wFile])

    stl.Done(start_time)


    # --- --- --- --- --- --- --- --- ---
    # --- Calculate Instrument Data ---
    # --- --- --- --- --- --- --- --- ---

    stl.prgMsg('\nCreating L1 instrument data')

    # 0 -> EEPAA
    # 1 -> LEESA
    # 2 -> IEPAA
    # 3 -> LP
    if wInstr[0] in [0, 1, 2]:
        data_dict = {**data_dict, **{'Energy':           [[0], {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2', 'UNITS': 'eV', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'log', 'LABLAXIS': 'Energy'}]}}
        data_dict = {**data_dict, **{'geometric_factor': [[0], {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2', 'UNITS': 'cm!U2!N str ev ev!U-1!N', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]}}
        data_dict = {**data_dict, **{'Epoch_esa':        [[],  {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2', 'UNITS': 'ns', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000', 'TIME_SCALE': 'Terrestrial Time', 'REFERENCE_POSITION': 'Rotating Earth Geoid', 'SCALETYP': 'linear'}]}}
        data_dict['Pitch_Angle'] = [[0], {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2', 'UNITS': 'deg', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear', 'LABLAXIS': 'Pitch_Angle'}]
        data_dict['28V_Monitor'][1]['DEPEND_0'] = 'Epoch_monitors'
        data_dict['Boom_Monitor'][1]['DEPEND_0'] = 'Epoch_monitors'
        sectorCounts = data_dict['Sector_Counts'][0]

        if wInstr[1] == 'eepaa': # EEPAA
            data_dict['Pitch_Angle'][0] = np.array(ACESII.ESA_instr_sector_to_pitch[0])
            pitches = data_dict['Pitch_Angle'][0]
            data_dict['Energy'][0] =  deepcopy(np.array(ACESII.ESA_instr_Energy[0]))
            Energy = deepcopy(data_dict['Energy'][0])
            data_dict['geometric_factor'][0] = np.array(ACESII.ESA_geometric_factor_TRACERS_ACE[0])

        elif wInstr[1] == 'leesa': # LEESA
            data_dict['Pitch_Angle'][0] = np.array(ACESII.ESA_instr_sector_to_pitch[1])
            pitches = data_dict['Pitch_Angle'][0]
            data_dict['Energy'][0] = deepcopy(np.array(ACESII.ESA_instr_Energy[1]))
            Energy = deepcopy(data_dict['Energy'][0])
            data_dict['geometric_factor'][0] = np.array(ACESII.ESA_geometric_factor_TRACERS_ACE[1])

        elif wInstr[1] == 'iepaa': # IEPAA
            data_dict['Pitch_Angle'][0] = np.array(ACESII.ESA_instr_sector_to_pitch[2][::-1])
            pitches = data_dict['Pitch_Angle'][0]
            data_dict['Energy'][0] = deepcopy(np.array(ACESII.ESA_instr_Energy[2]))
            Energy = deepcopy(data_dict['Energy'][0])
            data_dict['geometric_factor'][0] = np.array(ACESII.ESA_geometric_factor_TRACERS_ACE[2])

        # --- PROCESS ESA DATA ---
        sweepLength = len(data_dict['Energy'][0]) + 1 # must be a total of 50 values in a sweep

        # --- find index locations of the start and end point of the sweeps ---
        start_of_sweeps = np.where(data_dict['sweep_step'][0] == 49)[0] # before its turned on, sweep_steps are just zeros
        sweeps_start = start_of_sweeps[0] + 1 # index for the start of the first sweep, beginning on a retrace
        sweeps_end = start_of_sweeps[-1] + 1
        sectorCounts_reduced = sectorCounts[sweeps_start:sweeps_end]
        epoch_reduced = data_dict['Epoch'][0][sweeps_start:sweeps_end]
        no_of_sweeps = int(len(sectorCounts_reduced) / (sweepLength))

        # --- PLACE THE OUTPUT DATA ---
        counts = np.zeros(shape=(no_of_sweeps, len(pitches), len(data_dict['Energy'][0])))

        # loop through all the sweeps to cover the entire sector variable
        for i in tqdm(range(no_of_sweeps)):

            # --- downsample the ESA's epoch so it starts at the beginning of each sweep ---
            # NOTE: the sweeps for all instruments are low Energy to High Energy
            data_dict['Epoch_esa'][0].append(epoch_reduced[i*sweepLength])

            # --- take subset of sector counts and epoch data. This will include all pitch angles ---
            # NOTE: subSet is [[E_1 for all pitches], [E_2 for all pitches], ... [E_49 for all pitches]  ]
            subSet = sectorCounts_reduced[i * sweepLength: (1 + i) * sweepLength]

            # --- place the sector count data ---
            for j in (range(len(pitches))):

                # We don't want the retrace value, so we can ignore the first value in the subSet
                sweep = [int(subSet[k][j]) for k in (range(sweepLength)) if k not in [1]]
                counts[i][j] = sweep[::-1] # Invert the order so that highest energies are first. This is for convention

        # the iepaa had the sector number inverse to the pitch angle i.e. sector1 --> 180deg, sector 2 --> 150 deg, etc.
        # This is opposite to the eepaa and leesa. Here we correct it by reversing the order of the pitch.
        if wInstr[1] =='iepaa':
            counts = [counts[tme][::-1] for tme in range(len(counts))]
        # Adjust data if energy_adjust != 0:
        if wInstr[1] =='eepaa':
            energy_adjust = energy_adjusts[0]
        elif wInstr[1] =='iepaa':
            energy_adjust = energy_adjusts[1]
        elif wInstr[1] =='leesa':
            energy_adjust = energy_adjusts[2]
        else:
            energy_adjust = 0

        if energy_adjust != 0:
            stl.prgMsg('Adjusting Energy')
            data_dict['Energy'][0] = np.array(Energy[0:-energy_adjust]) # reduce the energy
            countsNew = [[[0 for i in range(len(data_dict['Energy'][0]))] for ptch in range(len(pitches))] for tme in range(len(counts))]
            counts = list(counts) # convert counts into a dynamic variable
            for tme in range(len(counts)):
                for ptch in range(len(pitches)):
                    countsNew[tme][ptch] = counts[tme][ptch][0:-energy_adjust]
            counts = np.array(countsNew)
            stl.Done(start_time)
        else:
            counts = np.array(counts)

        # --- Apply a Counts treshold to the data so that only counts > 3 exist ---
        finder = np.where(np.abs(counts) <= countNoiseThresh)
        counts[finder] = 0

        data_dict['Epoch_esa'][0] = np.array(data_dict['Epoch_esa'][0])
        data_dict = {**data_dict, **{ACESII.ESA_names_lower_case[wInstr[0]]:
                                         [counts, {'LABLAXIS': ACESII.ESA_names_lower_case[wInstr[0]],
                                                   'DEPEND_0': 'Epoch_esa', 'DEPEND_1': 'Pitch_Angle',
                                                   'DEPEND_2': 'Energy',
                                                   'FILLVAL': 65535, 'FORMAT': 'E12.2',
                                                   'UNITS': 'counts', 'VALIDMIN': 0,
                                                   'VALIDMAX': ACESII.ESA_max_counts, 'VAR_TYPE': 'data',
                                                   'SCALETYP': 'linear'}]}
                     }

        del data_dict['Sector_Counts'], data_dict['sfid'], data_dict['sweep_step'], data_dict['minor_frame_counter'], data_dict['major_frame_counter']

        # --- --- --- --- --- --- --- --- -
        # --- ADD ATTITUDE DATA TO FILE ---
        # --- --- --- --- --- --- --- --- -
        stl.prgMsg('Truncating Data')

        if truncateData:

            # --- --- --- --- --- ---
            # --- REDUCE DATASET ---
            # --- --- --- --- --- ---
            # Nothing interesting happens before 17:20 in the data, find this point and only keep data after it

            targetDate, targetTime  = dt.datetime(targetTruncDate[0], targetTruncDate[1], targetTruncDate[2]), dt.time(targetTruncTime[0], targetTruncTime[1], targetTruncTime[2], 0)
            targetDateTime_TT2000 = pycdf.lib.datetime_to_tt2000(targetDate.combine(targetDate, targetTime))
            Epoch_targetIndex = np.array([np.abs(pycdf.lib.datetime_to_tt2000(time) - targetDateTime_TT2000) for time in data_dict['Epoch'][0]]).argmin()
            Epoch_esa_targetIndex = np.array([np.abs(pycdf.lib.datetime_to_tt2000(time) - targetDateTime_TT2000) for time in data_dict['Epoch_esa'][0]]).argmin()
            Epoch_monitors_targetIndex = np.array([np.abs(pycdf.lib.datetime_to_tt2000(time) - targetDateTime_TT2000) for time in data_dict['Epoch_monitors'][0]]).argmin()

            noReduction = ['geometric_factor', 'Energy', 'Pitch_Angle', 'Sector_Number']
            needsReduction_esa = [ACESII.ESA_names_lower_case[wInstr[0]], 'Epoch_esa']
            needsReduction_monitors = ['28V_Monitor', 'Boom_Monitor', 'Epoch_monitors']

            for key, val in data_dict.items():
                if key in needsReduction_esa: # counts data
                    data_dict[key][0] = data_dict[key][0][Epoch_esa_targetIndex:]
                elif key in needsReduction_monitors: # Monitors data
                    data_dict[key][0] = data_dict[key][0][Epoch_monitors_targetIndex:]
                elif key not in noReduction: # Everything else
                    data_dict[key][0] = data_dict[key][0][Epoch_targetIndex:]

        # Change the name of the Epoch variable to ---> Epoch Housekeeping
        data_dict['Epoch_HouseKeeping'] = data_dict.pop('Epoch')

        # Change the name of the Epoch_esa variable to ---> Epoch
        data_dict['Epoch'] = data_dict.pop('Epoch_esa')

        # Update Dependencies
        data_dict[wInstr[1]][1]['DEPEND_0'] = 'Epoch'
        stl.Done(start_time)

        stl.prgMsg('Interpolating Attitude data into ESA files')
        inputFile_attitude = glob(f'C:\Data\ACESII\\attitude\{ACESII.fliers[wRocket-4]}\\*Attitude_Solution*')[0]
        data_dict_attitude = stl.loadDictFromFile(inputFile_attitude, wKeys=['Lat', 'Long', 'Alt', 'Lat_geom', 'Long_geom', 'Alt_geom', 'Epoch'])

        # the low flyer attitude data needs to be extended to 17:20:00
        if wRocket == 5:

            # find the average step size in the attitude data
            epochAttitude = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_attitude['Epoch'][0]])
            deltaT = np.average([epochAttitude[i+1] - epochAttitude[i] for i in range(len(epochAttitude)-1)])
            targetStart = pycdf.lib.datetime_to_tt2000(dt.datetime(2022,11,20,17,20,00,14453))
            NumOfPoints = int((epochAttitude[0]-targetStart)/deltaT)+1

            for key, val in data_dict_attitude.items():
                if key != 'Epoch':
                    newVals = [val[0][0] for i in range(NumOfPoints)]
                    data_dict_attitude[key][0] = np.array(newVals + list(data_dict_attitude[key][0]))
                else:
                    newTimes = [pycdf.lib.tt2000_to_datetime(int(targetStart + deltaT*i)) for i in range(NumOfPoints)]
                    data_dict_attitude[key][0] = np.array(newTimes + list(data_dict_attitude[key][0]))

        # interpolate attitude data onto ESA dataset
        data_dict_attitudeInterp = stl.InterpolateDataDict(InputDataDict=data_dict_attitude,
                                                       InputEpochArray=data_dict_attitude['Epoch'][0],
                                                       targetEpochArray=data_dict['Epoch'][0], wKeys=[])

        for key, val in data_dict_attitudeInterp.items():
            val[1]['VAR_TYPE'] = 'support_data'
            data_dict = {**data_dict, **{key: val}}
        stl.Done(start_time)

        # --- --- --- --- --- --- --- ---
        # --- WRITE OUT THE ESA DATA ---
        # --- --- --- --- --- --- --- ---
        if outputData:
            stl.prgMsg('Creating output ESA file')

            outputPath = f'{rocketFolderPath}L1\{ACESII.fliers[wRocket-4]}\\{fileoutName}'

            # stl.outputCDFdata(outputPath, data_dict,ModelData=L1ModelData, globalAttrsMod= globalAttrsMod, instrNam= wInstr[1])
            stl.outputCDFdata(outputPath, data_dict, globalAttrsMod=globalAttrsMod, instrNam=wInstr[1])

        stl.Done(start_time)



    elif wInstr[0] in [3]:

        # --- --- --- --- --- --- --- --- --- ---
        # --- SEPARATE AND COLLECT LP DATA ---
        # --- --- --- --- --- --- --- --- --- ---
        try:
            # --- Get L0 data ---
            L0epoch = [pycdf.lib.datetime_to_tt2000(data_dict['Epoch'][0][i]) for i in range(len(data_dict['Epoch'][0]))]
            L0ChannelCounts = data_dict['Channel_Counts'][0]
            L0sfid = data_dict['sfid'][0]
            samplePeriod = ACESII.LPSamplePeriod
            minorFrameTime = ACESII.MinorFrameTime # each LP value is really sampled 1 minorframe before its reported value in the PCM matrix

            del data_dict['Epoch'], data_dict['Channel_Counts'],data_dict['Channel_Number'], data_dict['minor_frame_counter'],data_dict['major_frame_counter'], data_dict['sfid']
            collect()

            # --- --- --- --- --- ---
            # --- REDUCE DATASET ---
            # --- --- --- --- --- ---

            if truncateData:
                # Nothing interesting happens before 17:20 in the data, find this point and only keep data after it
                targetDate, targetTime = dt.datetime(2022, 11, 20), dt.time(17, 20, 0, 0)
                targetDateTime_TT2000 = pycdf.lib.datetime_to_tt2000(targetDate.combine(targetDate, targetTime))
                Epoch_targetIndex = np.array([np.abs(time - targetDateTime_TT2000) for time in L0epoch]).argmin()

                # --- Backtrack to where sfid = 0 ---
                dataset_targetIndex = 0
                for i in range(0,len(L0sfid)):
                    if L0sfid[Epoch_targetIndex-1*i] == 0:
                        dataset_targetIndex = Epoch_targetIndex-1*i
                        break
            else:
                # if we don't truncate, just find where sfid = 0 and cut everything else
                dataset_targetIndex = np.where(L0sfid == 0)[0][0]

            # --- reduce the dataset ---
            L0epoch = L0epoch[dataset_targetIndex:]
            L0ChannelCounts = L0ChannelCounts[dataset_targetIndex:]
            L0sfid = L0sfid[dataset_targetIndex:]


            # create the 5 variables in the LP data_dict: deltaNdivN,ni,ne_swept,ni_swept,step
            for i in range(len(ACESII.LP_Variables)):
                data_dict = {**data_dict, **{f'Epoch_{ACESII.LP_Variables[i]}': [[], {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -9223372036854775808, 'FORMAT': 'E12.2', 'UNITS': 'ns', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000', 'TIME_SCALE': 'Terrestrial Time', 'REFERENCE_POSITION': 'Rotating Earth Geoid', 'SCALETYP': 'linear'}]}}
                data_dict = {**data_dict, **{f'{ACESII.LP_Variables[i]}': [[], {'LABLAXIS': ACESII.LP_Variables[i], 'DEPEND_0': f"Epoch_{ACESII.LP_Variables[i]}", 'DEPEND_1':None, 'DEPEND_2':None, 'FILLVAL': -1, 'FORMAT': 'E12.2', 'UNITS': 'Digital', 'VALIDMIN': 0, 'VALIDMAX': ACESII.esaMaxCounts, 'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}

            no_of_fillvals = []

            for j in tqdm(range(len(L0epoch))):

                LP_var_counter = L0sfid[j]%4 + 1 # sets which LP_var I'm on. (step, ne_swept, ni_swept, ni)

                # deltaNdivN takes the first 8 values
                data_dict[ACESII.LP_Variables[0]][0].append(L0ChannelCounts[j][0:8])

                if L0epoch[j] != ACESII.epoch_fillVal:

                    # assign epoch values for each point with a sample rate of 35.25kHz
                    data_dict[f'Epoch_{ACESII.LP_Variables[0]}'][0].append([(L0epoch[j] - minorFrameTime) + l * samplePeriod for l in range(8)])

                    if L0ChannelCounts[j][8] > 4095: # set all values above 4095 to the fillval
                        data_dict[ACESII.LP_Variables[LP_var_counter]][0].append(-1)
                    else:
                        data_dict[ACESII.LP_Variables[LP_var_counter]][0].append(L0ChannelCounts[j][8]) # LP_vars

                    data_dict[f'Epoch_{ACESII.LP_Variables[LP_var_counter]}'][0].append((L0epoch[j] - minorFrameTime) + 8 * samplePeriod)

                else: # if epoch is a fillval
                    no_of_fillvals.append([ACESII.LP_Variables[LP_var_counter],j])
                    data_dict[f'Epoch_{ACESII.LP_Variables[0]}'][0].append([ACESII.epoch_fillVal for l in range(8)]) # deltaNdivN
                    data_dict[ACESII.LP_Variables[LP_var_counter]][0].append(-1)  # LP_vars
                    data_dict[f'Epoch_{ACESII.LP_Variables[LP_var_counter]}'][0].append(ACESII.epoch_fillVal) # LP_vars


            # Converts all data to numpy arrays
            for key, var in data_dict.items():
                data_dict[key][0] = np.array(data_dict[key][0])

            data_dict['deltaNdivN'][0] = np.array(data_dict['deltaNdivN'][0]).flatten()
            data_dict['Epoch_deltaNdivN'][0] = np.array(data_dict['Epoch_deltaNdivN'][0]).flatten()
            dataFailed = False


        except Exception as e:
            dataFailed = True
            print(stl.color.RED + f"{e}" + stl.color.END)

        # --- --- --- --- --- ---
        # --- WRITE OUT DATA ---
        # --- --- --- --- --- ---

        if not dataFailed:

            stl.prgMsg('\nCreating output LP file')

            # LP data is too large to be dumped all at once, must break it up
            for i in range(len(ACESII.LP_Variables)):

                # determine the outputPath and the variables which need to be output'd
                outputPath = f'{rocketFolderPath}L1\{ACESII.fliers[wRocket-4]}\\{fileoutName.replace("lp", f"lp_{ACESII.LP_Variables[i]}")}'
                noWriteVarst = ["deltaNdivN","step", "ne_swept", "ni_swept","ni"]
                noWriteVarst.pop(i)
                noWriteVars = [thing for thing in noWriteVarst]

                for thing in noWriteVarst:
                    noWriteVars.append(f'Epoch_{thing}')

                writeOutVars = [key for key, val in data_dict.items() if key not in noWriteVars]

                # create another dictionary that contains only the needed variables
                data_dict_writeOut = {}
                data_dict_copy = deepcopy(data_dict)
                for key, val in data_dict_copy.items():
                    if key in writeOutVars:
                        data_dict_writeOut = {**data_dict_writeOut, **{key:val}}

                del data_dict_copy
                collect()

                # Convert everything to numpy array
                for key, val in data_dict_writeOut.items():
                    data_dict_writeOut[key][0] = np.array(data_dict_writeOut[key][0])

                # --- --- --- --- --- ---
                # --- WRITE OUT DATA ---
                # --- --- --- --- --- ---

                # --- delete output file if it already exists ---
                if os.path.exists(outputPath):
                    os.remove(outputPath)

                with pycdf.CDF(outputPath,'') as outputFile:
                    outputFile.readonly(False)

                    # --- write out global attributes ---
                    # inputGlobDic = L1ModelData.cdfFile.globalattsget()
                    inputGlobDic = {'globalAttributes':
                                    {'Source_name': f'MISSIONNAM_#####>MISSIONNAM RKT ##.###',
                                     'Data_type': 'K0>Key Parameter',
                                     'PI_name': 'EXAMPLE PI',
                                     'Logical_source': f'missionNam_#####_',
                                     'Logical_file_id': f'missionNam_#####_00000000_v01',
                                     'Logical_source_description': 'Raw Data from the MISSIONNAMII mission organized by minorframe.150 words per minor frame.40 minor frames to a major frame.',
                                     'TEXT': 'Raw Data from the MISSIONNAMII mission organized by minorframe.150 words per minor frame.40 minor frames to a major frame.'
                                     },
                                    }
                    for key, val in inputGlobDic.items():
                        if key == 'Descriptor':
                            globalAttrsMod[key] = ACESII.InstrNames_Full[wInstr[0]]

                        if key in globalAttrsMod:
                            outputFile.attrs[key] = globalAttrsMod[key]
                        else:
                            outputFile.attrs[key] = val


                    for varKey, varVal in data_dict_writeOut.items():
                        if 'Epoch' in varKey:
                            outputFile.new(varKey, data=varVal[0], type=33)
                        elif 'ni_swept' in varKey or 'ne_swept' in varKey:
                            outputFile.new(varKey, data=varVal[0], type=pycdf.const.CDF_INT4)
                        else:
                            outputFile.new(varKey, data=varVal[0])

                        # --- Write out the attributes and variable info ---
                        for attrKey, attrVal in data_dict_writeOut[varKey][1].items():
                            if attrKey == 'VALIDMIN':
                                outputFile[varKey].attrs[attrKey] = varVal[0].min()
                            elif attrKey == 'VALIDMAX':
                                outputFile[varKey].attrs[attrKey] = varVal[0].max()
                            elif attrVal != None:
                                outputFile[varKey].attrs[attrKey] = attrVal

            stl.Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 0:  # ACES II Integration High
    rocketFolderPath = DataPaths.Integration_data_folder
elif wRocket == 1: # ACES II Integration Low
    rocketFolderPath = DataPaths.Integration_data_folder
elif wRocket == 4:  # ACES II High
    rocketFolderPath = DataPaths.ACES_data_folder
elif wRocket == 5: # ACES II Low
    rocketFolderPath = DataPaths.ACES_data_folder

if len(glob(f'{rocketFolderPath}L0\{ACESII.fliers[wRocket-4]}\*.cdf')) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    if justPrintFileNames:
        L0_to_L1(wRocket, 0, rocketFolderPath, justPrintFileNames)
    elif not wFiles[wRocket-4]:
        for fileNo in (range(len(glob(f'{rocketFolderPath}L0\{ACESII.fliers[wRocket-4]}\*.cdf')))):
            L0_to_L1(wRocket, fileNo, rocketFolderPath, justPrintFileNames)
    else:
        for filesNo in wFiles[wRocket-4]:
            L0_to_L1(wRocket, filesNo, rocketFolderPath, justPrintFileNames)