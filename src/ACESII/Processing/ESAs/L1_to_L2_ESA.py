# --- L1_to_L2_ESA.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert electrostatic analyzer data from counts to differential energy flux

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
from src.ACESII.my_imports import *
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
wFiles = [[0, 2, 4], [0, 3]]

inputPath_modifier = 'L1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'L2' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# output the data
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import itertools
from warnings import filterwarnings # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.
filterwarnings("ignore")

def L1_to_L2_ESA(wRocket, wFile, rocketFolderPath, justPrintFileNames):

    # --- ACES II Flight/Integration Data ---
    L1Files = glob(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}\*.cdf')
    L2Files = glob(f'{rocketFolderPath}{outputPath_modifier}\{ACESII.fliers[wRocket-4]}\*.cdf')

    L1_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}\\', '') for ifile in L1Files]
    L2_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{ACESII.fliers[wRocket-4]}\\', '') for ofile in L2Files]

    L1_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace('l1_', '').replace('_v00', '') for ifile in L1_names]
    L2_names_searchable = [ofile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace('l2_', '').replace('_v00', '').replace('__', '_') for ofile in L2_names]

    dataFile_name = L1Files[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wRocket-4]}\\', '')
    fileoutName = dataFile_name.replace('l1', 'l2')

    if justPrintFileNames:
        for i, file in enumerate(L1Files):
            anws = ["yes" if L1_names_searchable[i].replace('.cdf', "") in L2_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB   Made L2: {:3s} '.format(i, L1_names_searchable[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
        return
    print('\n')
    print(stl.color.UNDERLINE + f'Converting to L2 data for {dataFile_name}' + stl.color.END)
    print('[' + str(wFile) + ']   ' + str(round(getsize(L1Files[wFile]) / (10 ** 6), 1)) + 'MiB')

    ##########################
    # --- LOAD THE L1 DATA ---
    ##########################
    stl.prgMsg('Loading data from L1Files')
    data_dict_L1 = stl.loadDictFromFile(L1Files[wFile])
    stl.Done(start_time)

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in dataFile_name:
            wInstr = [index, instr]

    # --- prepare the output ---
    ESA_dimensions = (len(data_dict_L1['counts'][0]),len(data_dict_L1['Pitch_Angle'][0]), len(data_dict_L1['Energy'][0]))
    data_dict_output = {'Differential_Number_Flux': [np.zeros(shape=ESA_dimensions), {'LABLAXIS': 'Differential_Number_Flux',
                                                  'DEPEND_0': 'Epoch', 'DEPEND_1': 'Pitch_Angle',
                                                  'DEPEND_2': 'Energy',
                                                  'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2',
                                                  'UNITS': 'cm!U-2!N str!U-1!N s!U-1!N eV!U-1!N',
                                                  'VAR_TYPE': 'data', 'SCALETYP': 'log'}],
                        'Differential_Energy_Flux': [np.zeros(shape=ESA_dimensions), {'LABLAXIS': 'Differential_Energy_Flux',
                                         'DEPEND_0': 'Epoch', 'DEPEND_1': 'Pitch_Angle',
                                         'DEPEND_2': 'Energy',
                                         'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2',
                                         'UNITS': 'cm!U-2!N str!U-1!N s!U-1!N eV/eV',
                                         'VAR_TYPE': 'data', 'SCALETYP': 'log'}],
                        'Differential_Number_Flux_stdDev': [np.zeros(shape=ESA_dimensions), {'LABLAXIS': 'Differential_Energy_Flux_stdDev',
                                                'DEPEND_0': 'Epoch', 'DEPEND_1': 'Pitch_Angle',
                                                'DEPEND_2': 'Energy',
                                                'FILLVAL': ACESII.epoch_fillVal,
                                                'FORMAT': 'E12.2',
                                                'UNITS': 'cm!U-2!N str!U-1!N s!U-1!N eV/eV',
                                                'VAR_TYPE': 'data', 'SCALETYP': 'log'}]
                        }

    ###################################
    # --- Calculate Instrument Data ---
    ###################################
    stl.prgMsg('Creating L2 instrument data')

    Energies = data_dict_L1['Energy'][0]
    counts = data_dict_L1['counts'][0]
    geo_factor = ACESII.ESA_geometric_factor_TRACERS_ACE[wInstr[0]]
    count_interval = (1E-3)*data_dict_L1['Count_Interval'][0] # amount of time BETWEEN energy steps that we actually collect counts (in seconds).

    # --- PROCESS ESA DATA ---
    for tme, ptch, engy in tqdm(itertools.product(*[range(ESA_dimensions[0]),range(ESA_dimensions[1]),range(ESA_dimensions[2])])):
        deltaT = (count_interval[tme]) - (counts[tme][ptch][engy] * ACESII.ESA_deadtime[wRocket-4])

        if counts[tme][ptch][engy] >= 0:
            try:
                data_dict_output['Differential_Number_Flux'][0][tme][ptch][engy] = int((counts[tme][ptch][engy]) / (Energies[engy] * geo_factor[ptch] * deltaT))
                data_dict_output['Differential_Number_Flux_stdDev'][0][tme][ptch][engy] = int(np.sqrt(counts[tme][ptch][engy]) / (Energies[engy] * geo_factor[ptch] * deltaT))
                data_dict_output['Differential_Energy_Flux'][0][tme][ptch][engy] = int((counts[tme][ptch][engy]) / (geo_factor[ptch] * deltaT))
            except:
                print(counts[tme][ptch][engy], counts[tme][ptch][engy]>=0)
        else:
            data_dict_output['Differential_Number_Flux'][0][tme][ptch][engy] = ACESII.epoch_fillVal
            data_dict_output['Differential_Number_Flux_stdDev'][0][tme][ptch][engy] = ACESII.epoch_fillVal
            data_dict_output['Differential_Energy_Flux'][0][tme][ptch][engy] = ACESII.epoch_fillVal

    stl.Done(start_time)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    stl.prgMsg('Creating output file')
    if outputData:
        del data_dict_L1['counts']
        data_dict_output = {**data_dict_output,
                            **data_dict_L1}


        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath= outputPath,data_dict= data_dict_output,instrNam=wInstr[1])
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

if len(glob(f'{rocketFolderPath}L1\{ACESII.fliers[wRocket-4]}\*.cdf')) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    if justPrintFileNames:
        L1_to_L2_ESA(wRocket, 0, rocketFolderPath, justPrintFileNames)
    elif not wFiles[wRocket-4]:
        for fileNo in (range(len(glob(f'{rocketFolderPath}L1\{ACESII.fliers[wRocket-4]}\*.cdf')))):
            L1_to_L2_ESA(wRocket, fileNo, rocketFolderPath, justPrintFileNames)
    else:
        for filesNo in wFiles[wRocket-4]:
            L1_to_L2_ESA(wRocket, filesNo, rocketFolderPath, justPrintFileNames)