# --- Interpolate_ephemeris_data_to_file.py ---
# Description: Given a list of .cdf files, interpolate the ephemeris data
# into the datafile. Always check to see if the ephemeris data is already added to the file
# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
from src.ACESII.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier

# --- OutputData ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from glob import glob
import spaceToolsLib as stl


FilePaths = [
    f'{DataPaths.ACES_data_folder}/L3/Langmuir',
    f'{DataPaths.ACES_data_folder}/L3/Energy_Flux',
    f'{DataPaths.ACES_data_folder}/L2',
    f'{DataPaths.ACES_data_folder}/L3/DERPA',
    f'{DataPaths.ACES_data_folder}/L3/RingCore',
    f'{DataPaths.ACES_data_folder}/L3/B_median',
]

interp_keys = ['L-Shell','Alt','Lat','Long']
T0_ref = dt.datetime(2022,11,20,17,20)

def Interpolate_ephemeris_data_to_file():

    # For each Rocket
    for i in range(2):
        # Load the L-Shell data
        data_path = glob(rf'{DataPaths.ACES_data_folder}/coordinates/Lshell/{ACESII.fliers[i]}/*.cdf*')[0]
        data_dict_LShell = stl.loadDictFromFile(data_path)
        T0_LShell = np.array(stl.EpochTo_T0_Rocket(data_dict_LShell['Epoch'][0],T0=T0_ref),dtype='float64')

        # Loop through each path
        for path in FilePaths:

            input_files = glob(path+f'/{ACESII.fliers[i]}/*.cdf*')

            # Loop through each file in the path
            for file in tqdm(input_files):
                data_dict = stl.loadDictFromFile(file)

                # check if Emphemeris is already in this file, if so skip this file
                # if not all(i in data_dict.keys() for i in interp_keys):
                try:
                    # Try to add the ephemeris
                    T0_data = np.array(stl.EpochTo_T0_Rocket(data_dict['Epoch'][0],T0=T0_ref),dtype='float64')

                    for iKey in interp_keys:
                        if iKey == 'Alt':
                            data_dict = {**data_dict,
                                         **{iKey: [np.interp(T0_data, T0_LShell, data_dict_LShell[iKey][0]), deepcopy(data_dict_LShell[iKey][1])]}
                                         }
                            data_dict[iKey][0] = data_dict[iKey][0]/1000
                            data_dict[iKey][1]['UNITS'] = 'km'
                        else:
                            data_dict = {**data_dict,
                                         **{iKey:[np.interp(T0_data,T0_LShell,data_dict_LShell[iKey][0]), deepcopy(data_dict_LShell[iKey][1])]}
                                         }

                    if outputData:
                        stl.outputDataDict(
                            outputPath=file,
                            data_dict=data_dict
                        )
                except Exception as err:
                    print(stl.color.RED + f'\nFAILED: '+stl.color.END + f'{file}')  # check
                    print('Reason ', err, '\n')
                    continue






Interpolate_ephemeris_data_to_file()