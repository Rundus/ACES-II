# --- Interpolate_ephemeris_data_to_file.py ---
# Description: Given a list of .cdf files, interpolate the ephemeris data
# into the datafile. Always check to see if the ephemeris data is already added to the file
# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
import time
start_time = time.time()
# --- --- --- --- ---

# --- OutputData ---
outputData = True
override_checker = True # forces a re-interpolation of the datafiles

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import datetime as dt
from src.ACESII.data_tools.my_imports import *


FilePaths = [
    f'{DataPaths.ACES_data_folder}/L3/',
    f'{DataPaths.ACES_data_folder}/L2/',
]

Instrs = [
    'EEPAA',
    'EFI',
    'ERPA',
    'IEPAA',
    'LEESA',
    'LP',
    'MAG'
]

interp_keys = ['L-Shell','Alt','Lat','Long', 'mLat','mLong']
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

            for instrNam in Instrs:

                input_files = glob(path+f'/{instrNam}/'+f'/{ACESII.fliers[i]}/*.cdf*')

                # Loop through each file in the path
                for file in tqdm(input_files):
                    data_dict = stl.loadDictFromFile(file)

                    # check if Emphemeris is already in this file, if so skip this file
                    if not all(key in data_dict.keys() for key in interp_keys):
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

                    # if you wish to re-set the ephemeris data even if it's already present in the file
                    elif override_checker:
                        print('OVERRIDE')
                        try:
                            # Try to add the ephemeris
                            T0_data = np.array(stl.EpochTo_T0_Rocket(data_dict['Epoch'][0], T0=T0_ref), dtype='float64')

                            for iKey in interp_keys:
                                if iKey == 'Alt':
                                    data_dict[iKey][0] = np.interp(T0_data, T0_LShell, data_dict_LShell[iKey][0])
                                    data_dict[iKey][0] = data_dict[iKey][0] / 1000
                                    data_dict[iKey][1]['UNITS'] = 'km'
                                else:
                                    data_dict[iKey][0] = np.interp(T0_data, T0_LShell, data_dict_LShell[iKey][0])

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