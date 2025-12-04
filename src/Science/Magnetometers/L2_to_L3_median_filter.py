# --- L2_to_L3_median_filter.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: median filter the ACES-II RingCore data

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
from src.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False # Just print the names of files

# --- Filter Toggles ---
fs = 128 # sample frequency of the data
low_cutoff = 1/40 # 20 seconds (or 1/20 freq) butterworth cutoff
order = 4
filter_window = 30*fs+1 # corresponds to 30 seconds

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

# --- OutputData ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.my_imports import *
from scipy.signal import savgol_filter

#######################
# --- MAIN FUNCTION ---
#######################
def L2_to_L3_median_filter(wRocket, justPrintFileNames):

    # Set the paths for the file names
    data_folder = rf'{DataPaths.ACES_data_folder}/L2/{ACESII.fliers[wRocket-4]}/'
    input_files = glob(data_folder+'*RingCore_auroral.cdf*')
    input_names = [ifile.replace(data_folder, '') for ifile in input_files]

    if justPrintFileNames:
        for i, file in enumerate(input_files):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names[i], round(os.path.getsize(file) / (10 ** 6), 1)))
        return

    # --- get the data from the L2 file ---
    stl.prgMsg(f'Loading data Files')
    data_dict = stl.loadDictFromFile(input_files[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}

    ############################
    # --- MEDIAN FILTER DATA ---
    ############################
    stl.prgMsg('Median Filtering Data')
    from scipy.signal import medfilt
    dat_keys = ['N', 'T', 'p']

    # Filter the data
    for key in dat_keys:
        # First butterworth filter the data
        filtered = stl.butterFilter().butter_filter(data=deepcopy(data_dict[f'B_{key}'][0]),
                                                    lowcutoff=low_cutoff,
                                                    highcutoff=low_cutoff,
                                                    order=order,
                                                    fs=fs,
                                                    filtertype='lowPass'
                                                    )

        # Median filter the data

        med_filtered = medfilt(volume=filtered, kernel_size=filter_window)
        # data_dict_mag[f'B_{key}'][0] = filtered

        # savitz-golay filter the data
        data_dict[f'B_{key}'][0] = savgol_filter(x=med_filtered,
                                                     window_length=filter_window,
                                                     polyorder=3)

        # rename the key
        data_dict[f'B_{key}_median'] = data_dict.pop(f'B_{key}')

    data_dict_output = {**data_dict_output,
                        **data_dict}
    stl.Done(start_time)


    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        stl.prgMsg('Creating output file')
        fileoutName_fixed = f'ACESII_{ACESII.payload_IDs[wRocket-4]}_L3_RingCore_auroral.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/L3/B_median/{ACESII.fliers[wRocket-4]}/{fileoutName_fixed}'
        stl.outputDataDict(outputPath, data_dict_output)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
L2_to_L3_median_filter(wRocket, justPrintFileNames)

