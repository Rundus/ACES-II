# --- RingCore_l2_median_filter ---
# Description: Median Filter the RingCore data that DOESN'T have the B-Field model in it
# in order to get the deltaB corresponding to currents. Include the L-Shell data


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
fs = 128 # sample frequency of the data

low_cutoff = 1/40 # 20 seconds (or 1/20 freq) butterworth cutoff
order = 4
filter_window = 30*fs+1 # corresponds to 30 seconds


# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4

# --- OutputData ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.my_imports import *
from src.mission_attributes import ACESII

#######################
# --- MAIN FUNCTION ---
#######################

def RingCore_l2_median_filter(wRocket):

    B_files = f'{DataPaths.ACES_data_folder}\\l2\\{ACESII.fliers[wRocket - 4]}\*RingCore_Field_Aligned.cdf*'
    L_files = f'{DataPaths.ACES_data_folder}\\coordinates\\Lshell\\{ACESII.fliers[wRocket - 4]}\*.cdf*'

    # --- FILE I/O ---
    stl.prgMsg(f'Loading data from L2 Files')
    data_dict_mag = stl.loadDictFromFile(glob(B_files)[0])
    data_dict_Lshell = stl.loadDictFromFile(glob(L_files)[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}


    ###################################
    # --- Perform the Median Filter ---
    ###################################
    stl.prgMsg('Median Filtering Data')
    from scipy.signal import medfilt
    dat_keys = ['r', 'e', 'p']

    # Filter the data
    for key in dat_keys:

        # First butterworth filter the data
        filtered = stl.butterFilter().butter_filter(data=deepcopy(data_dict_mag[f'B_{key}'][0]),
                                                    lowcutoff=low_cutoff,
                                                    highcutoff=low_cutoff,
                                                    order=order,
                                                    fs=fs,
                                                    filtertype='lowPass'
                                                    )

        # Median filter the data

        data_dict_mag[f'B_{key}'][0]= medfilt(volume=filtered, kernel_size=filter_window)
        # data_dict_mag[f'B_{key}'][0] = filtered


    data_dict_output = {**data_dict_output,
                        **data_dict_mag}
    stl.Done(start_time)


    ##########################################
    # --- Interpolate L-Shell into dataset ---
    ##########################################
    from scipy.interpolate import CubicSpline
    Epochtt2000_mag = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_mag['Epoch'][0]])
    Epochtt2000_Lshell = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_Lshell['Epoch'][0]])
    LShell_interp = np.zeros(shape=np.shape(data_dict_mag['Epoch'][0]))
    cs = CubicSpline(Epochtt2000_Lshell,data_dict_Lshell['L-Shell'][0])
    LShell_interp = np.array([cs(val)for val in Epochtt2000_mag])

    data_dict_output = {**data_dict_output,
                      **{'L-Shell':deepcopy(data_dict_Lshell['L-Shell'])}
                      }
    data_dict_output['L-Shell'][0] = LShell_interp

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:
        stl.prgMsg('Creating output file')
        file_out_name = f'ACESII_{ACESII.payload_IDs[wRocket-4]}_RingCore_Field_Aligned_median_filter.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\\l2\\{ACESII.fliers[wRocket - 4]}\{file_out_name}'
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

RingCore_l2_median_filter(wRocket)

