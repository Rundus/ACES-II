# --- EFI_rkt_to_ENU_despin.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:

# TODO: need to determine which time offset is required to de-spin the EFI data

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

# --- OutputData ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.my_imports import *

#######################
# --- MAIN FUNCTION ---
#######################
def EFI_rkt_to_ENU_despun(wRocket, justPrintFileNames):

    # --- FILE I/O ---
    data_folder_path = rf'{DataPaths.ACES_data_folder}\\L1\\{ACESII.fliers[wRocket-4]}\\'

    # get the EFI files
    input_files = glob(data_folder_path + '*EFI_rktFrm*')
    input_names = [ifile.replace(data_folder_path, '') for ifile in input_files]

    # get the attitude files
    data_folder_path = rf'{DataPaths.ACES_data_folder}\\attitude\\{ACESII.fliers[wRocket-4]}\\'
    input_files_attitude = glob(data_folder_path +'*.cdf*')[0]

    if justPrintFileNames:
        for i, file in enumerate(input_files):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names[i], round(os.path.getsize(file) / (10 ** 6), 1)))
        return

    # --- get the data from the L2 file ---
    stl.prgMsg(f'Loading data from L1 EFI Files')
    data_dict_EFI = stl.loadDictFromFile(input_files[0])
    data_dict_attitude = stl.loadDictFromFile(input_files_attitude)
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}

    ############################################################
    # --- Interpolate the Attitude DCM onto the EFI timebase ---
    ############################################################
    stl.prgMsg('Interpolating attitude DCM')
    T0 = dt.datetime(2022, 11, 20, 17, 20)
    T0_EFI = np.array([pycdf.lib.datetime_to_tt2000(data_dict_EFI['Epoch'][0][i]) for i in range(len(data_dict_EFI['Epoch'][0]))])
    T0_attitude = np.array([pycdf.lib.datetime_to_tt2000(data_dict_attitude['Epoch'][0][i]) for i in range(len(data_dict_attitude['Epoch'][0]))])
    TimeOffset = [127567241.37931032, 120789473.68421052]
    T0_attitude = np.array([T0_attitude[i] for i in range(len(T0_attitude))])

    # interpolate onto EFI timebase
    from scipy.interpolate import CubicSpline
    DCM_EFI = np.zeros(shape=(len(data_dict_EFI['Epoch'][0]),3,3))
    for i in range(1, 4):
        for j in range(1, 4):
            cs = CubicSpline(T0_attitude, data_dict_attitude[f'a{i}{j}'][0])
            DCM_EFI[:, i-1, j-1] = cs(T0_EFI)

    stl.Done(start_time)

    #############################
    # --- DESPIN THE EFI DATA ---
    #############################
    stl.prgMsg('Despinning EFI data')
    EFI_rkt_vec = np.array([data_dict_EFI['E_x'][0],data_dict_EFI['E_y'][0],data_dict_EFI['E_z'][0]]).T
    EFI_ENU_vec = np.array([DCM_EFI[i]@vec for vec in EFI_rkt_vec])
    data_dict_output = {**data_dict_output,
                        **{
                            'E_E' : [EFI_ENU_vec[:, 0], {'DEPEND_0':'Epoch','UNITS': 'V/m','LABLAXIS':'E_East'}],
                            'E_N': [EFI_ENU_vec[:, 1], {'DEPEND_0':'Epoch','UNITS': 'V/m','LABLAXIS':'E_North'}],
                            'E_Up': [EFI_ENU_vec[:, 2], {'DEPEND_0':'Epoch','UNITS': 'V/m','LABLAXIS':'E_up'}],
                            '|E|': [np.array([np.linalg.norm(vec) for vec in EFI_ENU_vec]), {'DEPEND_0':'Epoch','UNITS': 'V/m','LABLAXIS':'|E|'}],
                            'Epoch' : deepcopy(data_dict_EFI['Epoch'])
                        }
                        }
    stl.Done(start_time)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        stl.prgMsg('Creating output file')
        fileoutName = f'ACESII_{ACESII.payload_IDs[wRocket-4]}_l2_EFI_ENU.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket-4]}\\'+ fileoutName
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
EFI_rkt_to_ENU_despun(wRocket, justPrintFileNames)

