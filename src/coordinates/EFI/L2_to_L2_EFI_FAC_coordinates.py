# --- L2_to_L2_EFI_FAC_coordinates.py ---
# Description: convert EFI ENU coordinates to FAC and output as L2 data



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
from scipy.interpolate import CubicSpline

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 4 -> ACES-II High Flier
# 5 -> ACES-II Low Flier
wRocket = 5
wFiles = [0]
outputData = True

Plot_correction_term = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none

def L2_to_L2_EFI_FAC_coordinates(wRocket, justPrintFileNames):

    inputFiles_elec = glob(f'{DataPaths.ACES_data_folder}\\L1\\{ACESII.fliers[wRocket-4]}\*E_Field_ENU_no_corrections*')
    inputFiles_mag = glob(f'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket-4]}\\*RingCore_ENU.cdf*')
    input_names = [ifile.replace(f'{DataPaths.ACES_data_folder}\\L1\\{ACESII.fliers[wRocket-4]}\\', '') for ifile in inputFiles_elec]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles_elec):
            print('[{:.0f}] {:80s}'.format(i, input_names[i], round(getsize(file) / (10 ** 6))))
        return

    # --- get the data from the B-Field file ---
    stl.prgMsg(f'Loading data')
    data_dict_EFI =
    data_dict_transform = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\{ACESII.fliers[wRocket - 4]}\\*ECEF_to_FAC.cdf*')[0])
    stl.Done(start_time)

    ##############################################
    # --- Correct the Rocket Convection Effect ---
    ##############################################


    # Convert payload velocity to ENU




    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')

        fileoutName = rf'ACESII_{ACESII.payload_IDs[wRocket-4]}_l2_E_Field_ENU_fullCal.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output, globalAttrsMod=GlobalAttrs, instrNam='EFI')
        stl.Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
L2_to_L2_EFI_FAC_coordinates(wRocket, justPrintFileNames)