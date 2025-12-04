# --- convection_E_Field.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: use the E = - vxB formula to determine the convection velocity of the E-Field
# in the ACES-II mission.



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

from src.ACESII.my_imports import *
from scipy.interpolate import CubicSpline
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# ---------------------------
outputData = True
# ---------------------------

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none


def convection_E_Field():

    # --- --- --- --- --- -
    # --- LOAD THE DATA ---
    # --- --- --- --- --- -

    # --- get the data from the file ---
    stl.prgMsg(f'Loading data')
    data_dict_mag = stl.loadDictFromFile(r'C:\Data\ACESII\L2\low\ACESII_36364_l2_RingCore_Field_Aligned.cdf')
    data_dict_EFI = stl.loadDictFromFile(r'C:\Data\ACESII\L2\low\ACESII_36364_l2_E_Field_FAC_fullCal.cdf')
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}
    for key, val in data_dict_EFI.items():
        if key not in ['Alt_geom','Lat_geom','Long_geom']:
            data_dict_output = {**data_dict_output,
                                f'{key}':val}

    # --- --- --- --- --- --- --- --- --- ---
    # --- CALCULATE Convection Velocities ---
    # --- --- --- --- --- --- --- --- --- ---

    # [1] Interpolate the magnetic field strength onto EFI timebase
    T0 = pycdf.lib.datetime_to_tt2000(dt.datetime(2022,11,20,17,20))/1E9
    Epoch_seconds_MAG = np.array([pycdf.lib.datetime_to_tt2000(val)/1E9 - T0 for val in data_dict_mag['Epoch'][0]])
    Epoch_seconds_EFI = np.array([pycdf.lib.datetime_to_tt2000(val)/1E9 - T0 for val in data_dict_EFI['Epoch'][0]])
    cs = CubicSpline(Epoch_seconds_MAG, data_dict_mag['Bmag'][0])
    Bmag_EFI_timebase = (1E-9)*cs(Epoch_seconds_EFI)

    # [2] Calculate the convection velocities
    Vr = data_dict_EFI['E_e'][0] /Bmag_EFI_timebase
    Ve = -1*data_dict_EFI['E_r'][0] / Bmag_EFI_timebase

    data_dict_output = {
        **data_dict_output,
        **{'Vr':[Vr,{'DEPEND_0': 'Epoch', 'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2',
                        'UNITS': 'm/s', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data', 'SCALETYP': 'linear',
                        'LABLAXIS': 'Vr convection Velocity'} ],
           'Ve':[Ve,{'DEPEND_0': 'Epoch', 'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2',
                        'UNITS': 'm/s', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data', 'SCALETYP': 'linear',
                        'LABLAXIS': 'Ve convection Velocity'}]
           }
        }

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')
        fileoutName = f'ACESII_{ACESII.payload_IDs[1]}_E_Field_convection_velocity.cdf'
        outputPath = rf'C:\Data\ACESII\science\convection_velocities_E_Field\low\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)

# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = DataPaths.ACES_data_folder
convection_E_Field()
