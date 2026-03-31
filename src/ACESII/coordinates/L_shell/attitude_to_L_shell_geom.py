# --- attitude_to_L_shell_geom.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: calculate the dipolar L-Shell traversed by the rocket throughout its journey

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

from src.ACESII.data_tools.my_imports import *
start_time = time.time()
# --- --- --- --- ---



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

just_print_file_names_bool = False
rocket_str = 'low'
wInstr = 'attitude'
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['/', [[0],[0]]]
}
outputData = True

#################
# --- IMPORTS ---
#################
from src.ACESII.data_tools.my_imports import *
import spaceToolsLib as stl
stl.setupPYCDF()
from glob import glob
from spacepy import coordinates as coord
coord.DEFAULTS.set_values(use_irbem=False, itol=5)  # maximum separation, in seconds, for which the coordinate transformations will not be recalculated. To force all transformations to use an exact transform for the time, set ``itol`` to zero.
from spacepy.time import Ticktock #used to determine the time I'm choosing the reference geomagentic field
from scipy.interpolate import CubicSpline


def traj_to_L_shell(data_dicts):

    # --- get the data from the attitude file ---
    stl.prgMsg(f'Loading data')
    data_dict_attitude = deepcopy(data_dicts[0])
    Epoch = data_dict_attitude['Epoch'][0]
    stl.Done(start_time)

    ############################################
    # --- CONVERT TO GEOMAGNETIC COORDINATES ---
    ############################################
    stl.prgMsg('Converting to geomagnetic coordinates')

    # Trajectory Data
    pos = np.array([data_dict_attitude['Alt'][0]/stl.m_to_km, data_dict_attitude['Lat'][0], data_dict_attitude['Long'][0] ]).T
    ISOtime = [tme.isoformat() for tme in Epoch]
    cvals_GDZ = coord.Coords(pos, 'GDZ', 'sph')
    cvals_GDZ.ticks = Ticktock(ISOtime, 'ISO')
    cvals_GDZ_MAG = cvals_GDZ.convert('MAG', 'sph')

    mAlt = np.array([cvals_GDZ_MAG[i].radi[0] for i in range(len(Epoch))])
    mLat = np.array([cvals_GDZ_MAG[i].lati[0] for i in range(len(Epoch))])
    mLong = np.array([cvals_GDZ_MAG[i].long[0] for i in range(len(Epoch))])
    stl.Done(start_time)

    # ###################################################
    # --- Calculate L-Shell Parameter over the flight ---
    # ###################################################
    stl.prgMsg('Calculating L-Shell')
    L_shell = mAlt/np.square(np.cos(np.radians(mLat)))
    stl.Done(start_time)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('creating output file')

        data_dict_output = {'Epoch':deepcopy(data_dict_attitude['Epoch']),
                            'L-Shell':[L_shell,{'LABLAXIS': 'L-shell','DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None,
                                                   'FILLVAL': ACESII.epoch_fillVal,
                                                   'UNITS': None,
                                                   'VALIDMIN': L_shell.min(), 'VALIDMAX': L_shell.max(),
                                                   'VAR_TYPE': 'data', 'SCALETYP': 'linear'}],
                            'Alt': deepcopy(data_dict_attitude['Alt']),
                            'Long': deepcopy(data_dict_attitude['Long']),
                            'Lat': deepcopy(data_dict_attitude['Lat']),
                            'mLat': [mLat, {'DEPEND_0':'Epoch','UNITS':'Degrees','LABLAXIS':'mLat'}] ,
                            'mLong': [mLong, {'DEPEND_0':'Epoch','UNITS':'Degrees','LABLAXIS':'mLong'}],
                            }

        if outputData:
            file_name = f'ACESII_{ACESII.fliers_dict[rocket_str]}_LShell.cdf'
            outputPath = f'{DataPaths.ACES_data_folder}/coordinates/Lshell/{rocket_str}/' + file_name
            stl.outputDataDict(outputPath, data_dict_output)

        stl.Done(start_time)

#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(traj_to_L_shell,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)
