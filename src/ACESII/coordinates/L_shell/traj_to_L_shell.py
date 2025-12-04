# --- traj_to_L_shell.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: calculate the dipolar L-Shell traversed by the rocket throughout its journey



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
wRocket = [4, 5]
inputPath_modifier = 'attitude' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = '\coordinates\Lshell' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import spaceToolsLib as stl
stl.setupPYCDF()
from glob import glob
from spacepy import coordinates as coord
coord.DEFAULTS.set_values(use_irbem=False, itol=5)  # maximum separation, in seconds, for which the coordinate transformations will not be recalculated. To force all transformations to use an exact transform for the time, set ``itol`` to zero.
from spacepy.time import Ticktock #used to determine the time I'm choosing the reference geomagentic field
from scipy.interpolate import CubicSpline


def traj_to_L_shell(wRocket):

    # --- ACES II Flight/Integration Data ---

    rocketID = ACESII.payload_IDs[wRocket-4]

    inputFiles = glob(f'{DataPaths.ACES_data_folder}\\trajectories\\{DataPaths.fliers[wRocket-4]}\*_ECEF.cdf*')
    input_names = [ifile.replace(f'{DataPaths.ACES_data_folder}{inputPath_modifier}\{DataPaths.fliers[wRocket-4]}\\', '') for ifile in inputFiles]
    input_names_searchable = [ifile.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    fileoutName = f'ACESII_{rocketID}_Lshell.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return

    # --- get the data from the attitude file ---
    stl.prgMsg(f'Loading data from {inputPath_modifier} Files')
    data_dict_attitude = stl.loadDictFromFile(inputFiles[0])
    data_dict_attitude['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_attitude['Epoch'][0]])
    Epoch = data_dict_attitude['Epoch'][0]
    stl.Done(start_time)

    ############################################
    # --- CONVERT TO GEOMAGNETIC COORDINATES ---
    ############################################
    stl.prgMsg('Converting to geomagnetic coordinates')

    # Trajectory Data
    geodeticPos = np.array([[data_dict_attitude['Alt'][0][i], data_dict_attitude['Lat'][0][i], data_dict_attitude['Long'][0][i]] for i in range(len(Epoch))])
    ISOtime = [pycdf.lib.tt2000_to_datetime(tme).isoformat() for tme in Epoch]
    cvals_GDZ = coord.Coords(geodeticPos, 'GDZ', 'sph')
    cvals_GDZ.ticks = Ticktock(ISOtime, 'ISO')
    cvals_GDZ_MAG = cvals_GDZ.convert('MAG', 'sph')

    geomagAlt = np.array([cvals_GDZ_MAG[i].radi for i in range(len(Epoch))])
    geomagLat = np.array([cvals_GDZ_MAG[i].lati for i in range(len(Epoch))])
    geomagLong = np.array([cvals_GDZ_MAG[i].long for i in range(len(Epoch))])
    stl.Done(start_time)

    # ###################################################
    # --- Calculate L-Shell Parameter over the flight ---
    # ###################################################
    stl.prgMsg('Calculating L-Shell')
    L_shell = np.array([geomagAlt[i]/((np.cos(np.radians(geomagLat[i])))**2) for i in range(len(data_dict_attitude['Epoch'][0]))])
    stl.Done(start_time)

    # Interpolate onto the EEPAA timeframe
    stl.prgMsg('Interpolating on EEPAA timebase')
    data_dict_eepaa = stl.loadDictFromFile(f'C:\Data\ACESII\L2\{DataPaths.fliers[wRocket-4]}\ACESII_{rocketID}_l2_eepaa_fullCal.cdf')
    Epoch_eepaa_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_eepaa['Epoch'][0]])
    Epoch_attitude_tt2000 = data_dict_attitude['Epoch'][0]
    cs = CubicSpline(x=Epoch_attitude_tt2000, y=L_shell)
    L_shell = np.array([cs(val)[0] for val in Epoch_eepaa_tt2000])
    stl.Done(start_time)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('creating output file')

        varAttrs = {'LABLAXIS': 'L-shell','DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None,
                                                   'FILLVAL': ACESII.epoch_fillVal, 'FORMAT': 'E12.2',
                                                   'UNITS': None,
                                                   'VALIDMIN': L_shell.min(), 'VALIDMAX': L_shell.max(),
                                                   'VAR_TYPE': 'data', 'SCALETYP': 'linear'}
        Epoch_output = deepcopy(data_dict_attitude['Epoch'])
        Epoch_output[1]['VAR_TYPE'] = 'support_data'


        data_dict_output = {'Epoch':deepcopy(data_dict_eepaa['Epoch']),
                            'L-Shell':[L_shell,varAttrs],
                            'diffNFlux':data_dict_eepaa['Differential_Number_Flux'],
                            'Energy': data_dict_eepaa['Energy'],
                            'Pitch_Angle': data_dict_eepaa['Pitch_Angle'],
                            'Alt': data_dict_eepaa['Alt'],
                            'Long': data_dict_eepaa['Long'],
                            'Lat': data_dict_eepaa['Lat'],
                            }

        data_dict_output['diffNFlux'][1]['DEPEND_0'] = 'L-Shell'

        for key in data_dict_output.keys():
            data_dict_output[key][0] = np.array(data_dict_output[key][0])

            # adjust the L-Shell for the Low Flyer: we'd like to keep the same timebase as the High Flyer BUT the interpolation
            # goes crazy if you evaluate it too far. Just set these points to the beginning of the Low Flyer's L-Shell
            if wRocket == 5:
                target_idx = np.abs(data_dict_output['Epoch'][0] - dt.datetime(2022,11,20,17,21,40)).argmin()
                target_L_Shell = deepcopy(data_dict_output['L-Shell'][0][target_idx])
                data_dict_output['L-Shell'][0][:target_idx] = target_L_Shell

        outputPath = f'{DataPaths.ACES_data_folder}{outputPath_modifier}\{DataPaths.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output)

        stl.Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
for val in wRocket:
    traj_to_L_shell(val)
