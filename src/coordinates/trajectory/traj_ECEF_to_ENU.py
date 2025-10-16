# --- traj_ECEF_to_ENU.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: input the payload trajectory files and convert the ECEF velocities, positions into ENU coordinates (E,N,U)
# as well as determine the gradient in the distances

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

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = [4, 5]
coordKeys = ['E','N','U']

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import pyproj

def traj_ECEF_to_ENU(wRocket):

    # --- Load the trajectory Data ---
    inputFiles = glob(f'{DataPaths.ACES_data_folder}\\trajectories\{ACESII.fliers[wRocket-4]}\\*GPS_trajectory_ECEF*')[0]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles[i].replace(f'{DataPaths.ACES_data_folder}\{ACESII.fliers[wRocket-4]}',''), round(getsize(file) / (10 ** 6), 1)))
        return

    # --- get the data from the tmCDF file ---
    stl.prgMsg(f'Loading data')
    data_dict_traj = stl.loadDictFromFile(inputFiles)
    data_dict_transform = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\transforms\\{ACESII.fliers[wRocket-4]}\*ECEF_to_{"".join(coordKeys)}*')[0])
    data_dict_LShell = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\Lshell\\{ACESII.fliers[wRocket-4]}\*Lshell.cdf*')[0])
    stl.Done(start_time)

    # --- prepare output ---
    data_dict_output = {
        f'{coordKeys[0]}_POS': [np.zeros(shape=np.shape(data_dict_traj['ECEFXPOS'][0])), deepcopy(data_dict_traj['ECEFXPOS'][1])],
        f'{coordKeys[1]}_POS': [np.zeros(shape=np.shape(data_dict_traj['ECEFYPOS'][0])), deepcopy(data_dict_traj['ECEFYPOS'][1])],
        f'{coordKeys[2]}_POS': [np.zeros(shape=np.shape(data_dict_traj['ECEFZPOS'][0])), deepcopy(data_dict_traj['ECEFZPOS'][1])],
        f'{coordKeys[0]}_POS_GRAD': [np.zeros(shape=np.shape(data_dict_traj['ECEFZPOS'][0])), deepcopy(data_dict_traj['ECEFZPOS'][1])],
        f'{coordKeys[1]}_POS_GRAD': [np.zeros(shape=np.shape(data_dict_traj['ECEFZPOS'][0])), deepcopy(data_dict_traj['ECEFZPOS'][1])],
        f'{coordKeys[2]}_POS_GRAD': [np.zeros(shape=np.shape(data_dict_traj['ECEFZPOS'][0])), deepcopy(data_dict_traj['ECEFZPOS'][1])],
        f'{coordKeys[0]}_VEL': [np.zeros(shape=np.shape(data_dict_traj['ECEFXVEL'][0])), deepcopy(data_dict_traj['ECEFXVEL'][1])],
        f'{coordKeys[1]}_VEL': [np.zeros(shape=np.shape(data_dict_traj['ECEFYVEL'][0])), deepcopy(data_dict_traj['ECEFYVEL'][1])],
        f'{coordKeys[2]}_VEL': [np.zeros(shape=np.shape(data_dict_traj['ECEFZVEL'][0])), deepcopy(data_dict_traj['ECEFZVEL'][1])],
        'Epoch': deepcopy(data_dict_traj['Epoch']),
        'Alt': deepcopy(data_dict_traj['Alt']),
        'Lat': deepcopy(data_dict_traj['Lat']),
        'Long': deepcopy(data_dict_traj['Long']),
    }

    #############################################
    # --- CONVERT DATA TO AURORAL COORDINATES ---
    #############################################
    stl.prgMsg('Converting Coordinates')

    # form the rocket ECEF position vector
    rkt_pos_ECEF = np.array([data_dict_traj['ECEFXPOS'][0],data_dict_traj['ECEFYPOS'][0],data_dict_traj['ECEFZPOS'][0]]).T

    # for the position vector for the Launch site
    def gps_to_ecef_pyproj(lat, lon, alt):
        ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
        lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
        x, y, z = pyproj.transform(lla, ecef, lon, lat, alt, radians=False)

        return x, y, z

    # DEFINE THE X0 POINT for the position vector
    Lat_initial = data_dict_LShell['Lat'][0][0]
    Long_initial = data_dict_LShell['Long'][0][0]
    Alt_initial = data_dict_LShell['Alt'][0][0]

    reference_position_ECEF = np.array(gps_to_ecef_pyproj(lat=Lat_initial,lon=Long_initial, alt=Alt_initial))
    position_vector_ECEF = (rkt_pos_ECEF- reference_position_ECEF)

    # determine the in situ position vector and rotate it into auroral coordinates
    rkt_pos_auroral = np.array([np.matmul(data_dict_transform[f'ECEF_to_{"".join(coordKeys)}'][0][i], vec) for i,vec in enumerate(position_vector_ECEF)])
    data_dict_output[f'{coordKeys[0]}_POS'][0] = rkt_pos_auroral[:, 0]
    data_dict_output[f'{coordKeys[1]}_POS'][0] = rkt_pos_auroral[:, 1]
    data_dict_output[f'{coordKeys[2]}_POS'][0] = rkt_pos_auroral[:, 2]

    # update the attributes
    data_dict_output[f'{coordKeys[0]}_POS'][1]['LABLAXIS'] = 'East Position from X0'
    data_dict_output[f'{coordKeys[0]}_POS'][1]['UNITS'] = 'm'
    data_dict_output[f'{coordKeys[1]}_POS'][1]['LABLAXIS'] = 'North Position from X0'
    data_dict_output[f'{coordKeys[1]}_POS'][1]['UNITS'] = 'm'
    data_dict_output[f'{coordKeys[2]}_POS'][1]['LABLAXIS'] = 'Up Position from X0'
    data_dict_output[f'{coordKeys[2]}_POS'][1]['UNITS'] = 'm'

    # form the ECEF velocity vector
    rkt_vel_ECEF = np.array([data_dict_traj['ECEFXVEL'][0], data_dict_traj['ECEFYVEL'][0], data_dict_traj['ECEFZVEL'][0]]).T
    rkt_vel_ENU = np.array([np.matmul(data_dict_transform[f'ECEF_to_{"".join(coordKeys)}'][0][i], vec) for i,vec in enumerate(rkt_vel_ECEF)])
    data_dict_output[f'{coordKeys[0]}_VEL'][0] = rkt_vel_ENU[:, 0]
    data_dict_output[f'{coordKeys[1]}_VEL'][0] = rkt_vel_ENU[:, 1]
    data_dict_output[f'{coordKeys[2]}_VEL'][0] = rkt_vel_ENU[:, 2]

    # update the attributes
    data_dict_output[f'{coordKeys[0]}_VEL'][1]['LABLAXIS'] = 'East Velocity'
    data_dict_output[f'{coordKeys[1]}_VEL'][1]['LABLAXIS'] = 'North Velocity'
    data_dict_output[f'{coordKeys[2]}_VEL'][1]['LABLAXIS'] = 'Up Velocity'
    stl.Done(start_time)

    #######################################################
    # --- Calculate the Position Gradient from Velocity ---
    #######################################################

    # --- get the Epoch in seconds from Launch ---
    Launch_time = pycdf.lib.tt2000_to_datetime(ACESII.launch_T0_TT2000[wRocket-4])
    seconds_from_T0 = stl.EpochTo_T0_Rocket(InputEpoch=data_dict_output['Epoch'][0], T0=Launch_time)
    deltaT = np.array([0] + [seconds_from_T0[i+1] - seconds_from_T0[i] for i in range(len(seconds_from_T0)-1)])
    data_dict_output[f'{coordKeys[0]}_POS_GRAD'][0] = deepcopy(deltaT)*deepcopy(data_dict_output[f'{coordKeys[0]}_VEL'][0])
    data_dict_output[f'{coordKeys[0]}_POS_GRAD'][1]['UNITS'] = 'm'
    data_dict_output[f'{coordKeys[0]}_POS_GRAD'][1]['LABLAXIS'] = 'Tangent Position Gradient'

    data_dict_output[f'{coordKeys[1]}_POS_GRAD'][0] = deepcopy(deltaT) * deepcopy(data_dict_output[f'{coordKeys[1]}_VEL'][0])
    data_dict_output[f'{coordKeys[1]}_POS_GRAD'][1]['UNITS'] = 'm'
    data_dict_output[f'{coordKeys[1]}_POS_GRAD'][1]['LABLAXIS'] = 'Normal Position Gradient'

    data_dict_output[f'{coordKeys[2]}_POS_GRAD'][0] = deepcopy(deltaT) * deepcopy(data_dict_output[f'{coordKeys[2]}_VEL'][0])
    data_dict_output[f'{coordKeys[2]}_POS_GRAD'][1]['UNITS'] = 'm'
    data_dict_output[f'{coordKeys[2]}_POS_GRAD'][1]['LABLAXIS'] = 'Field-Aligned Position Gradient'


    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')

        fileoutName = rf'ACESII_{ACESII.payload_IDs[wRocket-4]}_GPS_trajectory_{"".join(coordKeys)}.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\\trajectories\\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)
        print('\n')




# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
for val in wRocket:
    traj_ECEF_to_ENU(val)