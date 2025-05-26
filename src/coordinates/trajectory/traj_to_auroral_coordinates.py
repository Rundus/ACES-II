# --- traj_to_auroral_coordinates.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: input the payload trajectory files and convert the ECEF velocities, positions into auroral coordinates (T,N,p)
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

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import pyproj

def traj_to_auroral_coordinates(wRocket):

    # --- Load the trajectory Data ---
    inputFiles = glob(f'{DataPaths.ACES_data_folder}\\trajectories\{ACESII.fliers[wRocket-4]}\\*GPS_trajectory_ECEF*')[0]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles[i].replace(f'{DataPaths.ACES_data_folder}\{ACESII.fliers[wRocket-4]}',''), round(getsize(file) / (10 ** 6), 1)))
        return


    # --- get the data from the tmCDF file ---
    stl.prgMsg(f'Loading data')
    data_dict_traj = stl.loadDictFromFile(inputFiles)
    data_dict_ECEF_auroral_transform = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\transforms\\{ACESII.fliers[wRocket-4]}\*ECEF_to_auroral*')[0])
    data_dict_LShell = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\Lshell\\{ACESII.fliers[wRocket-4]}\*Lshell.cdf*')[0])
    stl.Done(start_time)

    # --- prepare output ---
    data_dict_output = {
        'N_POS': [np.zeros(shape=np.shape(data_dict_traj['ECEFXPOS'][0])), deepcopy(data_dict_traj['ECEFXPOS'][1])],
        'P_POS': [np.zeros(shape=np.shape(data_dict_traj['ECEFYPOS'][0])), deepcopy(data_dict_traj['ECEFYPOS'][1])],
        'T_POS': [np.zeros(shape=np.shape(data_dict_traj['ECEFZPOS'][0])), deepcopy(data_dict_traj['ECEFZPOS'][1])],
        'T_POS_GRAD': [np.zeros(shape=np.shape(data_dict_traj['ECEFZPOS'][0])), deepcopy(data_dict_traj['ECEFZPOS'][1])],
        'N_POS_GRAD': [np.zeros(shape=np.shape(data_dict_traj['ECEFZPOS'][0])), deepcopy(data_dict_traj['ECEFZPOS'][1])],
        'P_POS_GRAD': [np.zeros(shape=np.shape(data_dict_traj['ECEFZPOS'][0])), deepcopy(data_dict_traj['ECEFZPOS'][1])],
        'N_VEL': [np.zeros(shape=np.shape(data_dict_traj['ECEFXVEL'][0])), deepcopy(data_dict_traj['ECEFXVEL'][1])],
        'P_VEL': [np.zeros(shape=np.shape(data_dict_traj['ECEFYVEL'][0])), deepcopy(data_dict_traj['ECEFYVEL'][1])],
        'T_VEL': [np.zeros(shape=np.shape(data_dict_traj['ECEFZVEL'][0])), deepcopy(data_dict_traj['ECEFZVEL'][1])],
        'Epoch': deepcopy(data_dict_traj['Epoch'])
    }

    #############################################
    # --- CONVERT DATA TO AURORAL COORDINATES ---
    #############################################
    stl.prgMsg('Converting to Auroral Coordinates')

    # form the rocket ECEF position vector
    rkt_pos_ECEF = np.array([data_dict_traj['ECEFXPOS'][0],data_dict_traj['ECEFYPOS'][0],data_dict_traj['ECEFZPOS'][0]]).T

    # for the position vector for the Launch site
    def gps_to_ecef_pyproj(lat, lon, alt):
        ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
        lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
        x, y, z = pyproj.transform(lla, ecef, lon, lat, alt, radians=False)

        return x, y, z


    # DEFINE THE X0 POINT for the position vector
    # Description: Use the L-Shell value right before the High Flyer reached the
    # Aurora as the initial "position" for the auroral coordinates on both flyers.
    L_shell_initial = 8.4
    L_shell_idx = np.abs(data_dict_LShell['L-Shell'][0] - L_shell_initial).argmin()
    Lat_initial = data_dict_LShell['Lat'][0][L_shell_idx]
    Long_initial = data_dict_LShell['Long'][0][L_shell_idx]
    Alt_initial = data_dict_LShell['Alt'][0][L_shell_idx]
    print(Alt_initial)

    launch_ECEF = np.array(gps_to_ecef_pyproj(lat=Lat_initial,lon=Long_initial, alt=Alt_initial))
    position_vector_ECEF = (rkt_pos_ECEF- launch_ECEF)

    # determine the in situ position vector and rotate it into auroral coordinates
    rkt_pos_auroral = np.array([np.matmul(data_dict_ECEF_auroral_transform['ECEF_to_auroral'][0][i], vec) for i,vec in enumerate(position_vector_ECEF)])
    data_dict_output['N_POS'][0] = rkt_pos_auroral[:, 0]
    data_dict_output['T_POS'][0] = rkt_pos_auroral[:, 1]
    data_dict_output['P_POS'][0] = rkt_pos_auroral[:, 2]

    # update the attributes
    data_dict_output['N_POS'][1]['LABLAXIS'] = 'Normal Position from Launch'
    data_dict_output['N_POS'][1]['UNITS'] = 'm'
    data_dict_output['T_POS'][1]['LABLAXIS'] = 'Tangent Position from Launch'
    data_dict_output['T_POS'][1]['UNITS'] = 'm'
    data_dict_output['P_POS'][1]['LABLAXIS'] = 'Field Aligned Position from Launch'
    data_dict_output['P_POS'][1]['UNITS'] = 'm'


    # form the ECEF velocity vector
    rkt_vel_ECEF = np.array([data_dict_traj['ECEFXVEL'][0], data_dict_traj['ECEFYVEL'][0], data_dict_traj['ECEFZVEL'][0]]).T
    rkt_vel_auroral = np.array([np.matmul(data_dict_ECEF_auroral_transform['ECEF_to_auroral'][0][i], vec) for i,vec in enumerate(rkt_vel_ECEF)])
    data_dict_output['N_VEL'][0] = rkt_vel_auroral[:, 0]
    data_dict_output['T_VEL'][0] = rkt_vel_auroral[:, 1]
    data_dict_output['P_VEL'][0] = rkt_vel_auroral[:, 2]

    # update the attributes
    data_dict_output['N_VEL'][1]['LABLAXIS'] = 'Normal Velocity'
    data_dict_output['P_VEL'][1]['LABLAXIS'] = 'Field Aligned Velocity'
    data_dict_output['T_VEL'][1]['LABLAXIS'] = 'Tangent Velocity'
    stl.Done(start_time)

    #######################################################
    # --- Calculate the Position Gradient from Velocity ---
    #######################################################

    # --- get the Epoch in seconds from Launch ---
    Launch_time = pycdf.lib.tt2000_to_datetime(ACESII.launch_T0_TT2000[wRocket-4])
    seconds_from_T0 = stl.EpochTo_T0_Rocket(InputEpoch=data_dict_output['Epoch'][0], T0=Launch_time)
    deltaT = np.array([0] + [seconds_from_T0[i+1] - seconds_from_T0[i] for i in range(len(seconds_from_T0)-1)])
    data_dict_output['T_POS_GRAD'][0] = deepcopy(deltaT)*deepcopy(data_dict_output['T_VEL'][0])
    data_dict_output['T_POS_GRAD'][1]['UNITS'] = 'm'
    data_dict_output['T_POS_GRAD'][1]['LABLAXIS'] = 'Tangent Position Gradient'

    data_dict_output['N_POS_GRAD'][0] = deepcopy(deltaT) * deepcopy(data_dict_output['N_VEL'][0])
    data_dict_output['N_POS_GRAD'][1]['UNITS'] = 'm'
    data_dict_output['N_POS_GRAD'][1]['LABLAXIS'] = 'Normal Position Gradient'

    data_dict_output['P_POS_GRAD'][0] = deepcopy(deltaT) * deepcopy(data_dict_output['P_VEL'][0])
    data_dict_output['P_POS_GRAD'][1]['UNITS'] = 'm'
    data_dict_output['P_POS_GRAD'][1]['LABLAXIS'] = 'Field-Aligned Position Gradient'


    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')

        fileoutName = rf'ACESII_{ACESII.payload_IDs[wRocket-4]}_GPS_trajectory_auroral.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\\trajectories\\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)
        print('\n')




# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
for val in wRocket:
    traj_to_auroral_coordinates(val)