# --- L2_EFI_to_integratedPotential.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: use the



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
import numpy as np
from src.my_imports import *
from scipy.interpolate import CubicSpline
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# ---------------
outputData = True
# ---------------

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none


def L2_EFI_to_integratedPotential():

    # --- --- --- --- --- -
    # --- LOAD THE DATA ---
    # --- --- --- --- --- -

    # --- get the data from the file ---
    stl.prgMsg(f'Loading data')
    data_dict_EFI = stl.loadDictFromFile(r'C:\Data\ACESII\L2\low\ACESII_36364_l2_E_Field_auroral_fullCal.cdf')
    data_dict_traj = stl.loadDictFromFile(r'C:\Data\ACESII\trajectories\low\ACESII_36364_GPS_trajectory_auroral.cdf')
    data_dict_LShell_low = stl.loadDictFromFile('C:\Data\ACESII\coordinates\Lshell\low\ACESII_36364_Lshell.cdf')
    data_dict_LShell_high = stl.loadDictFromFile('C:\Data\ACESII\coordinates\Lshell\high\ACESII_36359_Lshell.cdf')
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
        'Potential': [np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))), {}],
        'Epoch': deepcopy(data_dict_EFI['Epoch']),
        'L-Shell':[np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))), deepcopy(data_dict_LShell_low['L-Shell'][1])],
        'Lat':[np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))), deepcopy(data_dict_LShell_low['Lat'][1])],
        'Long': [np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))), deepcopy(data_dict_LShell_low['Long'][1])],
        'Alt': [np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))), deepcopy(data_dict_LShell_low['Alt'][1])],
    }

    # --- --- --- --- --- --- --- --- --- --- --- --- --- --- -
    # ---Interpolate Trajectory information on EFI timebase ---
    # --- --- --- --- --- --- --- --- --- --- --- --- --- --- -
    # interpolate L-Shell onto EFI timebase
    stl.prgMsg('Interpolating Gradient')
    from scipy.interpolate import CubicSpline
    T0 = dt.datetime(2022, 11, 20, 17, 20)
    time_EFI = stl.EpochTo_T0_Rocket(data_dict_EFI['Epoch'][0], T0=T0)
    stl.Ep
    time_traj = stl.EpochTo_T0_Rocket(data_dict_traj['Epoch'][0], T0=T0)

    for key in data_dict_traj.keys():
        if key != 'Epoch':
            cs = CubicSpline(x=time_traj, y=data_dict_traj[key][0])
            data_dict_traj[key][0] = cs(time_EFI)
    data_dict_traj['Epoch'][0] = deepcopy(data_dict_EFI['Epoch'][0])
    stl.Done(start_time)

    # --- --- --- --- --- --- --- --- --- --- --- --- --- --- -
    # ---Interpolate L-Shell information on EFI timebase ---
    # --- --- --- --- --- --- --- --- --- --- --- --- --- --- -
    stl.prgMsg('Interpolating LShell')

    # --- low flyer ---
    time_EFI = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_EFI['Epoch'][0]])
    time_LShell = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_LShell_low['Epoch'][0]])
    cs = CubicSpline(x=time_LShell, y=data_dict_LShell_low['L-Shell'][0])
    Lshell_EFI = cs(time_EFI)
    data_dict_output['L-Shell'][0] = Lshell_EFI
    stl.Done(start_time)

    # --- --- --- --- --- --- --- --- --- --- --- --- -
    # --- Reduce data to Region JUST outside Aurora ---
    # --- --- --- --- --- --- --- --- --- --- --- --- -

    # pick the initial integration point (i.e. Vref = 0) via high flyer altitude, but
    # correlated via L-Shell value
    high_flyer_alt_limit = 350*stl.m_to_km # in km
    limit_indicies = np.where(data_dict_LShell_high['Alt'][0]>=high_flyer_alt_limit)[0]
    L_lower, L_higher = data_dict_LShell_high['L-Shell'][0][limit_indicies][0],data_dict_LShell_high['L-Shell'][0][limit_indicies][-1]

    Vref_low = np.abs(data_dict_output['L-Shell'][0] - L_lower).argmin()
    Vref_high = np.abs(data_dict_output['L-Shell'][0] - L_higher).argmin() +1

    for key in data_dict_EFI.keys():
        data_dict_EFI[key][0] = deepcopy(data_dict_EFI[key][0][Vref_low:Vref_high])

    for key in data_dict_traj.keys():
        data_dict_traj[key][0] = deepcopy(data_dict_traj[key][0][Vref_low:Vref_high])

    for key in data_dict_output.keys():
        data_dict_output[key][0] = deepcopy(data_dict_output[key][0][Vref_low:Vref_high])

    for key in ['Lat', 'Long', 'Alt']:
        data_dict_output[key] = deepcopy(data_dict_traj[key])

    # --- --- --- --- --- --- --- --- --- --- --- --- ---
    # --- Line Integrate the E-Field to get potential ---
    # --- --- --- --- --- --- --- --- --- --- --- --- ---
    stl.prgMsg('Calculating Line Integral')

    # form the E-Field vector - set the z-component to zero since we dont actually know what this was
    E_Field = np.array([data_dict_EFI['E_N'][0],data_dict_EFI['E_T'][0], np.zeros(shape=(len(data_dict_EFI['Epoch'][0])))]).T

    # Form the position vector
    Pos_vec = np.array([data_dict_traj['N_POS'][0], data_dict_traj['T_POS'][0], data_dict_traj['P_POS'][0]]).T

    # line integrate the result
    from scipy.integrate import simpson

    for i in tqdm(range(len(data_dict_EFI['Epoch'][0]))):

        # N-direction
        N_vals = simpson(x=Pos_vec[:, 0][:i+1], y=E_Field[:, 0][:i+1])

        # T-direction
        T_vals = simpson(x=Pos_vec[:, 1][:i + 1], y=E_Field[:, 1][:i + 1])

        # P-direction
        P_vals = simpson(x=Pos_vec[:, 2][:i + 1], y=E_Field[:, 2][:i + 1])

        # total
        data_dict_output['Potential'][0][i] = -1*(N_vals + T_vals + P_vals)

    data_dict_output['Potential'][1]['UNITS'] = 'V'
    data_dict_output['Potential'][1]['LABLAXIS'] = 'Volts'
    data_dict_output['Potential'][1]['VAR_TYPE'] = 'data'
    data_dict_output['Potential'][1]['DEPEND_0'] = 'Epoch'



    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')
        fileoutName = f'ACESII_{ACESII.payload_IDs[1]}_integrated_potential.cdf'
        outputPath = rf'C:\Data\ACESII\science\integrated_potential\low\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)

# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
L2_EFI_to_integratedPotential()
