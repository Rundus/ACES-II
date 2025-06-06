# --- L1_to_L1_rocket_convection_cal.py ---
# Description: Create an L2 calibrated EFI dataset
# [1] Interpolate the EFI data onto the RingCore timebase. There's a temporal correction to the E-Field data that Roger provided that must
# be added before interpolation
# [2] Use the CHAOS Bgeo model to subtract the convection electric field from the rocket
# [3] Process the data up to L2



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

def L1_to_L2_rocket_convection_cal(wRocket, justPrintFileNames):

    inputFiles_elec = glob(f'{DataPaths.ACES_data_folder}\\L1\\{ACESII.fliers[wRocket-4]}\*E_Field_ENU_no_corrections*')
    inputFiles_mag = glob(f'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket-4]}\\*RingCore_ENU.cdf*')
    input_names = [ifile.replace(f'{DataPaths.ACES_data_folder}\\L1\\{ACESII.fliers[wRocket-4]}\\', '') for ifile in inputFiles_elec]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles_elec):
            print('[{:.0f}] {:80s}'.format(i, input_names[i], round(getsize(file) / (10 ** 6))))
        return

    # --- get the data from the B-Field file ---
    stl.prgMsg(f'Loading data')
    data_dict_mag = stl.loadDictFromFile(inputFiles_mag[0])
    data_dict_output, GlobalAttrs = stl.loadDictFromFile(inputFiles_elec[0], getGlobalAttrs=True)
    data_dict_traj = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\trajectories\\{ACESII.fliers[wRocket-4]}\\*_ECEF.cdf*')[0])
    stl.Done(start_time)

    ########################################
    # --- ADD IN ROGER'S TIME CORRECTION ---
    ########################################
    stl.prgMsg('Adding In Rogers timebase correction')
    timeCorrection = (0.1157 * 1E9) # in ns
    data_dict_output['Epoch'][0] = deepcopy(np.array([pycdf.lib.tt2000_to_datetime(int(pycdf.lib.datetime_to_tt2000(tme) + timeCorrection)) for tme in data_dict_output['Epoch'][0]]))
    stl.Done(start_time)

    #########################################################
    # --- interpolate E-Field data onto RingCore TimeBase ---
    #########################################################
    stl.prgMsg('Interpolating E-Field Data')
    Epoch_EFI_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_output['Epoch'][0]])
    Epoch_mag_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_mag['Epoch'][0]])

    for key in data_dict_output.keys():
        if key != 'Epoch':
            cs = CubicSpline(Epoch_EFI_tt2000,data_dict_output[key][0])
            data_dict_output[key][0] = deepcopy(cs(Epoch_mag_tt2000))

    data_dict_output['Epoch'][0] = deepcopy(data_dict_mag['Epoch'][0])
    stl.Done(start_time)


    ##############################################
    # --- Correct the Rocket Convection Effect ---
    ##############################################
    stl.prgMsg('Correcting Rocket Convection Effect')

    # Get the payload velocity vector
    rkt_VEL_ECEF = np.array([data_dict_traj['ECEFXVEL'][0], data_dict_traj['ECEFYVEL'][0], data_dict_traj['ECEFZVEL'][0]]).T

    # Convert payload velocity to ENU
    data_dict_transform = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\{ACESII.fliers[wRocket-4]}\\*ECEF_to_ENU.cdf*')[0])
    rkt_VEL_ENU = np.array([np.matmul(data_dict_transform['ECEF_to_ENU'][0][i],vec) for i,vec in enumerate(rkt_VEL_ECEF)])

    # Get the CHAOS geomagnetic field
    B_model = 1E-9*stl.CHAOS(lat=data_dict_traj['Lat'][0],
                        long=data_dict_traj['Long'][0],
                        alt=data_dict_traj['Alt'][0],
                        times=data_dict_traj['Epoch'][0])  # CHAOS in ENU coordinates

    # Calculate the -vxB electric field
    vxB_term = np.array([np.cross(vec,B_model[i]) for i, vec in enumerate(rkt_VEL_ENU)])

    # plot the calibration term
    if Plot_correction_term:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3)

        ax[0].plot(data_dict_output['Epoch'][0], data_dict_output['E_East'][0],color='blue')
        ax[0].plot(data_dict_transform['Epoch'][0], vxB_term[:,0], color='red')

        ax[1].plot(data_dict_output['Epoch'][0], data_dict_output['E_North'][0], color='blue')
        ax[1].plot(data_dict_transform['Epoch'][0], vxB_term[:,1], color='red')

        ax[2].plot(data_dict_output['Epoch'][0], data_dict_output['E_Up'][0], color='blue')
        ax[2].plot(data_dict_transform['Epoch'][0], vxB_term[:,2], color='red')

        for i in range(3):
            ax[i].set_ylim(-0.12,0.12)

        plt.show()
    stl.Done(start_time)

    # interpolate the -vxB term onto the magnetometer timebase
    Epoch_traj_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_traj['Epoch'][0]])

    # X
    cs = CubicSpline(Epoch_traj_tt2000,vxB_term[:,0])
    vxB_E = cs(Epoch_mag_tt2000)

    # Y
    cs = CubicSpline(Epoch_traj_tt2000, vxB_term[:, 1])
    vxB_N = cs(Epoch_mag_tt2000)

    # Z
    cs = CubicSpline(Epoch_traj_tt2000, vxB_term[:, 2])
    vxB_U = cs(Epoch_mag_tt2000)

    # apply the -vxB correction and output the data
    data_dict_output['E_East'][0] = deepcopy( data_dict_output['E_East'][0]- vxB_E)
    data_dict_output['E_North'][0] = deepcopy(data_dict_output['E_North'][0] - vxB_N)
    data_dict_output['E_Up'][0] = deepcopy(data_dict_output['E_Up'][0] - vxB_U)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')

        fileoutName = rf'ACESII_{ACESII.payload_IDs[wRocket-4]}_l1_E_Field_ENU_no_rkt_convec.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\\calibration\\EFI_rkt_convection_calibration\\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output, globalAttrsMod=GlobalAttrs, instrNam='EFI')
        stl.Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
L1_to_L2_rocket_convection_cal(wRocket, justPrintFileNames)