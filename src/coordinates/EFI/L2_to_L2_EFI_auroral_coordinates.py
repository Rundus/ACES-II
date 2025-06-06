# --- L2_to_L2_EFI_auroral_coordinates.py ---
# Description: convert EFI ENU coordinates to auroral and output as L2 data



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

# --- Select the Rocket ---
# 4 -> ACES-II High Flier
# 5 -> ACES-II Low Flier
wRocket = 5
wFiles = [0]
remove_AC_flucuations = False
order= 4
freq_cutoff = 0.35
Plot_removed_AC = False

detrend_data = True
Plot_detrend = False
outputData = True


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none

def L2_to_L2_EFI_auroral_coordinates(wRocket):

    # --- get the data from the B-Field file ---
    stl.prgMsg(f'Loading data')
    data_dict_EFI = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket - 4]}\\*E_Field_ENU_fullCal.cdf*')[0])
    data_dict_transform_ENU = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\transforms\\{ACESII.fliers[wRocket - 4]}\\*ECEF_to_ENU.cdf*')[0])
    data_dict_transform_auroral = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\transforms\\{ACESII.fliers[wRocket - 4]}\\*ECEF_to_auroral.cdf*')[0])
    data_dict_LShell = stl.loadDictFromFile(glob(f'{DataPaths.ACES_data_folder}\\coordinates\\Lshell\\{ACESII.fliers[wRocket - 4]}\\*Lshell.cdf*')[0])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
        'E_N' : [np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))),data_dict_EFI['E_East'][1]],
        'E_T': [np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))),data_dict_EFI['E_North'][1]],
        'E_p': [np.zeros(shape=(len(data_dict_EFI['Epoch'][0]))),data_dict_EFI['E_Up'][1]],
        'E_mag':deepcopy(data_dict_EFI['E_mag']),
        'Epoch':deepcopy(data_dict_EFI['Epoch'])
    }

    ###########################################
    # --- Interpolate trasnformation matrix ---
    ###########################################

    Epoch_EFI_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_EFI['Epoch'][0]])
    Epoch_ENU_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_transform_ENU['Epoch'][0]])
    Epoch_FAC_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_transform_auroral['Epoch'][0]])

    # interpolate ENU_to_ECEF matrix onto EFI timebase
    interp_keys = ['a11','a12','a13','a21','a22','a23','a31','a32','a33']
    for key in interp_keys:
        cs = CubicSpline(Epoch_ENU_tt2000, data_dict_transform_ENU[key][0])
        data_dict_transform_ENU[key][0] = deepcopy(cs(Epoch_EFI_tt2000))

    ECEF_to_ENU_matrix = np.array([
        [[data_dict_transform_ENU['a11'][0][i], data_dict_transform_ENU['a12'][0][i], data_dict_transform_ENU['a13'][0][i]],
        [data_dict_transform_ENU['a21'][0][i], data_dict_transform_ENU['a22'][0][i], data_dict_transform_ENU['a23'][0][i]],
        [data_dict_transform_ENU['a31'][0][i], data_dict_transform_ENU['a32'][0][i], data_dict_transform_ENU['a33'][0][i]]]
        for i in range(len(data_dict_EFI['Epoch'][0]))
    ])

    # interpolate ENU_to_ECEF matrix onto EFI timebase
    interp_keys = ['a11', 'a12', 'a13', 'a21', 'a22', 'a23', 'a31', 'a32', 'a33']
    for key in interp_keys:
        cs = CubicSpline(Epoch_FAC_tt2000, data_dict_transform_auroral[key][0])
        data_dict_transform_auroral[key][0] = deepcopy(cs(Epoch_EFI_tt2000))

    ECEF_to_auroral_matrix = np.array([
        [[data_dict_transform_auroral['a11'][0][i], data_dict_transform_auroral['a12'][0][i], data_dict_transform_auroral['a13'][0][i]],
         [data_dict_transform_auroral['a21'][0][i], data_dict_transform_auroral['a22'][0][i], data_dict_transform_auroral['a23'][0][i]],
         [data_dict_transform_auroral['a31'][0][i], data_dict_transform_auroral['a32'][0][i], data_dict_transform_auroral['a33'][0][i]]]
        for i in range(len(data_dict_EFI['Epoch'][0]))
    ])

    ###################################
    # --- Transform the Coordinates ---
    ###################################

    # form the EFI ENU vector
    EFI_ENU = np.array([data_dict_EFI['E_East'][0],data_dict_EFI['E_North'][0],data_dict_EFI['E_Up'][0]]).T
    EFI_ECEF = np.array([np.matmul(ECEF_to_ENU_matrix[i].T, vec) for i, vec in enumerate(EFI_ENU)])
    EFI_auroral = np.array([np.matmul(ECEF_to_auroral_matrix[i], vec) for i, vec in enumerate(EFI_ECEF)])

    data_dict_output['E_N'][0] = EFI_auroral[:, 0]
    data_dict_output['E_N'][1]['LABLAXIS'] = 'Normal Component'

    data_dict_output['E_T'][0] = EFI_auroral[:, 1]
    data_dict_output['E_T'][1]['LABLAXIS'] = 'Tangent Component'

    data_dict_output['E_p'][0] = EFI_auroral[:, 2]
    data_dict_output['E_p'][1]['LABLAXIS'] = 'Field-Aligned Component'

    # --- --- --- --- --- --- --- --- --- --- --
    # --- REDUCE DATA TO ONLY RELEVANT PARTS ---
    # --- --- --- --- --- --- --- --- --- --- --

    # # find the nearest epoch to this LShell
    # target_LShell = 8.4
    # target_time_idx = np.abs(data_dict_LShell['L-Shell'][0] - target_LShell).argmin()
    # target_time = data_dict_LShell['Epoch'][0][target_time_idx]
    #
    # idx = np.abs(data_dict_output['Epoch'][0] - target_time).argmin()
    #
    # for key in data_dict_output.keys():
    #     data_dict_output[key][0] = deepcopy(data_dict_output[key][0][idx:])


    # --- --- --- --- --- --- --- --- ---
    # --- OPTIONALLY Remove AC Trends ---
    # --- --- --- --- --- --- --- --- ---
    if remove_AC_flucuations:
        fs = round(
            1/ (1E-9*(pycdf.lib.datetime_to_tt2000(data_dict_output['Epoch'][0][1]) - pycdf.lib.datetime_to_tt2000(data_dict_output['Epoch'][0][0])))
                    )

        data_keys = ['E_T','E_N','E_p']
        non_filt_T = []
        non_filt_N =[]
        for val in data_keys:

            if val == 'E_T':
                non_filt_T = deepcopy(data_dict_output[val][0])
            elif val == 'E_N':
                non_filt_N = deepcopy(data_dict_output[val][0])
            data_dict_output[val][0] = stl.butter_filter(data=deepcopy(data_dict_output[val][0]),
                          lowcutoff=freq_cutoff,
                          highcutoff=freq_cutoff,
                          filtertype='LowPass',
                          order=order,
                          fs=fs)

        if Plot_removed_AC:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(2)
            ax[0].plot(data_dict_output['Epoch'][0], non_filt_T, color='tab:blue', label='Original',linewidth=3)
            ax[0].plot(data_dict_output['Epoch'][0], data_dict_output['E_T'][0], color='red', label='Filtered',linewidth=1)
            ax[0].set_ylabel('E_T [V/m]')
            ax[0].set_xlabel('Epoch')
            ax[0].set_ylim(-0.1, 0.1)
            ax[0].legend()

            ax[1].plot(data_dict_output['Epoch'][0], non_filt_N, color='tab:blue', label='Original', linewidth=3)
            ax[1].plot(data_dict_output['Epoch'][0], data_dict_output['E_N'][0], color='red', label='Filtered', linewidth=1)
            ax[1].set_ylabel('E_N [V/m]')
            ax[1].set_xlabel('Epoch')
            ax[1].set_ylim(-0.1, 0.1)
            ax[1].legend()

            plt.show()

    # --- --- --- --- --- --- --- --- ---
    # --- OPTIONALLY DETREND E_tangent ---
    # --- --- --- --- --- --- --- --- ---

    if detrend_data:
        from scipy.signal import detrend

        non_detrend_T = deepcopy(data_dict_output['E_T'][0])
        data_dict_output['E_T'][0] = detrend(
            data=data_dict_output['E_T'][0],
            type='linear',
        )


        non_detrend_N = deepcopy(data_dict_output['E_N'][0])
        data_dict_output['E_N'][0] = detrend(
            data=data_dict_output['E_N'][0],
            type='constant',
        )

        if Plot_detrend:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(2)
            ax[0].plot(data_dict_output['Epoch'][0], data_dict_output['E_T'][0], color='red', label='Detrended')
            ax[0].plot(data_dict_output['Epoch'][0], non_detrend_T, color='tab:blue', label='Original')
            ax[0].set_ylabel('E_T [V/m]')
            ax[0].set_xlabel('Epoch')
            ax[0].set_ylim(-0.1, 0.1)
            ax[0].legend()

            ax[1].plot(data_dict_output['Epoch'][0], data_dict_output['E_N'][0], color='red', label='Detrended')
            ax[1].plot(data_dict_output['Epoch'][0], non_detrend_N, color='tab:blue', label='Original')
            ax[1].set_ylabel('E_N [V/m]')
            ax[1].set_xlabel('Epoch')
            ax[1].set_ylim(-0.1, 0.1)
            ax[1].legend()

            plt.show()


    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')

        fileoutName = rf'ACESII_{ACESII.payload_IDs[wRocket-4]}_l2_E_Field_auroral_fullCal.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output, instrNam='EFI')
        stl.Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
L2_to_L2_EFI_auroral_coordinates(wRocket)