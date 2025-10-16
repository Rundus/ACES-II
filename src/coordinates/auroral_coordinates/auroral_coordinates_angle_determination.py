# --- auroral_coordinates_angle_determination.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:
# [1] Develop an algorthim to center the E-East around 0 mV/m
# [2] develop a method to rotate E-Field, field-aligned coordinate data in order to minimize
# the perpendicular E-Field. This is the "slab" approximation. Then, plot over All Sky imager to compare

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

from src.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False # Just print the names of files

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [0]

input_file_path = 'C:\Data\ACESII\L2\low'

# --- OutputData ---
outputData = False

# --- Plots ---
plot_interactive_slider = True
plot_best_fit_rotation = True


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from src.my_imports import *
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#######################
# --- MAIN FUNCTION ---
#######################
def auroral_coordinates_angle_determination(wflyer, wFile, justPrintFileNames):

    # --- FILE I/O ---
    rocket_folder_path = DataPaths.ACES_data_folder
    rocketID = ACESII.payload_IDs[wflyer]
    globalAttrsMod = deepcopy(ACESII.global_attributes[wflyer])
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'LangmuirData'

    # Set the paths for the file names
    data_repository = f'{input_file_path}\\'
    input_files = glob(data_repository+'*EFI_FAC_fullCal.cdf*')
    input_names = [ifile.replace(data_repository, '') for ifile in input_files]

    if justPrintFileNames:
        for i, file in enumerate(input_files):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names[i], round(os.path.getsize(file) / (10 ** 6), 1)))
        return

    # --- get the data from the L2 file ---
    stl.prgMsg(f'Loading data')
    data_dict_EFI = stl.loadDictFromFile(input_files[wFile])
    data_dict_LShell = stl.loadDictFromFile('C:\Data\ACESII\coordinates\Lshell\low\ACESII_36364_Lshell.cdf')
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}

    ####################################
    # --- [1] Savgol Filter the Data ---
    ####################################
    for key in ['E_e', 'E_r', 'E_p']:
        data_dict_EFI[key][0] = savgol_filter(x=data_dict_EFI[key][0],
                                            window_length=500,
                                            polyorder=3)


    #######################################
    # --- [2] Create the E-Field Vector ---
    #######################################

    # Define the range of L-Shells to denote the "Auroral Region"
    L_Shell_range = [8.4, 9.75]

    # interpolate L-Shell onto EFI timebase
    stl.prgMsg('Interpolating LShell')
    from scipy.interpolate import CubicSpline
    T0 = dt.datetime(2022, 11, 20, 17, 20)
    time_EFI = stl.EpochTo_T0_Rocket(data_dict_EFI['Epoch'][0], T0=T0)
    time_LShell = stl.EpochTo_T0_Rocket(data_dict_LShell['Epoch'][0], T0=T0)
    cs = CubicSpline(x=time_LShell, y=data_dict_LShell['L-Shell'][0])
    Lshell_EFI = cs(time_EFI)
    low_idx = np.abs(Lshell_EFI - L_Shell_range[0]).argmin()
    high_idx = np.abs(Lshell_EFI - L_Shell_range[1]).argmin()
    stl.Done(start_time)

    # Reduce data only to the Auroral regions
    for key, val in data_dict_EFI.items():
        data_dict_EFI[key][0] = data_dict_EFI[key][0][low_idx:high_idx]

    # Form the E-Field vector
    E_Field = np.array([data_dict_EFI['E_r'][0], data_dict_EFI['E_e'][0], data_dict_EFI['E_p'][0]]).T

    ###########################################################
    # --- [3] Rotate E-East until Maximum in hodogram slope ---
    ###########################################################

    # define set of rotations and calculate deviation from fitted line
    stl.prgMsg('Rotating Fields')
    N = 5
    angles = np.linspace(-20, 20, N)

    # fit a linear line to the E_e compoennt
    def fitFunc(x, a, b):
        return a*x+b

    rotation_statistics = np.zeros(shape=(N, 3))
    for idx, ang in tqdm(enumerate(angles)):

        # rotate the vectors
        rotated = np.array([np.matmul(stl.Rz(ang), vec) for vec in E_Field])

        # sort the hodogram data based on E_e
        E_e_rotated, E_r_rotated = zip(*sorted(zip(rotated[:, 1], rotated[:, 0])))

        params, cov = curve_fit(fitFunc, xdata=E_e_rotated, ydata=E_r_rotated)

        rotation_statistics[idx][0] = ang
        rotation_statistics[idx][1] = params[0]
        rotation_statistics[idx][2] = params[1]
    stl.Done(start_time)

    best_fit_idx = np.abs(rotation_statistics[:, 1]).argmax()
    best_angle = rotation_statistics[:, 0][best_fit_idx]
    best_slope = rotation_statistics[:, 1][best_fit_idx]
    best_intercept = rotation_statistics[:, 2][best_fit_idx]

    choice_idx = np.abs(rotation_statistics[:, 0]-(-9)).argmin()
    print(rotation_statistics[:, 1][choice_idx],rotation_statistics[:, 2][choice_idx])
    print(best_slope,best_intercept)
    yData_rotated = np.array([np.matmul(stl.Rz(best_angle), vec) for vec in np.array([data_dict_EFI['E_r'][0], data_dict_EFI['E_e'][0], data_dict_EFI['E_p'][0]]).T])

    if plot_interactive_slider:
        from matplotlib.widgets import Slider
        fig, ax = plt.subplots(2)
        ax[0].plot(data_dict_EFI['Epoch'][0], data_dict_EFI['E_r'][0], color='tab:blue', label='E_r')
        # ax[0].plot(data_dict_EFI['Epoch'][0], yData_fit, color='tab:red', label='E_e_fit')
        line1_r, = ax[0].plot(data_dict_EFI['Epoch'][0], data_dict_EFI['E_r'][0], color='tab:purple', label='E_r_rotated')
        ax[0].legend()
        ax[0].set_ylim(-0.15, 0.15)
        ax[0].set_ylabel('E_r')

        ax[1].plot(data_dict_EFI['Epoch'][0], data_dict_EFI['E_e'][0], color='tab:blue', label='Non-rotated E_e')
        line2_e, = ax[1].plot(data_dict_EFI['Epoch'][0], data_dict_EFI['E_e'][0], color='tab:purple', label='E_e_rotated')
        ax[1].axhline(y=0, color='black')
        ax[1].legend()
        ax[1].set_ylim(-0.15, 0.15)
        ax[1].set_ylabel('E_e')

        # Slider
        axrot = fig.add_axes([0, 0.25, 0.0225, 0.63])
        rotation_slider = Slider(
            ax=axrot,
            label='Angle',
            valmin=-180,
            valmax=180,
            orientation='vertical'
        )

        def update(val):
            rotated = np.array([np.matmul(stl.Rz(rotation_slider.val), vec) for vec in E_Field])
            line1_r.set_ydata(rotated[:, 0])
            line2_e.set_ydata(rotated[:, 1])
            fig.canvas.draw_idle()

        rotation_slider.on_changed(update)
        plt.show()

    if plot_best_fit_rotation:
        # Plot the rotating test statistics AND the rotated E-Field data at the lowest angle
        fig, ax = plt.subplots(nrows=2, ncols=2)
        fig.set_figwidth(14)
        fig.set_figheight(10)

        # ax[0, 0].plot(data_dict_EFI['Epoch'][0], yData_fit, color='tab:red', label='fit E_e')
        ax[0, 0].plot(data_dict_EFI['Epoch'][0], data_dict_EFI['E_e'][0], color='tab:blue', label='Non-rotated E_e')
        ax[0, 0].plot(data_dict_EFI['Epoch'][0], yData_rotated[:, 1], color='tab:purple', label='Non-rotated E_e')
        ax[0, 0].set_title('Detrend Rotation Minimization')
        ax[0, 0].legend()

        ax[0, 1].plot(data_dict_EFI['Epoch'][0], data_dict_EFI['E_r'][0], color='tab:blue', label='Non-rotated E_r')
        ax[0, 1].plot(data_dict_EFI['Epoch'][0], yData_rotated[:, 0], color='tab:purple', label=f'Rotated E_r ({best_angle} deg)')
        ax[0, 1].axhline(y=0, color='black')
        ax[0, 1].legend()

        ax[1,0].scatter(rotation_statistics[:, 0], np.abs(rotation_statistics[:, 1]))
        ax[1,0].set_ylabel('Test Statistic')
        ax[1,0].set_xlabel('Angle [deg]')

        ax[1, 1].scatter(yData_rotated[:,1], yData_rotated[:,0], color='blue', label='rotated')
        ax[1, 1].scatter(data_dict_EFI['E_e'][0], data_dict_EFI['E_r'][0], color='black',alpha=0.5, label='non-rotated')
        ax[1, 1].axline((0,best_intercept), slope=best_slope,color='red', label='Fit', linestyle='--')
        ax[1, 1].set_ylabel('E_r_rotated [V/m]')
        ax[1, 1].set_xlabel('E_e_rotated [V/m]')
        ax[1, 1].set_ylim(-0.05,0.15)
        ax[1, 1].set_xlim(-0.1, 0.1)
        ax[1, 1].grid()
        ax[1, 1].axvline(0)
        ax[1, 1].axhline(0)
        ax[1,1].legend()

        fig.savefig(r'C:\Data\ACESII\coordinates\auroral_coordinates\low\rotation_analysis.png')

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        exampleVar = {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -9223372036854775808,
                      'FORMAT': None, 'UNITS': None, 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data',
                      'SCALETYP': 'linear', 'LABLAXIS': None}


        data_dict_output = {**data_dict_output,
                            **{
                                'Epoch': deepcopy(data_dict_EFI['Epoch']),
                                'E_tangent':[np.array(yData_rotated[:,0]) ,deepcopy(data_dict_EFI['E_e'][1])],
                                'E_p':[np.array(yData_rotated[:,1]),deepcopy(data_dict_EFI['E_p'][1])],
                                'E_normal':[np.array(yData_rotated[:,2]),deepcopy(data_dict_EFI['E_r'][1])],
                                'rotation_Angle':[np.array([best_angle]), deepcopy(exampleVar)],
                                'L-Shell' : [np.array(Lshell_EFI),deepcopy(data_dict_LShell['L-Shell'][1])]
                            }
                            }

        data_dict_output['E_tangent'][1]['LABLAXIS'] = 'E_tangent'
        data_dict_output['E_normal'][1]['LABLAXIS'] = 'E_normal'
        data_dict_output['rotation_Angle'][1]['UNITS'] = 'Deg'


        stl.prgMsg('Creating output file')
        fileoutName_fixed = f'ACESII_{rocketID}_auroral_coordinates_angle.cdf'
        outputPath = fr'C:\Data\ACESII\coordinates\auroral_coordinates\low\{fileoutName_fixed}'
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
target_files = f'C:\Data\ACESII\L2\low\*.cdf'

if len(glob(target_files)) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    for file_idx in wFiles:
        auroral_coordinates_angle_determination(wRocket-4, file_idx, justPrintFileNames)

