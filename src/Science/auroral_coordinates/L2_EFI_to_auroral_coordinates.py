# --- L2_EFI_to_auroral_coordinates.py ---
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
outputData = True

# --- Plots ---
plot_detrend = False
plot_best_fit_rotation = True
plot_interactive_slider = False

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
def L2_to_auroral_coordinates(wflyer, wFile, justPrintFileNames):

    # --- FILE I/O ---
    rocket_folder_path = DataPaths.ACES_data_folder
    rocketID = ACESII.payload_IDs[wflyer]
    globalAttrsMod = deepcopy(ACESII.global_attributes[wflyer])
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'LangmuirData'

    # Set the paths for the file names
    data_repository = f'{input_file_path}\\'
    input_files = glob(data_repository+'*ACESII_36364_l2_E_Field_Field_Aligned*')
    input_names = [ifile.replace(data_repository, '') for ifile in input_files]

    if justPrintFileNames:
        for i, file in enumerate(input_files):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names[i], round(os.path.getsize(file) / (10 ** 6), 1)))
        return

    # --- get the data from the L2 file ---
    stl.prgMsg(f'Loading data')
    data_dict = stl.loadDictFromFile(input_files[wFile])
    data_dict_LShell = stl.loadDictFromFile('C:\Data\ACESII\science\L_shell\low\ACESII_36364_Lshell.cdf')
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {}

    ####################################
    # --- [1] Savgol Filter the Data ---
    ####################################
    for key in ['E_e','E_r','E_p']:
        data_dict[key][0] = savgol_filter(x=data_dict[key][0],
                                            window_length=500,
                                            polyorder=3)


    #####################################################
    # --- [1] Detrend the E-East to remove convection ---
    #####################################################

    # Detrend the E-East componet
    low_idx, high_idx = np.abs(data_dict['Epoch'][0] - dt.datetime(2022, 11, 20, 17, 24, 13)).argmin(), np.abs(data_dict['Epoch'][0] - dt.datetime(2022, 11, 20, 17, 27, 0)).argmin()
    xData_fit = np.array([(pycdf.lib.datetime_to_tt2000(val) - pycdf.lib.datetime_to_tt2000(data_dict['Epoch'][0][0])) / 1E9 for val in data_dict['Epoch'][0][low_idx:high_idx]])
    yData_fit = data_dict['E_e'][0][low_idx:high_idx]
    def fitfunc_detrend(x, a, b):
        return a * x + b

    # get the detrended fit
    params_detrend, cov = curve_fit(f=fitfunc_detrend,
                            xdata=xData_fit,
                            ydata=yData_fit
                            )

    if plot_detrend:
        fig, ax = plt.subplots()
        ax.plot(xData_fit,yData_fit, color='tab:blue')
        ax.plot(xData_fit, fitfunc_detrend(xData_fit,*params_detrend), color='tab:red')
        plt.show()

    #######################################
    # --- [2] Create the E-Field Vector ---
    #######################################

    # Reduce data only to the Auroral regions
    low_idx, high_idx = np.abs(data_dict['Epoch'][0] - dt.datetime(2022,11,20,17,24,48)).argmin(),np.abs(data_dict['Epoch'][0] - dt.datetime(2022, 11, 20, 17, 26,25)).argmin()

    for key, val in data_dict.items():
        data_dict[key][0] = data_dict[key][0][low_idx:high_idx]

    # calculate detrend E-East
    xData_fit = np.array([(pycdf.lib.datetime_to_tt2000(val) - pycdf.lib.datetime_to_tt2000(data_dict['Epoch'][0][0])) / 1E9 for val in data_dict['Epoch'][0]])
    E_East_detrend = deepcopy(data_dict['E_e'][0] - fitfunc_detrend(xData_fit,*params_detrend))


    # Form the E-Field vector
    E_Field = np.array([E_East_detrend, data_dict['E_p'][0], data_dict['E_r'][0]]).T

    ################################
    # --- [3] Fit Detrend E-East ---
    ################################

    def fitfunc(x, a, b):
        return a * x + b

    params, cov = curve_fit(f=fitfunc,
                                    xdata=xData_fit,
                                    ydata=E_East_detrend
                                    )
    # yData_fit = fitfunc(xData_fit,*params)
    yData_fit = fitfunc(xData_fit, *[0,0.01])


    #########################################
    # --- [4] Rotate E-East until Minimum ---
    #########################################

    # define set of rotations and calculate deviation from fitted line
    stl.prgMsg('Rotating Fields')
    N = 100
    angles = np.linspace(-15, 3, N)

    rotation_statistics = np.zeros(shape=(N, 2))
    for idx, ang in tqdm(enumerate(angles)):
        rotated = np.array([np.matmul(stl.Ry(ang), vec) for vec in E_Field])
        E_e_rotated = rotated[:,0]
        test_statistic = np.sum(np.power(E_e_rotated- yData_fit,2))
        rotation_statistics[idx][0] = ang
        rotation_statistics[idx][1] = test_statistic
    stl.Done(start_time)

    lowest_angle = rotation_statistics[:, 0][rotation_statistics[:, 1].argmin()]
    yData_rotated_detrend = np.array([np.matmul(stl.Ry(lowest_angle), vec) for vec in E_Field])
    yData_rotated = np.array([np.matmul(stl.Ry(lowest_angle), vec) for vec in np.array([data_dict['E_e'][0], data_dict['E_p'][0], data_dict['E_r'][0]]).T])

    if plot_interactive_slider:
        from matplotlib.widgets import Slider
        fig, ax = plt.subplots(2)
        ax[0].plot(data_dict['Epoch'][0], data_dict['E_e'][0], color='tab:blue', label='E_e')
        ax[0].plot(data_dict['Epoch'][0], E_E_fit, color='tab:red', label='E_e_fit')
        line1, = ax[0].plot(data_dict['Epoch'][0], data_dict['E_e'][0], color='tab:purple', label='E_e_rotated')
        ax[0].legend()

        ax[1].plot(data_dict['Epoch'][0], data_dict['E_r'][0], color='tab:blue', label='Non-rotated E_r')
        line2, = ax[1].plot(data_dict['Epoch'][0], data_dict['E_r'][0], color='tab:purple', label='E_r_rotated')
        ax[1].axhline(y=0, color='black')
        ax[1].legend()

        # Slider
        axrot = fig.add_axes([0, 0.25, 0.0225, 0.63])
        rotation_slider = Slider(
            ax=axrot,
            label='Angle',
            valmin=-25,
            valmax=5,
            orientation='vertical'
        )

        def update(val):
            rotated = np.array([np.matmul(stl.Ry(rotation_slider.val), vec) for vec in E_Field])
            line1.set_ydata(rotated[:, 0])
            line2.set_ydata(rotated[:, 2])
            fig.canvas.draw_idle()

        rotation_slider.on_changed(update)
        plt.show()

    if plot_best_fit_rotation:
        # Plot the rotating test statistics AND the rotated E-Field data at the lowest angle
        fig, ax = plt.subplots(nrows=2,ncols=2)
        fig.set_figwidth(14)
        fig.set_figheight(10)
        ax[1,0].scatter(rotation_statistics[:, 0], rotation_statistics[:, 1])
        ax[1,0].set_ylabel('Test Statistic')
        ax[1,0].set_xlabel('Angle [deg]')

        ax[1, 1].scatter(yData_rotated[:,0],yData_rotated[:,2], color='black')
        ax[1, 1].set_ylabel('E_r_rotated [mV/m]')
        ax[1, 1].set_xlabel('E_e_rotated [mV/m]')
        ax[1, 1].set_ylim(-0.2,0.2)
        ax[1, 1].set_xlim(-0.2, 0.2)
        ax[1, 1].grid()
        ax[1, 1].axvline(0)
        ax[1, 1].axhline(0)

        ax[0, 0].plot(data_dict['Epoch'][0], yData_fit, color='tab:red', label='fit E_e')
        ax[0, 0].plot(data_dict['Epoch'][0], E_East_detrend, color='tab:blue', label='Non-rotated E_e')
        ax[0, 0].plot(data_dict['Epoch'][0], yData_rotated_detrend[:, 0], color='tab:purple', label=f'Rotated E_e ({lowest_angle} deg)')
        ax[0, 0].set_title('Detrend Rotation Minimization')
        ax[0, 0].legend()

        ax[0,1].plot(data_dict['Epoch'][0], data_dict['E_r'][0], color='tab:blue', label='Non-rotated E_r')
        ax[0,1].plot(data_dict['Epoch'][0], yData_rotated[:, 2], color='tab:purple', label=f'Rotated E_r ({lowest_angle} deg)')
        ax[0,1].axhline(y=0, color='black')
        ax[0,1].legend()
        fig.savefig(r'C:\Data\ACESII\science\auroral_coordinates\low\rotation_analysis.png')


    # --- Interpolate L-Shell onto auroral coordinate timebase ---
    stl.prgMsg('Interpolating LShell')
    from scipy.interpolate import CubicSpline
    T0 = dt.datetime(2022, 11, 20, 17, 20)
    time_EFI = stl.EpochTo_T0_Rocket(data_dict['Epoch'][0], T0=T0)
    time_LShell = stl.EpochTo_T0_Rocket(data_dict_LShell['Epoch'][0], T0=T0)
    cs = CubicSpline(x=time_LShell, y=data_dict_LShell['L-Shell'][0])
    Lshell_EFI = cs(time_EFI)
    stl.Done(start_time)

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:

        exampleVar = {'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -9223372036854775808,
                      'FORMAT': None, 'UNITS': None, 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data',
                      'SCALETYP': 'linear', 'LABLAXIS': None}


        data_dict_output = {**data_dict_output,
                            **{
                                'Epoch': deepcopy(data_dict['Epoch']),
                                'E_tangent':[np.array(yData_rotated[:,0]) ,deepcopy(data_dict['E_e'][1])],
                                'E_p':[np.array(yData_rotated[:,1]),deepcopy(data_dict['E_p'][1])],
                                'E_normal':[np.array(yData_rotated[:,2]),deepcopy(data_dict['E_r'][1])],
                                'rotation_Angle':[np.array([lowest_angle]), deepcopy(exampleVar)],
                                'L-Shell' : [np.array(Lshell_EFI),deepcopy(data_dict_LShell['L-Shell'][1])]
                            }
                            }

        data_dict_output['E_tangent'][1]['LABLAXIS'] = 'E_tangent'
        data_dict_output['E_normal'][1]['LABLAXIS'] = 'E_normal'
        data_dict_output['rotation_Angle'][1]['UNITS'] = 'Deg'


        stl.prgMsg('Creating output file')
        fileoutName_fixed = f'ACESII_{rocketID}_E_Field_Auroral_Coordinates.cdf'
        outputPath = fr'C:\Data\ACESII\science\auroral_coordinates\low\{fileoutName_fixed}'
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
        L2_to_auroral_coordinates(wRocket-4, file_idx, justPrintFileNames)

