# --- L1_to_L2_offset_fit_cal.py ---
# Description: For the auroral coordinate vxB calibration values,
# [1] determine the GAIN value (y=mx) needed to bring the Tangent component of the
# E-Field into alignment with the vxB term.
# [2] THEN apply the gain correction term to ALL components
# [3] THEN subtract the vxB term to get a fully calibrated E-Field
# [4] Output the EFI in ENU coordinates as L2 data



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

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
outputData = False

# --- TOGGLES ---
smooth_fit_data = False
Plot_fit_data = True
Plot_fitted_data = True
Plot_corrected_data = True


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

def L1_to_L2_offset_fit_cal(wRocket, justPrintFileNames):

    # --- Load the Data ---
    stl.prgMsg(f'Loading data')
    data_dict_vxB = stl.loadDictFromFile(glob('C:\Data\ACESII\calibration\EFI_rkt_convection_calibration\low\*.cdf*')[0])
    stl.Done(start_time)

    # --- --- --- --- --- --- ---
    # --- COLLECTION FIT DATA ---
    # --- --- --- --- --- --- ---
    # Description: Collect a subset of data to calibrate the E-Field data to
    target_time_low = np.abs(data_dict_vxB['Epoch'][0] - dt.datetime(2022, 11, 20, 17, 23, 46) ).argmin()
    target_time_high = np.abs(data_dict_vxB['Epoch'][0] - dt.datetime(2022, 11, 20, 17, 29, 00)).argmin()


    # define a working data dictonary
    data_dict_fit_data = {}
    for key, val in data_dict_vxB.items():
        data_dict_fit_data = {**data_dict_fit_data,
                              f'{key}':[deepcopy(val[0][target_time_low:target_time_high]),deepcopy(val[1])]}

    if smooth_fit_data:
        keys = ['N', 'T', 'p']
        for idx, key in enumerate(keys):
            dat = data_dict_fit_data[f'E_{key}_raw'][0]
            data_dict_fit_data[f'E_{key}_raw'][0] = savgol_filter(dat, window_length=500, polyorder=3)


    if Plot_fit_data:

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3)
        keys = ['N', 'T', 'p']

        for idx,key in enumerate(keys):
            if key == 'N':
                ax[idx].plot(data_dict_fit_data['Epoch'][0], data_dict_fit_data[f'E_{key}_raw'][0], label=f'E_{key}')
                ax[idx].plot(data_dict_fit_data['Epoch'][0], data_dict_fit_data[f'vxB_{key}'][0], label=f'vxB_{key}',
                             color='red')
            elif key == 'T':
                ax[idx].plot(data_dict_fit_data['Epoch'][0], data_dict_fit_data[f'E_{key}_raw'][0], label=f'E_{key}')
                ax[idx].plot(data_dict_fit_data['Epoch'][0], data_dict_fit_data[f'vxB_{key}'][0], label=f'vxB_{key}',
                             color='red')
            else:
                ax[idx].plot(data_dict_fit_data['Epoch'][0], data_dict_fit_data[f'E_{key}_raw'][0], label=f'E_{key}')
                ax[idx].plot(data_dict_fit_data['Epoch'][0], data_dict_fit_data[f'vxB_{key}'][0], label=f'vxB_{key}',
                             color='red')

            ax[idx].legend()

        plt.show()



    # METHOD 1: define the gain function and use scipy.curvefit to make the best fit on the Tangent component
    def fitFunc(x, a):
        return a*x

    # xData = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_fit_data['Epoch'][0]])
    xData = np.array(deepcopy(data_dict_fit_data['E_T_raw'][0]))
    yData = np.array(data_dict_fit_data['vxB_T'][0])


    params, cov = curve_fit(f=fitFunc,
                            xdata=xData,
                            ydata=yData
                           )

    # METHOD 2: METHOD OF LEAST SQUARES
    # Description: Loop over many slope values to see which gives the smallest reduced ChiSquare value for error=1 on all values
    N = 1000
    ChiSquares = []
    slopes = np.linspace(0.1,1.1,N)
    nu = 1/(len(xData) - 1)
    for idx, val in enumerate(slopes):
        ChiSquares.append(nu*np.nansum( np.power(yData-fitFunc(xData,val),2)/np.power(np.std(yData),2)))

    ChiSquares = np.array(ChiSquares)

    # find the best chi square value
    best_fit_idx = np.abs(ChiSquares).argmin()
    best_chival_Chi = ChiSquares[best_fit_idx]
    best_slope_Chi = slopes[best_fit_idx]


    if Plot_fitted_data:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3)


        # get the fit data
        yData_fit = fitFunc(xData,*params)
        yData_Chi = fitFunc(xData,best_slope_Chi)

        # plot the fit and data
        ax[0].plot(xData, yData)
        ax[0].plot(xData,yData_fit,color='red', linewidth=2, label=f'Slope (curve_fit Fitted): {round(params[0],3)}')
        # ax.plot(xData, yData_Chi, linewidth=1,color='blue', label=f'Slope (Chi): {round(best_slope_Chi, 3)}, Chi: {round(best_chival_Chi,8)}')
        ax[0].legend()
        ax[0].set_xlabel('E_T_raw [V/m]')
        ax[0].set_ylabel('vxB_T [V/m]')

        # plot the ratio of vxB_T/E_T
        ratio = yData/xData
        ax[1].plot(data_dict_fit_data['Epoch'][0],ratio,label='vxB_T/E_T')
        ax[1].axhline(y=np.average(ratio),label='(vxB_T/E_T)$_{avg} =$' + f' {round(np.average(ratio),3)}',color='red')
        ax[1].axhline(y=1/np.sqrt(2), label=r'$\frac{1}{\sqrt{2}}$', color='tab:orange')
        ax[1].set_xlabel('Epoch')
        ax[1].set_ylim(0.57,0.73)
        ax[1].set_ylabel('vxB_T/E_T')
        ax[1].legend()

        ratio = data_dict_fit_data['vxB_N'][0]/data_dict_fit_data['E_N_raw'][0]
        ax[2].plot(data_dict_fit_data['Epoch'][0], ratio, label='vxB_N/E_N')
        ax[2].axhline(y=np.nanmean(ratio), label='(vxB_N/E_N)$_{avg} =$' + f' {round(np.average(ratio), 3)}', color='red')
        ax[2].axhline(y=1 / np.sqrt(2), label=r'$\frac{1}{\sqrt{2}}$', color='tab:orange')
        ax[2].set_xlabel('Epoch')
        ax[2].set_ylim(0.5, 8)
        ax[2].set_ylabel('vxB_N/E_N')
        ax[2].legend()

        plt.show()

    # E_N_correction = np.nanmean(data_dict_fit_data['vxB_N'][0] / data_dict_fit_data['E_N_raw'][0])

    E_N_correction = 3.31
    E_T_correction = params[0]

    if Plot_corrected_data:

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3)
        keys = ['N', 'T', 'p']

        for idx, key in enumerate(keys):
            if key == 'N':
                ax[idx].plot(data_dict_fit_data['Epoch'][0], (E_N_correction) * data_dict_fit_data[f'E_{key}_raw'][0], label=f'E_{key}')
            else:
                ax[idx].plot(data_dict_fit_data['Epoch'][0], (E_T_correction) * data_dict_fit_data[f'E_{key}_raw'][0], label=f'E_{key}')

            ax[idx].plot(data_dict_fit_data['Epoch'][0], data_dict_fit_data[f'vxB_{key}'][0], label=f'vxB_{key}',color='red')
            ax[idx].legend()

        plt.show()

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')

        # Apply the gain-correction term
        data_dict_output = {}
        gain_corrections = [E_N_correction,E_T_correction,E_T_correction]
        keys = ['N', 'T', 'p']
        for idx, key in enumerate(keys):
            correct_data = deepcopy(data_dict_vxB[f'E_{key}_raw'][0])*gain_corrections[idx] - deepcopy(data_dict_vxB[f'vxB_{key}'][0])

            data_dict_output = {**data_dict_output, **{f'E_{key}':[correct_data,deepcopy(data_dict_vxB[f'E_{key}_raw'][1])]}}

        data_dict_output = {**data_dict_output,
                            **{'Epoch':deepcopy(data_dict_vxB['Epoch'])}}

        data_dict_output['E_N'][1]['LABLAXIS'] = 'E_Normal'
        data_dict_output['E_T'][1]['LABLAXIS'] = 'E_Tangent'
        data_dict_output['E_p'][1]['LABLAXIS'] = 'E_Field_Aligned'

        E_Field = np.array([deepcopy(data_dict_output['E_N'][0]),deepcopy(data_dict_output['E_T'][0]),deepcopy(data_dict_output['E_p'][0])]).T
        Emag = np.array([np.linalg.norm(val) for val in E_Field])
        data_dict_output = {**data_dict_output,
                            **{'Emag': [Emag,deepcopy(data_dict_output['E_N'][1])]}}
        data_dict_output['Emag'][1]['LABLAXIS'] = 'Emag'

        fileoutName = rf'ACESII_{ACESII.payload_IDs[wRocket-4]}_l2_E_Field_auroral_fullCal.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}\\L2\\{ACESII.fliers[wRocket-4]}\\{fileoutName}'
        stl.outputCDFdata(outputPath, data_dict_output, instrNam='EFI')
        stl.Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
L1_to_L2_offset_fit_cal(wRocket, justPrintFileNames)