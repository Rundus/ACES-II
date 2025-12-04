#--- L1_to_L2_langmuir_fixed.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert the engineering Langmuir data to scientifically useful units.
# Also renames "Boom_Monitor_1/2" --> "Fixed_Boom_Monitor_1/2" etc


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import itertools
# --- --- --- --- ---

import time
from ACESII_code.class_var_func import Done, setupPYCDF

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
wRocket = 4

# --- FIXED PROBE TOGGLES ---
fixedProbe = False

# --- SWEPT PROBE TOGGLES ---
sweptProbe = True

# --- DATA OUTPUT ---
outputData = False


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy

from warnings import filterwarnings # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.
filterwarnings("ignore")
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, L2_ACES_Quick,L2_TRICE_Quick, prgMsg
from glob import glob
from os.path import getsize
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)

from scipy import signal



def LangmuirAnalysis(wRocket, rocketFolderPath, justPrintFileNames, wflyer):

    outputFolderPath = rocketFolderPath + r'science\\Langmuir\\'

    # --- Get ACESII rocket Attributes ---
    rocketAttrs,b,c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'Langmuir'
    L2ModelData = L2_ACES_Quick(wflyer)


    # Set the paths for the file names
    inputFiles = glob(f'{rocketFolderPath}L2\{fliers[wflyer]}\*langmuir*')
    outputFiles = glob(f'{outputFolderPath}\{fliers[wflyer]}\*LangmuirSciData*')

    input_names = [ifile.replace(f'{rocketFolderPath}L2\{fliers[wflyer]}\\', '') for ifile in inputFiles]
    output_names = [ofile.replace(f'{outputFolderPath}{fliers[wflyer]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace('_v00', '') for ifile in input_names]
    output_names_searchable = [ofile.replace('ACESII_', '').replace('36359_', '').replace('36364_', '').replace('LangmuirSciData_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = input_names_searchable[0].replace(f'{rocketFolderPath}L2\{fliers[wflyer]}\\', '').replace('langmuir_','')
    fileoutName = rf'ACESII_{rocketAttrs.rocketID[wRocket-4]}_LangmuirSciData_{dataFile_name}'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Creating langmuir science data for {fliers[wRocket-4]} flyer' + color.END)

        # --- get the data from the tmCDF file ---
        prgMsg('Loading data from input Langmuir Files')
        data_dict = {}

        # Collect the LP data except deltaNdivN into a data dict
        with pycdf.CDF(inputFiles[0]) as inputDataFile:
            for key, val in inputDataFile.items():
                if key not in data_dict:
                    data_dict = {**data_dict, **{key : [inputDataFile[key][...], {key:val for key,val in inputDataFile[key].attrs.items()  }  ]  }  }

        Done(start_time)


        #####################
        # --- Fixed Probe ---
        #####################


        #####################
        # --- Swept Probe ---
        #####################
        if sweptProbe:
            # Apply the calibration function curve
            sweptEpoch = data_dict['Epoch_swept_Current'][0]
            sweptVoltage = data_dict['step_Voltage'][0]
            sweptCurrent = data_dict['swept_Current'][0]
            stepEpoch = data_dict['Epoch_step'][0]

            if wRocket == 4:
                start_plot = 300300
                end_plot =   301300
                # start_plot = 300000
                # end_plot = 308000
            elif wRocket == 5:
                # start_plot = 201900
                # end_plot = 203460
                start_plot = 301500
                end_plot = 303200

            # Look at the data in time before plotting
            plotVsEpoch = False
            if plotVsEpoch:
                xData = np.array(sweptEpoch[start_plot:end_plot])
                yData = np.array(sweptCurrent[start_plot:end_plot])
                I_sat = yData.min() - 1
                yData = np.log(yData - I_sat)
                plt.subplot(2,1,1)
                plt.plot(xData,yData)
                plt.subplot(2, 1, 2)
                plt.plot(stepEpoch[start_plot:end_plot],sweptVoltage[start_plot:end_plot])
                plt.show()

            smoothData = True
            if smoothData:
                # Get the probe data for my region and determine I_probe - I_ion_sat
                xData = np.array(sweptVoltage[start_plot:end_plot])
                yData = np.array(sweptCurrent[start_plot:end_plot])
                I_sat = yData.min()
                yData = yData - I_sat

                #Smooth the data a bit
                w_size = 60
                polyOrder = 3
                smoothing = signal.savgol_filter(yData, window_length=w_size, polyorder=polyOrder, mode="nearest")

                yDataLog = np.log(smoothing)
                yData_smooth = signal.savgol_filter(yDataLog, window_length=w_size, polyorder=polyOrder, mode="nearest")

                # Base data plot
                plt.subplot(1,2,1)
                plt.plot(xData, yData, color='black', label='Base Data')
                plt.plot(xData, smoothing, color='blue', label='Smoothed Base Data')
                plt.legend()

                plt.subplot(1,2,2)
                plt.plot(xData,yDataLog,color='black', label='Log(Smoothed Base Data)')
                plt.plot(xData,yData_smooth,color='red', label='Smoothed Log(Smoothed Data)')
                plt.legend()
                plt.show()

            twoPlotFitTemp = False
            if twoPlotFitTemp:
                # Fit the temperature line
                voltageLow = np.abs(xData - 0.1).argmin()
                voltageHigh = np.abs(xData - 0.7).argmin()

                def line(x, a, b):
                    y = a * x + b
                    return y

                parameters, covariance = scipy.optimize.curve_fit(line, xData[voltageLow:voltageHigh], yData[voltageLow:voltageHigh], maxfev=100000)
                temperXData = np.linspace(xData[voltageLow], xData[voltageHigh], 2000)
                temperYData = [parameters[0] * temperXData[i] + parameters[1] for i in range(len(temperXData))]

                # Base data plot
                plt.subplot(1,2,1)
                plt.plot(xData,yData,color='black',label='Probe Current')
                plt.plot(xData, yData_smooth, color='blue', label='smoothed Current')

                # # Sav-Gol filter plot
                plt.subplot(1, 2, 2)
                plt.xlim(-1, max(xData))
                plt.plot(xData,yData_smooth,color='blue',label='smoothed Current')

                # Linear Temperature Plot
                plt.plot(temperXData, temperYData, color='red',label=f'slope={round(parameters[0],3)},Inter={round(parameters[1],3)}')

                # Electron Saturation Plot

                # Plot labels
                plt.xlabel('Applied Voltage [V]')
                plt.ylabel('$Ln(I_{probe} - I_{sat})$ [Digital]')
                plt.legend()
                plt.show()




        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:

            prgMsg('Creating output file')
            outputPath = f'{outputFolderPath}\\' + f'{fliers[wflyer]}\{fileoutName}'

            # --- delete output file if it already exists ---
            if os.path.exists(outputPath):
                os.remove(outputPath)

            # --- open the output file ---
            with pycdf.CDF(outputPath, '') as LangmuirFile:
                LangmuirFile.readonly(False)

                # --- write out global attributes ---
                inputGlobDic = L2ModelData.cdfFile.globalattsget()
                for key, val in inputGlobDic.items():
                    if key in globalAttrsMod:
                        LangmuirFile.attrs[key] = globalAttrsMod[key]
                    else:
                        LangmuirFile.attrs[key] = val

                # --- WRITE OUT DATA ---
                for varKey, varVal in data_dict.items():
                    if 'Epoch' in varKey: # epoch data
                        LangmuirFile.new(varKey, data=varVal[0], type=33)
                    else: # other data
                        LangmuirFile.new(varKey, data=varVal[0],type=pycdf.const.CDF_REAL8)

                    # --- Write out the attributes and variable info ---
                    for attrKey, attrVal in data_dict[varKey][1].items():
                        if attrKey == 'VALIDMIN':
                            LangmuirFile[varKey].attrs[attrKey] = varVal[0].min()
                        elif attrKey == 'VALIDMAX':
                            LangmuirFile[varKey].attrs[attrKey] = varVal[0].max()
                        elif attrVal != None:
                            LangmuirFile[varKey].attrs[attrKey] = attrVal

                Done(start_time)










# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 0:  # ACES II Integration High
    rocketFolderPath = Integration_data_folder
    wflyer = 0
elif wRocket == 1: # ACES II Integration Low
    rocketFolderPath = Integration_data_folder
    wflyer = 1
elif wRocket == 2:  # TRICE II High
    rocketFolderPath = TRICE_data_folder
    wflyer = 0
elif wRocket == 3: # TRICE II Low
    rocketFolderPath = TRICE_data_folder
    wflyer = 1
elif wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}L2\{fliers[wflyer]}\*langmuir*')) == 0:
    print(color.RED + 'There are no Langmuir files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        LangmuirAnalysis(wRocket, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        LangmuirAnalysis(wRocket, rocketFolderPath, justPrintFileNames,wflyer)