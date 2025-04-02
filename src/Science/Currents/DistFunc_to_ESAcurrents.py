# --- L0_to_L1.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Convert electrostatic analyzer data from diffNFlux to Distribution Function



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"


import itertools
# --- --- --- --- ---
import time
from ACESII_code.myImports import *

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
wRocket = 4

# select which files to convert
# [] --> all files
# [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = []


# Number of points to interpolate between coordinate values in the velocity space coordinate change
N = 500
fillvalue = 0 # value to place
EngyRange = [12, 40] # nominally [12,40]. Reduces the energy range to interpolate over since we don't really reach the higher energies ever

pitchBound = [-10,10] # boundaries of pitch angle for acceptence to calculate J_parallel
energyBound = 7 # same as above but energy. In eV. Nominally 7eV

useSpecificLocations = False
locations = [i for i in range(0,9000,100)]


reduceData=True
targetTimes= [dt.datetime(2022,11,20,17,22,30,00),dt.datetime(2022,11,20,17,27,30,00)]

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from warnings import filterwarnings # USED TO IGNORE WARNING ABOUT "UserWarning: Invalid dataL1 type for dataL1.... Skip warnings.warn('Invalid dataL1 type for dataL1.... Skip')" on Epoch High dataL1.
filterwarnings("ignore")
from tqdm import tqdm
from scipy.interpolate import LinearNDInterpolator
from ACESII_code.missionAttributes import ACES_mission_dicts, TRICE_mission_dicts
from ACESII_code.data_paths import Integration_data_folder, ACES_data_folder, TRICE_data_folder, fliers
from ACESII_code.class_var_func import color, L2_ACES_Quick,L2_TRICE_Quick, prgMsg, cm_to_m,q0,IonMasses,m_e, outputCDFdata
from glob import glob
from os.path import getsize
setupPYCDF()
from spacepy import pycdf
pycdf.lib.set_backward(False)


def DistFunc_to_ESAcurrents(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):
    rocketFolderPath = rocketFolderPath + r'science\\'
    outputFolderPath = rocketFolderPath + r'ESA_currents\\'

    # --- ACES II Flight/Integration Data ---
    rocketAttrs,b,c = ACES_mission_dicts()
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'ESA_currents'
    L2ModelData = L2_ACES_Quick(wflyer)

    # Set the paths for the file names
    inputFiles = glob(f'{rocketFolderPath}DistFunc\{fliers[wflyer]}\*.cdf')
    outputFiles = glob(f'{outputFolderPath}ESA_Currents\{fliers[wflyer]}\*.cdf')

    input_names = [file.replace(f'{rocketFolderPath}DistFunc\{fliers[wflyer]}\\', '') for file in inputFiles]
    output_names = [file.replace(f'{outputFolderPath}ESA_Currents\{fliers[wflyer]}\\', '') for file in outputFiles]

    input_names_searchable = [file.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace('distFunc_', '').replace('_v00', '').replace('__', '_') for file in input_names]
    output_names_searchable = [file.replace('ACES_', '').replace('36359_', '').replace('36364_', '').replace('esaCurrents_', '').replace('_v00', '') for file in output_names]

    dataFile_name = input_names[wFile].replace(f'{rocketFolderPath}\{fliers[wflyer]}\\', '')
    fileoutName = dataFile_name.replace('distFunc', 'jParallel')

    if justPrintFileNames:
            for i, file in enumerate(inputFiles):
                anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
                print('[{:.0f}] {:70s}{:5.1f} MB   Made ESACurrents: {:3s} '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1), anws[0]))
    else:
        print('\n')
        print(color.UNDERLINE + f'Calculating J_parallel for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the tmCDF file ---
        prgMsg('Loading data from distribution function Files')

        data_dict = loadDictFromFile(inputFiles[wFile],{},reduceData=reduceData,targetTimes=targetTimes,wKeys=['Distribution_Function','Epoch'])


        # --- --- --- --- --- --- --- --
        # --- CALCULATE ESA CURRENTS ---
        # --- --- --- --- --- --- --- --

        Done(start_time) ; prgMsg('Performing Distribution Function Coordinate Change')

        # The distribution function data must be converted into a new coordinate space, v_perp & v_para
        # Since we ultimately wish to integrate the data, we should interpolate between the non-uniform spacing of v_perp & v_para to get a nice grid
        # We can then integrate the interpolated distribution function

        # --- organize/convert the data ---
        distFunc = data_dict['Distribution_Function'][0]
        DistFuncs= []
        Pitch = data_dict['Pitch_Angle'][0]
        Energy = data_dict['Energy'][0]

        # The vperp and vparallel coordinates never change, so pre-calculate them here for use later
        Vpars = [np.cos(np.radians(Pitch[ptch])) * np.sqrt(2 * q0 * Energy[engy] / m_e) for engy in range(EngyRange[0], EngyRange[1]) for ptch in range(len(distFunc[0]))]
        Vperps = [np.sin(np.radians(Pitch[ptch])) * np.sqrt(2 * q0 * Energy[engy] / m_e) for engy in range(EngyRange[0], EngyRange[1]) for ptch in range(len(distFunc[0]))]

        for tme in range(len(distFunc)):
            tempDists = [distFunc[tme][ptch][engy] for engy in range(EngyRange[0],EngyRange[1]) for ptch in range(len(distFunc[0]))]
            DistFuncs.append(tempDists)

        Vpars = np.array(Vpars) ; Vperps = np.array(Vperps) ; DistFuncs = np.array(DistFuncs)

        Done(start_time) ; prgMsg('Interpolating Dist Func')

        # interpolate Data
        X = np.linspace(Vpars.min(), Vpars.max(), N)
        Y = np.linspace(Vperps.min(), Vperps.max(), N)
        X, Y = np.meshgrid(X, Y, indexing='ij')
        Z = []

        # interpolate all the data

        if useSpecificLocations:
            iterateThis = locations
        else:
            iterateThis = range(len(distFunc))

        for i in tqdm(iterateThis):
            interp = LinearNDInterpolator(list(zip(Vpars, Vperps)), DistFuncs[i], fill_value=fillvalue)
            Z.append(interp(X, Y))

        # --- --- --- --- --- --- --- ---
        # --- Trapezoidal Integration ---
        # --- --- --- --- --- --- --- ---
        Done(start_time) ; prgMsg('Calculating J_parallel \n')

        j_para = []
        for tme in tqdm(iterateThis):
            zData = Z[tme]

            ##########################################################################
            # --- Filter out z-data that is within the pitch and energy boundaries ---
            ##########################################################################
            par_Strips = []
            perp_Strips = []
            z_Strips = []

            for parIndex in range(N): # loop over x axis
                par_temp = []
                perp_temp = []
                z_temp = []

                for perpIndex in range(N): # loop over y axis:

                    particlePitch = np.degrees(np.arctan(Y[0][perpIndex] /X[parIndex][0] )) # particle's pitch angle
                    particleEnergy =  (0.5*m_e*((Y[0][perpIndex]**2) + (X[parIndex][0]**2)))/q0 # particle's energy in eV

                    if Y[0][perpIndex] < 0: # if you're below Vperp = 0
                        # if you're pitch is outside the -10 to 190deg range
                        if (particlePitch <= pitchBound[0]) or (particlePitch >= pitchBound[1]):
                            par_temp.append(X[parIndex][0])
                            perp_temp.append(Y[0][perpIndex])
                            z_temp.append(0)
                        else: # if you're within the correct, pitch
                            if particleEnergy >= energyBound:  # place the data if it's good
                                par_temp.append(X[parIndex][0])
                                perp_temp.append(Y[0][perpIndex])
                                z_temp.append(zData[parIndex][perpIndex])
                            else:  # if its too low energy
                                par_temp.append(X[parIndex][0])
                                perp_temp.append(Y[0][perpIndex])
                                z_temp.append(0)
                    else:
                        if particleEnergy >= energyBound: # place the data if it's good
                            par_temp.append(X[parIndex][0])
                            perp_temp.append(Y[0][perpIndex])
                            z_temp.append(zData[parIndex][perpIndex])
                        else: # if it's bad do a fillvalue
                            par_temp.append(X[parIndex][0])
                            perp_temp.append(Y[0][perpIndex])
                            z_temp.append(0)

                par_Strips.append(par_temp)
                perp_Strips.append(perp_temp)
                z_Strips.append(z_temp)

            ###############################
            # --- Integrate over Strips ---
            ###############################

            # Method: Trapezoid rule: the more efficent form:G_n =  deltaX /2 (f(x0) + 2f(x1) + ... + 2f(X_N-1) + f(x_N))

            G = []

            deltaVperp = (perp_Strips[0][1] - perp_Strips[0][0])/2
            deltaVpar = (par_Strips[1][0] - par_Strips[0][0])/2

            for n in range(len(par_Strips)): # for each par strip

                firstVal = perp_Strips[n][0]*par_Strips[n][0]*z_Strips[n][-1]
                lastVal = perp_Strips[n][-1]*par_Strips[n][-1]*z_Strips[n][-1]

                # Adjust the first val if it's within the 190deg or -10deg region
                if perp_Strips[n][0] < 0:
                    firstVal = -1 * firstVal

                # loop over the rest of the values, changing the multiplication if I'm in the V_perp <0 region
                middleSum = sum([perp_Strips[n][k]*par_Strips[n][k]*z_Strips[n][k] if perp_Strips[n][k] >= 0 else -1*perp_Strips[n][k]*par_Strips[n][k]*z_Strips[n][k]  for k in range(1,len(par_Strips[n])-1)])

                G.append(deltaVperp * (firstVal + middleSum + lastVal))

            # calculate j_parallel by one more trapizoidal integration_tad_files
            j_para_temp = -2*np.pi*q0*sum([deltaVpar*(G[i+1] + G[i]) for i in range(len(G) - 1)])
            if np.abs(j_para_temp) >= 10**(-6):
                j_para_temp = rocketAttrs.epoch_fillVal

            j_para.append(j_para_temp)

        j_para = np.array(j_para)

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        Done(start_time) ; prgMsg('Creating output file')

        outputPath = f'{outputFolderPath}{fliers[wflyer]}\\{fileoutName}'

        data_dict = {**data_dict, **{'J_parallel':
                                         [j_para, {'LABLAXIS': 'J_parallel',
                                                   'DEPEND_0': 'Epoch_esa',
                                                   'DEPEND_1': None,
                                                   'DEPEND_2': None,
                                                   'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                   'UNITS': '!N A!N m!U-2!N',
                                                   'VALIDMIN': j_para.min(), 'VALIDMAX': j_para.max(),
                                                   'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}

        del data_dict['Distribution_Function']

        if useSpecificLocations:
            newEpoch = [data_dict['Epoch_esa'][0][i] for i in iterateThis]
            data_dict['Epoch_esa'][0] = np.array(newEpoch)



        # ADD IN THE traj DATA
        prgMsg('Interpolating trajectory data into output file')
        inputFiles_traj = glob(f'{rocketFolderPath}trajectories\{fliers[wflyer]}\*.cdf')
        data_dict_traj = loadDictFromFile(inputFiles_traj[0],{},reduceData=False,targetTimes=[data_dict['Epoch'][0][0],data_dict['Epoch'][0][-1]])
        data_dict_trajInterp = InterpolateDataDict(InputDataDict=data_dict_traj,
                                                       InputEpochArray=data_dict_traj['Epoch'][0],
                                                       wKeys=['geomagLat','geomagLong','geomagAlt'],
                                                       targetEpochArray=data_dict['Epoch'][0])

        data_dict = {**data_dict, **{'geomagLat': data_dict_trajInterp['geomagLat']}}
        data_dict = {**data_dict, **{'geomagLong': data_dict_trajInterp['geomagLong']}}
        data_dict = {**data_dict, **{'geomagAlt': data_dict_trajInterp['geomagAlt']}}
        outputCDFdata(outputPath, data_dict, L2ModelData, globalAttrsMod)


        Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
if wRocket == 4:  # ACES II High
    wflyer = 0
elif wRocket == 5: # ACES II Low
    wflyer = 1

if len(glob(f'{rocketFolderPath}science\DistFunc\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        DistFunc_to_ESAcurrents(wRocket, 0, rocketFolderPath, justPrintFileNames, wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}science\DistFunc\{fliers[wflyer]}\*.cdf')))):
            DistFunc_to_ESAcurrents(wRocket, fileNo, rocketFolderPath, justPrintFileNames, wflyer)
    else:
        for filesNo in wFiles:
            DistFunc_to_ESAcurrents(wRocket, filesNo, rocketFolderPath, justPrintFileNames, wflyer)