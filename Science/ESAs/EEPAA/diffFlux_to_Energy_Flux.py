# --- diffFlux_to_Energy_Flux.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: using the specs of the ACESII ESAs, convert from differential Energy Flux
# to just energy flux as described in EEPAA_Flux_Conversion.pdf document in Overleaf



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

from myImports import *
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False
wRocket = 4
wFiles = [0]
modifier = ''
inputPath_modifier = 'l2' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier = 'l3\Energy_Flux' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

# ---------------------------
outputData = True
# ---------------------------
downwardPitchRange = [1,2+1] # what pitch indicies to consider when calculating parallel (downward) # 0-90deg
upwardPitchRange = [18,19+1] # 90-180deg

# solid angle contributions of each pitch bin (starting at 0deg to -190 deg)
solidAngleHfactor =[
    204.4, # -10deg
    25.714, # 0 deg
    204.4, # 10 deg
    402.9, # 20 deg
    589, # 30
    757.1, # 40
    902.3, # 50
    1020.1, # 60
    1106.9, # 70
    1160, # 80
    1177.9, # 90
    1160, # 100
    1106.9, # 110
    1020.1, # 120
    902.3, # 130
    757.1, # 140
    589, # 150 deg
    402.9, # 160 deg
    204.4, # 170deg
    25.714, # 180deg
    204.4 # 190 deg
]



# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from numpy import trapz


def diffFlux_to_Energy_Flux(wRocket, rocketFolderPath, justPrintFileNames, wflyer, wfile):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    ModelData = L2_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf')
    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    input_names_searchable = [ifile.replace(inputPath_modifier.lower() +'_', '') for ifile in input_names]

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in inputFiles[wfile]:
            descriptiorNam = ['EEPAA', 'LEESA', 'IEPAA', 'Langmuir']
            wInstr = [index, instr, descriptiorNam[index]]
            break
        else:
            wInstr = [0, 'attitude', '']

    fileoutName = f'ACESII_{rocketID}_{wInstr[1].lower()}_Energy_Flux.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')

        # --- --- --- --- --- -
        # --- LOAD THE DATA ---
        # --- --- --- --- --- -

        # --- get the data from the file ---
        prgMsg(f'Loading data from ESA Files')
        # load data here
        data_dict, globalAttrs = loadDictFromFile(inputFilePath=inputFiles[wfile], getGlobalAttrs=True)
        Done(start_time)

        # --- --- --- --- --- --- --- -
        # --- INTEGRATE OVER ENERGY ---
        # --- --- --- --- --- --- --- -

        prgMsg('Calculating Energy Flux')

        Epoch = data_dict['Epoch'][0]
        Pitch = data_dict['Pitch_Angle'][0]
        Energy = data_dict['Energy'][0][::-1]
        diffEflux = data_dict['Differential_Energy_Flux'][0]

        diffEFlux_avg = [ [[] for ptch in Pitch] for tme in Epoch]
        Eflux = [ [[] for ptch in Pitch] for tme in Epoch]

        # perform the integration
        for tme in tqdm(range(len(Epoch))):
            for ptch in range(len(Pitch)):
                diffFLux_data = diffEflux[tme][ptch][::-1]

                if all(np.array(diffFLux_data) >= 0): # see if there's any fillvals or negative values in array. If not, then you can integrate
                    integratedValue = trapz(y=diffFLux_data, x=Energy)
                    EFluxVal = (1/erg_to_eV) * integratedValue * solidAngleHfactor[ptch]
                else:
                    EFluxVal = rocketAttrs.epoch_fillVal
                    integratedValue = rocketAttrs.epoch_fillVal

                diffEFlux_avg[tme][ptch].append(integratedValue)
                diffEFlux_avg[tme][ptch] = diffEFlux_avg[tme][ptch][0] # reduce each element from a list to a single value
                Eflux[tme][ptch].append(EFluxVal)
                Eflux[tme][ptch] = Eflux[tme][ptch][0] # reduce each element from a list to a single value

        diffEFlux_avg = np.array(diffEFlux_avg)
        Eflux = np.array(Eflux)
        Done(start_time)


        # Calculate Upward and Downward Particle Energy Flux
        EFlux_downward = np.zeros(shape=(len(Eflux)))
        EFlux_Upward = np.zeros(shape=(len(Eflux)))
        cosineContribution = np.array([np.cos(np.radians(10*(i-1))) for i in range(21)])

        for tme in range(len(Epoch)):

            # remove fillval contributions
            EFluxArray = Eflux[tme]

            EFluxArray[EFluxArray ==rocketAttrs.epoch_fillVal ] = 0
            EFlux_contribution = np.array(Eflux[tme]) * cosineContribution
            EFlux_Upward[tme] = -1*sum(EFlux_contribution[upwardPitchRange[0]:upwardPitchRange[1]])
            EFlux_downward[tme] = sum(EFlux_contribution[downwardPitchRange[0]:downwardPitchRange[1]])



        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:
            prgMsg('Creating output file')

            data_dict_output = {'Energy_Flux': [Eflux, deepcopy(data_dict['Differential_Energy_Flux'][1])],
                                'Differential_Energy_Flux_avg': [diffEFlux_avg, deepcopy(data_dict['Differential_Energy_Flux'][1])],
                                'Pitch_Angle': data_dict['Pitch_Angle'],
                                'Energy': data_dict['Energy'],
                                'Epoch': data_dict['Epoch']}
            data_dict_output['Energy_Flux'][1]['LABLAXIS'] = 'Energy_Flux'
            data_dict_output['Energy_Flux'][1]['UNITS'] = 'erg cm!A-2!N'
            data_dict_output['Differential_Energy_Flux_avg'][1]['LABLAXIS'] = 'Differential_Energy_Flux_avg'

            data_dict_output = {**data_dict_output,**{
                'Energy_Flux_Downward':[EFlux_downward, deepcopy(data_dict_output['Energy_Flux'][1])]
            }}

            data_dict_output = {**data_dict_output, **{
                'Energy_Flux_Upward': [EFlux_Upward, deepcopy(data_dict_output['Energy_Flux'][1])]
            }}


            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

            outputCDFdata(outputPath, data_dict_output, globalAttrsMod=globalAttrs, instrNam=wInstr[1])

            Done(start_time)





# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        diffFlux_to_Energy_Flux(wRocket, rocketFolderPath, justPrintFileNames, wflyer, 0)
    elif wFiles == []:
        for i, wfile in enumerate(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')):
            diffFlux_to_Energy_Flux(wRocket, rocketFolderPath, justPrintFileNames, wflyer, wfile)
    else:
        for wfile in wFiles:
            diffFlux_to_Energy_Flux(wRocket, rocketFolderPath, justPrintFileNames, wflyer, wfile)
