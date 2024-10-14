# --- EEPAA_LEESA_superset.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Combine the High Flyer EEPAA and LEESA data into a single, large electron dataset



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False
# ---------------------------
outputData = True
# ---------------------------

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
# none


def EEPAA_LEESA_superset(wRocket, rocketFolderPath):
    wflyer = wRocket-4

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    ModelData = L2_TRICE_Quick(wflyer)

    inputFile_EEPAA_Flux = glob(f'{rocketFolderPath}L2\{fliers[wflyer]}\*eepaa_fullCal*')[0]
    inputFile_LEESA_Flux = glob(f'{rocketFolderPath}L2\{fliers[wflyer]}\*leesa_fullCal*')[0]

    inputFile_EEPAA_DistFunc = glob(f'{rocketFolderPath}L3\distFunc\{fliers[wflyer]}\*eepaa_fullCal*')[0]
    inputFile_LEESA_DistFunc = glob(f'{rocketFolderPath}L3\distFunc\{fliers[wflyer]}\*leesa_fullCal*')[0]

    fileoutName_flux = f'ACESII_{rocketID}_l2_eSuperSet.cdf'
    fileoutName_Dist = f'ACESII_{rocketID}_distFunc_eSuperSet.cdf'

    # --- --- --- --- --- -
    # --- LOAD THE DATA ---
    # --- --- --- --- --- -

    # --- get the data from the file ---
    prgMsg(f'Loading Electron Files Files')
    # load data here
    data_dict_EEPAA_flux = loadDictFromFile(inputFile_EEPAA_Flux)
    data_dict_LEESA_flux = loadDictFromFile(inputFile_LEESA_Flux)
    data_dict_EEPAA_distFunc = loadDictFromFile(inputFile_EEPAA_DistFunc)
    data_dict_LEESA_distFunc = loadDictFromFile(inputFile_LEESA_DistFunc)
    Done(start_time)

    # match the time bases:
    prgMsg('Matching Time Bases')
    # the EEPAA dataset is longer than the LEESA dataset. Reduce EEPAA to match LEESA
    reducedEpochIndicies = np.array([np.abs(data_dict_EEPAA_flux['Epoch'][0] - data_dict_LEESA_flux['Epoch'][0][i]).argmin() for i in range(len(data_dict_LEESA_flux['Epoch'][0]))])
    data_dict_EEPAA_flux['Differential_Energy_Flux'][0] = data_dict_EEPAA_flux['Differential_Energy_Flux'][0][reducedEpochIndicies]
    data_dict_EEPAA_flux['Differential_Number_Flux'][0] = data_dict_EEPAA_flux['Differential_Number_Flux'][0][reducedEpochIndicies]
    data_dict_EEPAA_flux['oneCountLevel'][0] = data_dict_EEPAA_flux['oneCountLevel'][0][reducedEpochIndicies]
    data_dict_EEPAA_distFunc['Distribution_Function'][0] = data_dict_EEPAA_distFunc['Distribution_Function'][0][reducedEpochIndicies]
    Done(start_time)

    newEnergy = np.append(data_dict_EEPAA_flux['Energy'][0],data_dict_LEESA_flux['Energy'][0])
    newEFlux = np.zeros(shape=(len(data_dict_EEPAA_flux['Differential_Energy_Flux'][0]),len(data_dict_EEPAA_flux['Differential_Energy_Flux'][0][0]),len(newEnergy)))
    newNFlux = np.zeros(shape=(len(data_dict_EEPAA_flux['Differential_Energy_Flux'][0]), len(data_dict_EEPAA_flux['Differential_Energy_Flux'][0][0]), len(newEnergy)))
    newDistFunc = np.zeros(shape=(len(data_dict_EEPAA_flux['Differential_Energy_Flux'][0]), len(data_dict_EEPAA_flux['Differential_Energy_Flux'][0][0]), len(newEnergy)))
    newoneCountLevel = np.zeros(shape=(len(data_dict_EEPAA_flux['oneCountLevel'][0]), len(data_dict_EEPAA_flux['oneCountLevel'][0][0]), len(newEnergy)))

    from itertools import product
    for tme, ptch in tqdm(product(*[range(len(data_dict_EEPAA_flux['Differential_Energy_Flux'][0])),range(len(data_dict_EEPAA_flux['Differential_Energy_Flux'][0][0]))])):
        newEFlux[tme][ptch] = np.append(data_dict_EEPAA_flux['Differential_Energy_Flux'][0][tme][ptch], data_dict_LEESA_flux['Differential_Energy_Flux'][0][tme][ptch])
        newNFlux[tme][ptch] = np.append(data_dict_EEPAA_flux['Differential_Number_Flux'][0][tme][ptch], data_dict_LEESA_flux['Differential_Number_Flux'][0][tme][ptch])
        newDistFunc[tme][ptch] = np.append(data_dict_EEPAA_distFunc['Distribution_Function'][0][tme][ptch], data_dict_LEESA_distFunc['Distribution_Function'][0][tme][ptch])
        newoneCountLevel[tme][ptch] = np.append(data_dict_EEPAA_distFunc['oneCountLevel'][0][tme][ptch], data_dict_LEESA_distFunc['oneCountLevel'][0][tme][ptch])

    # create a new flux data_dict
    newFluxVals = {'Differential_Energy_Flux': newEFlux,
                   'Differential_Number_Flux': newNFlux,
                   'Energy': newEnergy}
    data_dict_newFlux = deepcopy(data_dict_LEESA_flux)
    for key, val in newFluxVals.items():
        data_dict_newFlux[key][0] = val

    # create a new distribution data_dict
    newDistVals = {
        'Distribution_Function': newDistFunc,
        'oneCountLevel': newoneCountLevel,
        'Energy': newEnergy}

    data_dict_newDist = deepcopy(data_dict_LEESA_distFunc)

    for key, val in newDistVals.items():
        data_dict_newDist[key][0] = newDistVals[key]

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        prgMsg('Creating output file')

        outputPath_Flux = f'{rocketFolderPath}L2\{fliers[0]}\\{fileoutName_flux}'
        outputPath_Dist = f'{rocketFolderPath}L3\DistFunc\{fliers[0]}\\{fileoutName_Dist}'

        outputCDFdata(outputPath_Flux, data_dict_newFlux, ModelData= ModelData, globalAttrsMod= globalAttrsMod, instrNam= 'electronSuperSet')
        outputCDFdata(outputPath_Dist, data_dict_newDist,ModelData= ModelData, globalAttrsMod= globalAttrsMod, instrNam= 'electronSuperSet')

        Done(start_time)


# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
EEPAA_LEESA_superset(4, rocketFolderPath)
