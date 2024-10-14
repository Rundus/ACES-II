# --- mSSA_to_calcComponents.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:  The output are SSA component files found in \L3\SSAcomponents_B
# it this code must be able to do the following SIMPLY
# [0] handle Electric or magnetic data
# [1] perform mSSA on an input data_dict
# [2] be able to output components for the whole time series, a subset or the whole series broken into subsets
import math

import numpy as np

# --- bookkeeping ---
from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

# --- Select the DataSet ---
inputPath_modifier = 'science\PoyntingFlux'
# inputPath_modifier = 'L2'
wFiles = [0]


#################
# --- TOGGLES ---
#################

# --- Mirror Data ---
# for calculating the components, the data is mirrored by "mirrorPercentage" amount on either end
MirrorData = True
mirrorPercentage = 0.2

# --- Window Size ---
SSA_window_Size = 501
# SSA_window_Size = 5

# ---- SSA Components ----
computeSSAcomponents = True
breakIntoThisManyFiles = 5
generate_only_these_files = [] # if I only want to re-generate specific component subsets
# ------------------------
outputData = True
# ------------------------


def mSSA_to_calcComponents(wRocket, rocketFolderPath, justPrintFileNames,wFile):

    # --- ACES-II Attributes Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wRocket-4]
    inputFiles, input_names, input_names_searchable = getInputFiles(rocketFolderPath=rocketFolderPath,wRocket=wRocket,inputPath_modifier=inputPath_modifier)

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return


    # --- Get the Input Data ---
    prgMsg('Loading Data for ' + input_names_searchable[wFile])
    data_dict, globalAttrs = loadDictFromFile(inputFiles[wFile], getGlobalAttrs=True)
    Done(start_time)

    # --- Determine All the information about the input file ---
    # component names
    compNames, coordSys, coordSet = getCoordinateKeys(data_dict)

    # Where the componentSSA files need to go
    varName = compNames[0].replace(coordSet[0], '')
    pathModifier = 'SSAcomponents' + f'_{varName}'

    wInstr = ''
    for vari in ['RingCore','PoyntingFlux','mBeEr','BrEe','E_Field']:
        if vari in inputFiles[wFile]:
            wInstr = vari
            break

    if computeSSAcomponents:
        ##############################
        # --- CALC mSSA Components ---
        ##############################
        # --- break data into segments based on "BreakIntoThisManyFiles" variable ---
        # find the inidcies to break up the dataset into "breakIntoThisManyPieces" parts
        rangeLimits = [[ar[0], ar[-1]+1] for ar in np.array_split([i for i in range(len(data_dict['Epoch'][0]))], breakIntoThisManyFiles)]

        # modify rangeLimits to only include those in "generate_these_files"
        if generate_only_these_files != []:
            rangeLimits = [rangeLimits[num] for num in generate_only_these_files]


        # --- DO THE mSSA component Calculation ---
        from ACESII_code.class_var_func import mSSA_components

        for i, subset in enumerate(rangeLimits):

            prgMsg(f'Calculating SSA components for Indicies {subset[0]} to {subset[1]}')

            # get the data_dict subset
            temp_dict = deepcopy(data_dict)
            data_dict_chunk = {key: [val[0][subset[0]:subset[1]], val[1]] for key, val in temp_dict.items()}

            data_dict_output = mSSA_components(data_dict_input=data_dict_chunk, compNames=compNames,
                                               SSA_window_Size=SSA_window_Size, mirrorData=MirrorData,
                                               mirrorPercent=mirrorPercentage)

            # --- --- --- --- ---
            # --- OUTPUT DATA ---
            # --- --- --- --- ---
            if outputData:

                outputFilePath = rf'{rocketFolderPath}L3\\{pathModifier}\\{fliers[wRocket-4]}\\ACESII_{rocketID}_{wInstr}_SSAComponents_{coordSys}_WL{SSA_window_Size}_subset_{i}'
                if MirrorData:
                    outputFilePath = outputFilePath + f'_mirrored{int(mirrorPercentage*100)}'
                outputCDFdata(outputFilePath+".cdf", data_dict_output, globalAttrsMod=globalAttrs, instrNam=wInstr)
                Done(start_time)





# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        mSSA_to_calcComponents(wRocket, rocketFolderPath, justPrintFileNames, 0)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}\*.cdf')))):
            mSSA_to_calcComponents(wRocket, rocketFolderPath, justPrintFileNames, fileNo)
    else:
        for filesNo in wFiles:
            mSSA_to_calcComponents(wRocket, rocketFolderPath, justPrintFileNames, filesNo)

