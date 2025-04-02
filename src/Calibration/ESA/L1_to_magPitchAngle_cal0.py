# --- L1_to_magPitchAngle_cal0.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: For each rocket and each insrument, use the spun integration_tad_files data
# IN THE ROCKET FRAME as well as the orientation of the ESAs to determine a pitch angle
# for every time, pitch pad and energy of the instrument. Will use this to calibrate
# the ESAs.


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt

# --- --- --- --- ---
from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files. Both of these must be false to continue
justPrintFileNames = False
justPrintMAGFileNames = False

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

# select which files to convert
# [] --> all files, [#0,#1,#2,...etc] --> only specific files. Follows python indexing. use justPrintFileNames = True to see which files you need.
wFiles = [1, 4]
wMagFiles = [0] # use spun data!
inputPath_modifier = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_mag = 'l1'
outputPath_modifier = 'calibration\ESA_magPitch_calibration' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
modifier = ''

plot3DUnitVectors = False
plotMagVector = False

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.class_var_func import Rx, Ry, Rz


def L1_to_L1ESAmagCal(wFile, wMagFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L1'
    outputModelData = L1_TRICE_Quick(wflyer)

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\*.cdf')

    inputFiles_mag = glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wflyer]}{modifier}\*RingCore*')

    outputFiles = glob(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]
    input_names_mag = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wflyer]}\\', '') for ifile in inputFiles_mag]

    output_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\', '') for ofile in outputFiles]

    input_names_searchable = [ifile.replace('ACESII_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names]
    input_names_searchable_mag = [ifile.replace('ACESII_', '').replace(inputPath_modifier.lower() + '_','').replace('_v00', '') for ifile in input_names_mag]
    output_names_searchable = [ofile.replace('ACESII_', '').replace(outputPath_modifier.lower() +'_', '').replace('_v00', '').replace('__', '_') for ofile in output_names]

    dataFile_name = inputFiles[wFile].replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}{modifier}\\', '')

    # determine which instrument the file corresponds to:
    for index, instr in enumerate(['eepaa', 'leesa', 'iepaa', 'lp']):
        if instr in dataFile_name:
            wInstr = [index, instr]

    fileoutName = f'ACESII_{rocketID}_l1_{wInstr[1]}_magPitch.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            anws = ["yes" if input_names_searchable[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
    elif justPrintMAGFileNames:
        for i, file in enumerate(inputFiles_mag):
            anws = ["yes" if input_names_searchable_mag[i].replace('.cdf', "") in output_names_searchable else "no"]
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable_mag[i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Converting to {outputPath_modifier} data for {dataFile_name}' + color.END)
        print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

        # --- get the data from the l1 ESA file ---
        prgMsg(f'Loading data from {inputPath_modifier} Files')
        data_dict_esa = loadDictFromFile(inputFiles[wFile], {},reduceData=False,targetTimes=[],wKeys=[])
        data_dict_esa['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_esa['Epoch'][0][i]) for i in range(len(data_dict_esa['Epoch'][0]))])
        Done(start_time)

        # --- get the data from the l1 ESA file ---
        prgMsg(f'Loading data from Mag Files')
        data_dict_mag = loadDictFromFile(inputFiles_mag[wMagFile], {},reduceData=False,targetTimes=[],wKeys=[])
        data_dict_mag['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_mag['Epoch'][0][i]) for i in (range(len(data_dict_mag['Epoch'][0])))])
        Done(start_time)

        ################################################
        # --- Define EEPAA Unit Vectors in Mag Frame ---
        ################################################

        # Define the SPHERICAL description of the unit vectors IN THE MAG FRAME
        if wInstr[0] in [0, 2]:  # EEPAA and # IEPAA
            ThetaPolar = np.radians(90)
            unit_vect_temp = [[np.sin(ThetaPolar) * np.cos(np.radians(180-pitch)), np.sin(ThetaPolar) * np.sin(np.radians(180-pitch)), np.cos(ThetaPolar)] for pitch in data_dict_esa['Pitch_Angle'][0]]
            unit_vect = [np.matmul(Rx(45), temp_vect) for temp_vect in unit_vect_temp]

        elif wInstr[0] == 1:  # LEESA
            ThetaPolar = np.radians(90)
            unit_vect_temp = np.array([[np.sin(ThetaPolar) * np.cos(np.radians(180-ptch)), np.sin(ThetaPolar) * np.sin(np.radians(180-ptch)), np.cos(ThetaPolar)] for ptch in data_dict_esa['Pitch_Angle'][0]])
            unit_vect = [np.matmul(Rx(180), temp_vect) for temp_vect in unit_vect_temp]

        # 3D plot of vectors to verify their orientation
        if plot3DUnitVectors:
            ax = plt.figure().add_subplot(projection='3d')
            origin = [0, 0, 0]

            if wInstr[1] != 'iepaa':
                for i, vec in enumerate(unit_vect):
                    if i == 1:
                        ax.quiver(*origin, *vec, color='cyan',label='0deg')
                    elif i == 19:
                        ax.quiver(*origin, *vec, color='purple',label='180deg')
                    else:
                        ax.quiver(*origin, *vec)
            else:
                for i, vec in enumerate(unit_vect):
                    if i == 0:
                        ax.quiver(*origin, *vec, color='cyan', label='0deg')
                    elif i == len(data_dict_esa['Pitch_Angle'][0])-1:
                        ax.quiver(*origin, *vec, color='purple', label='180deg')
                    else:
                        ax.quiver(*origin, *vec)

            # define the X,Y,Z mag axis
            ax.quiver(*[0, 0, 0], *[1, 0, 0], color='red', label='x-axis')
            ax.quiver(*[0, 0, 0], *[0, 1, 0], color='blue', label='y-axis')
            ax.quiver(*[0, 0, 0], *[0, 0, 1], color='green', label='z-axis')

            ax.set_title(wInstr[1])
            ax.set_ylim(-1,1)
            ax.set_xlim(-1,1)
            ax.set_zlim(-1,1)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            plt.legend()
            plt.show()

        ########################
        # --- DOWNSAMPLE MAG ---
        ########################

        # --- organize Mag data into 3x1 column vectors ---
        reOrg_magData = np.array(
            [
                [data_dict_mag['Bx'][0][i] for i in range(len(data_dict_mag['Bx'][0]))],
                [data_dict_mag['By'][0][i] for i in range(len(data_dict_mag['By'][0]))],
                [data_dict_mag['Bz'][0][i] for i in range(len(data_dict_mag['Bz'][0]))]
            ]).transpose()

        # --- Downsample the Mag data ---
        prgMsg('Creating Mag fullDataSet data\n')
        esaLens = [len(data_dict_esa[wInstr[1]][0]), len(data_dict_esa[wInstr[1]][0][0]), len(data_dict_esa[wInstr[1]][0][0][0])]
        esaRanges = [range(esaLens[0]), range(esaLens[1]), range(esaLens[2])]

        epochFullDataSet = np.zeros(shape=(esaLens[0], esaLens[2]))
        magFullDataSet = np.zeros(shape=(esaLens[0], esaLens[2], 3))

        # --- populate epochFullDataSet ---
        # NOTE: the eepaa only has 41 energy bins since the first 8 were removed. Take this into account

        energyAdjust = 8 if wInstr[1] =='eepaa' else 0

        # since the Epoch is the timestamp of the lowest energy step (7.25eV)
        # the time of the step should be: Timestamp + (#ofSteps - 1 + account for missing 8 steps if needed - index of next step)
        # this format has the highest energy in first and decreases in time to the lowest step

        for tme, engy in tqdm(itertools.product(*[esaRanges[0], esaRanges[2]])):
            epochFullDataSet[tme][engy] = data_dict_esa['Epoch'][0][tme] + (esaLens[2] - 1 + energyAdjust - engy) * rocketAttrs.timeBetweenSteps_in_ns

        # --- assign a integration_tad_files Bx, By, Bz to each tme and engy ---
        # NOTE: don't need to do it for ptch since they're all sampled at the same time
        for tme in tqdm(esaRanges[0]):

            # find the closest value between the B-epoch and the assigned epoch value for the ESA data. Data is sorted from highest energy value to lowest
            B_vectors = [reOrg_magData[np.abs(data_dict_mag['Epoch'][0] - epochFullDataSet[tme][engy]).argmin()] for engy in esaRanges[2] ]

            # normalize all the B_vectors
            magFullDataSet[tme] = np.array(list([vector/np.linalg.norm(vector) for vector in B_vectors]))

        # magData should look like: [ [[x,y,z],[x,y,z],...], [[x,y,z],[x,y,z],...], ...   ]

        Done(start_time)

        # --- 3D plot of vectors to verify their orientation
        wTimes = np.linspace(100, len(magFullDataSet), 30)

        if plotMagVector:

            for j in range(len(wTimes)):

                bvec = magFullDataSet[int(wTimes[j])]

                for k, val in enumerate(bvec):

                    if k == 0:

                        ax = plt.figure().add_subplot(projection='3d')
                        origin = [0, 0, 0]
                        zeroDegvec= origin

                        for i, vec in enumerate(unit_vect):
                            if i == 1:
                                ax.quiver(*origin, *vec, color='cyan', label='$0^{\circ}$')
                                zeroDegvec = vec
                            elif i == 19:
                                ax.quiver(*origin, *vec, color='purple', label='$180^{\circ}$')
                            else:
                                ax.quiver(*origin, *vec)

                        # define the X,Y,Z mag axis
                        ax.quiver(*origin, *[1, 0, 0], color='red',label='x-axis')
                        ax.quiver(*origin, *[0, 1, 0], color='blue',label='y-axis')
                        ax.quiver(*origin, *[0, 0, 1], color='green',label='z-axis')

                        # angle between 0deg and B
                        zeroDegDot = np.degrees(np.arccos(np.dot(val,zeroDegvec)))

                        ax.quiver(*[0, 0, 0], *val, color='lightpink', label=fr'B, $n_{0} \dot B: {zeroDegDot}')

                        ax.set_title(wInstr[1] + f'Time Index {int(wTimes[j])}, Energy Index {k}')
                        ax.set_ylim(-1, 1)
                        ax.set_xlim(-1, 1)
                        ax.set_zlim(-1, 1)
                        ax.set_xlabel('x')
                        ax.set_ylabel('y')
                        ax.set_zlabel('z')
                        plt.legend()
                        plt.show()

        ####################################
        # --- Calculate TRUE pitch angle ---
        ####################################

        pitch_Angle_magCalc = np.zeros(shape=(len(data_dict_esa[wInstr[1]][0]), len(data_dict_esa[wInstr[1]][0][0]), len(data_dict_esa[wInstr[1]][0][0][0])))

        for tme, ptch, engy in tqdm(itertools.product(*esaRanges)):
            pitch_Angle_magCalc[tme][ptch][engy] = np.degrees(
                np.arccos(
                    np.dot(unit_vect[ptch], magFullDataSet[tme][engy])/np.linalg.norm(magFullDataSet[tme][engy])
                )
            )

        # --- --- --- --- --- --- --- ---
        # --- PREPARE DATA FOR OUTPUT ---
        # --- --- --- --- --- --- --- ---
        data_dict = {}
        data_dict = {**data_dict, **{'Mag_Calculated_Pitch_Angle':
                                         [pitch_Angle_magCalc, {'LABLAXIS': 'Mag_Calculated_Pitch_Angle',
                                                      'DEPEND_0': 'Epoch', 'DEPEND_1': 'Pitch_Angle',
                                                      'DEPEND_2': 'Energy',
                                                      'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                      'UNITS': 'deg',
                                                      'VALIDMIN': pitch_Angle_magCalc.min(), 'VALIDMAX': pitch_Angle_magCalc.max(),
                                                      'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}
        # get Energy and Pitch angle into data_dict
        for key, val in data_dict_esa.items():
            if key in ['Energy', 'Pitch_Angle', 'Epoch']:
                data_dict = {**data_dict, **{key:data_dict_esa[key]}}

        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---
        prgMsg('Creating output file')

        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'
        globalAttrsMod['Descriptor'] = rocketAttrs.InstrNames_Full[wInstr[0]]
        outputCDFdata(outputPath, data_dict, outputModelData, globalAttrsMod, wInstr[1])

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
        L1_to_L1ESAmagCal(0, wMagFiles[0], rocketFolderPath, justPrintFileNames, wflyer)
    elif wFiles == []:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wflyer]}\*.cdf')))):
            L1_to_L1ESAmagCal(fileNo, wMagFiles[0], rocketFolderPath, justPrintFileNames, wflyer)
    else:
        for filesNo in wFiles:
            L1_to_L1ESAmagCal(filesNo, wMagFiles[0], rocketFolderPath, justPrintFileNames, wflyer)