# --- template.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Determine the regions during flight that the payloads were on
# a similar magnetic footprint. The attitude data is loaded in so I can plot things



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

# just print
justPrintFileNames = False # just prints the attitude data
wRocketPrint = 0 #0 - high flyer, 1- low flyer

# Just print the names of files of the magnetic field
justPrintMagFileNames = False
wRocketMagPrint = 1 # which magnetic field data to print. 0 - high flyer, 1- low flyer
wBFile = [0, 0]  # choice of magnetic field data to use [high flyer, low flyer]

modifier = ''
inputPath_modifier = 'attitude'
inputPath_modifier_mag = 'l2'
outputPath_modifier = 'science'

# angle the two magnetic fields need to be within to be "aligned"
angleThresh = 5

outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---


def magneticFootprintAlignment(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)


    inputFiles = [glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[0]}{modifier}\*.cdf'),glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[1]}{modifier}\*.cdf')]
    inputFiles_mag = [glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[0]}{modifier}\*Ringcore*'), glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[1]}{modifier}\*Ringcore*')]
    outputFiles = glob(f'{rocketFolderPath}{outputPath_modifier}\*MagneticAlignment*')

    input_names = [[ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[0]}{modifier}\\', '') for ifile in inputFiles[0]],[ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{fliers[1]}{modifier}\\', '') for ifile in inputFiles[1]]]
    input_names_mag = [[ifile.replace(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[0]}{modifier}\\', '') for ifile in inputFiles_mag[0]],[ifile.replace(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[1]}{modifier}\\', '') for ifile in inputFiles_mag[1]]]
    output_names = [ofile.replace(f'{rocketFolderPath}{outputPath_modifier}\\', '') for ofile in outputFiles]

    input_names_searchable = [[ifile.replace('ACES_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names[0]],[ifile.replace('ACES_', '').replace(inputPath_modifier.lower() +'_', '').replace('_v00', '') for ifile in input_names[1]]]
    input_names_searchable_mag = [[ifile.replace('ACES_', '').replace(inputPath_modifier_mag.lower() + '_', '').replace('_v00', '') for ifile in input_names_mag[0]],[ifile.replace('ACES_', '').replace(inputPath_modifier_mag.lower() + '_', '').replace('_v00', '') for ifile in input_names_mag[1]]]

    fileoutName = f'ACESII_MagneticAlignment.cdf'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles[wRocketPrint]):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[wRocketPrint][i], round(getsize(file) / (10 ** 6), 1)))
    elif justPrintMagFileNames:
        for i, file in enumerate(inputFiles_mag[wRocketMagPrint]):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable_mag[wRocketMagPrint][i], round(getsize(file) / (10 ** 6), 1)))
    else:
        print('\n')
        print(color.UNDERLINE + f'Aligning Magnetic Fields' + color.END)

        # --- get the data from the ACS file ---
        prgMsg(f'Loading Attitude and Magnetometer Data')
        attitudeDicts = []
        magDicts = []

        for i in range(2):

            data_dict_attitude = {}
            with pycdf.CDF(inputFiles[i][0]) as inputDataFile:
                for key, val in inputDataFile.items():
                    data_dict_attitude = {**data_dict_attitude, **{key : [inputDataFile[key][...], {key:val for key,val in inputDataFile[key].attrs.items()  }  ]  }  }

            data_dict_attitude['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_attitude['Epoch'][0][i]) for i in (range(len(data_dict_attitude['Epoch'][0])))])
            attitudeDicts.append(data_dict_attitude)

            # --- get the data from the mag file ---
            data_dict_mag = {}
            with pycdf.CDF(inputFiles_mag[i][wBFile[i]]) as inputDataFile:
                for key, val in inputDataFile.items():
                    data_dict_mag = {**data_dict_mag, **{key: [inputDataFile[key][...], {key: val for key, val in inputDataFile[key].attrs.items()}]}}

            data_dict_mag['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(data_dict_mag['Epoch'][0][i]) for i in (range(len(data_dict_mag['Epoch'][0])))])
            magDicts.append(data_dict_mag)

        Done(start_time)

        ######################################
        # --- Determine Magnetic Alignment ---
        ######################################

        # find where low flyer t=0,t=end occurs in high flyer mag data
        LF_epoch_start_index = np.abs(magDicts[1]['Epoch'][0] - magDicts[1]['Epoch'][0][0]).argmin()
        LF_epoch_end_index = np.abs(magDicts[1]['Epoch'][0] - magDicts[1]['Epoch'][0][-1]).argmin() + 1

        # reduce HF mag data
        if wBFile[0] == 0 or wBFile[1] == 0:
            Blabels = ['B_east', 'B_north', 'B_up']
        else:
            Blabels = ['Bx', 'By', 'Bz']

        B_Field_HF = np.array([magDicts[0][Blabels[0]][0][LF_epoch_start_index:LF_epoch_end_index], magDicts[0][Blabels[1]][0][LF_epoch_start_index:LF_epoch_end_index], magDicts[0][Blabels[2]][0][LF_epoch_start_index:LF_epoch_end_index]]).transpose()
        B_Field_LF = np.array([magDicts[1][Blabels[0]][0], magDicts[1][Blabels[1]][0], magDicts[1][Blabels[2]][0]]).transpose()
        Epoch_LF = np.array(magDicts[1]['Epoch'][0])

        magAlignment = []

        prgMsg('Calculating MagneticAlignment')
        # Calculate the magAlignment
        for i in range(len(Epoch_LF)):
            magAlignment.append(np.degrees(np.arccos(np.dot(B_Field_LF[i], B_Field_HF[i])/ (np.linalg.norm(B_Field_LF[i])*np.linalg.norm(B_Field_HF[i])))))

        magAlignment = np.array(magAlignment)

        Done(start_time)
        # --- --- --- --- --- --- ---
        # --- WRITE OUT THE DATA ---
        # --- --- --- --- --- --- ---

        if outputData:

            prgMsg('Creating output file')

            data_dict = {}
            data_dict = {**data_dict, **{f'Epoch':
                                             [Epoch_LF, {'LABLAXIS': f'Epoch',
                                                     'DEPEND_0': 'Epoch',
                                                     'DEPEND_1': None,
                                                     'DEPEND_2': None,
                                                     'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                     'UNITS': None,
                                                     'VALIDMIN': Epoch_LF.min(),
                                                     'VALIDMAX': Epoch_LF.max(),
                                                     'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]}}

            data_dict = {**data_dict, **{'MagneticAlignment':
                                             [magAlignment, {'LABLAXIS': f'magAlignment',
                                                     'DEPEND_0': 'Epoch',
                                                     'DEPEND_1': None,
                                                     'DEPEND_2': None,
                                                     'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2',
                                                     'UNITS': 'deg',
                                                     'VALIDMIN': magAlignment.min(),
                                                     'VALIDMAX': magAlignment.max(),
                                                     'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]}}




            outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fileoutName}'
            outputCDFdata(outputPath, data_dict, outputModelData,globalAttrsMod, 'MagneticAlignment')

        Done(start_time)




# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
magneticFootprintAlignment(0, 0, rocketFolderPath, justPrintFileNames, 0)
