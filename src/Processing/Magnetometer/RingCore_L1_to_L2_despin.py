# # --- RingCore_L1_to_L2_despin.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: Just interpolate the attitude solution data and apply the wallops DCM
# to the integration_tad_files data, getting it in ENU coordinates. Apply it over the entire journey of the flight



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
from src.my_imports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5
wFiles = [0]

modifier = ''
inputPath_modifier = 'l1' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
inputPath_modifier_attitude = 'attitude' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
outputPath_modifier_despin = 'l2' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder

include_CHAOS_model = False
find_time_offset = False
find_DC_offset = False
replaceNANS = True
outputData = True # plots RKT XYZ and CHAOS model Spun-up XYZ

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from scipy.interpolate import CubicSpline
from src.Processing.Magnetometer.RingCore_despin_toggles import DespinToggles

def RingCore_L1_to_L2_Despin(wRocket, wFile, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---

    rocketID = ACESII.payload_IDs[wflyer]
    globalAttrsMod = ACESII.global_attributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wflyer]}{modifier}\*RingCore_rktFrm*')
    inputFiles_attitude = glob(f'{rocketFolderPath}{inputPath_modifier_attitude}\{ACESII.fliers[wflyer]}{modifier}\*.cdf')

    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wflyer]}{modifier}\\', '') for ifile in inputFiles]

    input_names_searchable = [ifile.replace('_v00', '') for ifile in input_names]

    fileoutName_despin = f'ACESII_{rocketID}_l2_RingCore_ENU'

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return

    print('\n')
    print(stl.color.UNDERLINE + f'DeSpining RingCore Data' + stl.color.END)
    print('[' + str(wFile) + ']   ' + str(round(getsize(inputFiles[wFile]) / (10 ** 6), 1)) + 'MiB')

    # --- get the data from the Magnetometer file ---
    stl.prgMsg(f'Loading data from {inputPath_modifier} RingCore Files')
    data_dict_mag = stl.loadDictFromFile(inputFiles[wFile],targetVar=[DespinToggles.reduceTimes[wRocket-4],'Epoch'])
    stl.Done(start_time)

    # --- get the data from the attitude file ---
    stl.prgMsg(f'Loading data from {inputPath_modifier_attitude} Files')
    data_dict_attitude = stl.loadDictFromFile(inputFiles_attitude[0],targetVar=[DespinToggles.reduceTimes[wRocket-4],'Epoch'])
    data_dict_attitude['Alt'][0] = data_dict_attitude['Alt'][0]/1000
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = deepcopy(data_dict_mag)

    ########################
    # --- Reduce dataset ---
    ########################
    # --- prepare some variables for later ---
    # create B_rkt variable
    B_rkt = np.array([[data_dict_mag['Bx'][0][i], data_dict_mag['By'][0][i], data_dict_mag['Bz'][0][i]] for i in range(len(data_dict_mag['Epoch'][0]))])

    # convert attitude epoch to tt2000
    Epoch_attitude_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_attitude['Epoch'][0]])
    Epoch_mag_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_mag['Epoch'][0]])

    #################################################
    # --- Determine CHAOS model on Attitude Epoch ---
    #################################################
    stl.prgMsg('Calculating CHAOS model on Attitude Epoch')
    B_CHAOS_ENU_attitude = stl.CHAOS(
                    lat=data_dict_attitude['Lat'][0],
                    long=data_dict_attitude['Long'][0],
                    alt=data_dict_attitude['Alt'][0],
                    times=data_dict_attitude['Epoch'][0])
    stl.Done(start_time)

    ############################################
    # --- Spin Up CHAOS model into RKT frame ---
    ############################################
    stl.prgMsg('Convert CHAOS to rkt frame')
    # --- get the DCM ---
    DCM = np.array(
        [
            [[data_dict_attitude['a11'][0][i], data_dict_attitude['a12'][0][i], data_dict_attitude['a13'][0][i]],
             [data_dict_attitude['a21'][0][i], data_dict_attitude['a22'][0][i], data_dict_attitude['a23'][0][i]],
             [data_dict_attitude['a31'][0][i], data_dict_attitude['a32'][0][i], data_dict_attitude['a33'][0][i]]]
            for i in range(len(data_dict_attitude['Epoch'][0]))
        ]
    )

    DCMinv = np.array([np.linalg.inv(mat) for mat in DCM])
    B_CHAOS_rkt = np.array([np.matmul(DCMinv[i], B_CHAOS_ENU_attitude[i]) for i in range(len(data_dict_attitude['Epoch'][0]))])
    stl.Done(start_time)

    ##################################################
    # --- Find Time offset between CHAOS and B_rkt ---
    ##################################################
    if find_time_offset:

        # reduce the datasets to a small region where you want these calibrations to be right
        calibrationRegion = [dt.datetime(2022, 11, 20, 17, 22, 30, 00), dt.datetime(2022, 11, 20, 17, 26, 00, 00)]

        # attitude
        lowTime, highTime = np.abs(data_dict_attitude['Epoch'][0] - calibrationRegion[0]).argmin(), np.abs(data_dict_attitude['Epoch'][0] - calibrationRegion[1]).argmin()
        for key,val in data_dict_attitude.items():
            data_dict_attitude[key][0] = data_dict_attitude[key][0][lowTime:highTime]

        B_CHAOS_rkt = B_CHAOS_rkt[lowTime:highTime]
        Epoch_attitude_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_attitude['Epoch'][0]])

        # mag
        lowTime, highTime = np.abs(data_dict_mag['Epoch'][0] - calibrationRegion[0]).argmin(), np.abs(data_dict_mag['Epoch'][0] - calibrationRegion[1]).argmin()
        for key,val in data_dict_mag.items():
            data_dict_mag[key][0]=data_dict_mag[key][0][lowTime:highTime]
        B_rkt = np.array([[data_dict_mag['Bx'][0][i], data_dict_mag['By'][0][i], data_dict_mag['Bz'][0][i]] for i in range(len(data_dict_mag['Epoch'][0]))])
        Epoch_mag_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_mag['Epoch'][0]])

        fig, ax = plt.subplots(3)
        for i in range(3):
            ax[i].plot(data_dict_attitude['Epoch'][0], B_CHAOS_rkt[:, i], label='B_CHAOS_attitude')
            ax[i].plot(data_dict_mag['Epoch'][0], B_rkt[:, i], label='B_rkt')
        plt.show()

        stl.prgMsg('Finding best temporal offset')

        N = 30
        timeOffsets = np.linspace(0.12755, 0.12757, N) * 1E9

        bestChi = [1, 1E30]

        for deltaT in timeOffsets:
            newAttitudeEpoch = np.array([ int(tme + deltaT) for tme in Epoch_attitude_tt2000])
            newMag = []

            for i in range(3):
                # --- cubic interpolation ---
                splCub = CubicSpline(newAttitudeEpoch, B_CHAOS_rkt[:, i])

                # --- evaluate the interpolation at all the epoch_mag points ---
                newMag.append(np.array([splCub(timeVal) for timeVal in Epoch_mag_tt2000]))

            # create B_CHAOS_RKT variable that's on mag Epoch
            B_CHAOS_new = np.array([[newMag[0][i], newMag[1][i], newMag[2][i]] for i in range(len(data_dict_mag['Epoch'][0]))])

            # determine if its a good choice

            ChiSquare = []

            for i in range(len(B_CHAOS_new)):
                B1 = (B_CHAOS_new[i][0] - B_rkt[i][0])**2
                B2 = (B_CHAOS_new[i][1] - B_rkt[i][1])**2
                B3 = (B_CHAOS_new[i][2] - B_rkt[i][2])**2

                if np.isnan(B1) or np.isnan(B2) or np.isnan(B3):
                    ChiSquare.append((0) / len(B_CHAOS_new))
                else:
                    ChiSquare.append((B1 + B2 + B3) / len(B_CHAOS_new))

            ChiSquare = sum(ChiSquare)

            print(deltaT/1E9, ChiSquare)
            if ChiSquare < bestChi[1]:
                bestChi = [deltaT, ChiSquare]

        print('The best offset was:', bestChi)
        stl.Done(start_time)
    else:  # apply the temporal offset
        TimeOffset = [127567241.37931032, 120789473.68421052]

    Epoch_interpolateThis = np.array([int(tme + TimeOffset[wRocket-4]) for tme in Epoch_attitude_tt2000])

    ##########################################
    # --- interpolate CHAOS onto mag Epoch ---
    ##########################################
    stl.prgMsg('Interpolating CHAOS onto integration_tad_files Epoch')
    B_CHAOS_rkt_magTime = []

    for i in range(3):

        # --- cubic interpolation ---
        splCub = CubicSpline(Epoch_interpolateThis, B_CHAOS_rkt[:, i])

        # --- evaluate the interpolation at all the epoch_mag points ---
        B_CHAOS_rkt_magTime.append(np.array([splCub(timeVal) for timeVal in Epoch_mag_tt2000]))

    # create B_CHAOS_RKT variable that's on mag Epoch
    B_CHAOS_rkt_magTime = np.array([[B_CHAOS_rkt_magTime[0][i], B_CHAOS_rkt_magTime[1][i], B_CHAOS_rkt_magTime[2][i]] for i in range(len(Epoch_mag_tt2000))])
    stl.Done(start_time)


    ###################################################
    # --- Find DC offset between CHAOS and B_rkt ---
    ###################################################

    if find_DC_offset:
        # the offset primarily exists for the B_X variable

        fig, ax = plt.subplots(3)
        for i in range(3):
            ax[i].plot(data_dict_mag['Epoch'][0], B_CHAOS_rkt_magTime[:, i], label='B_CHAOS_attitude')
            ax[i].plot(data_dict_mag['Epoch'][0], B_rkt[:, i], label='B_rkt')
        plt.show()

        N = 10000
        magOffsets = np.linspace(-100, 100, N)
        bestChi = [1, 1E30]

        # --- B_DC_offset ---
        wComp = 2
        for deltaB in magOffsets:
            newB = np.array(B_CHAOS_rkt_magTime[:,wComp] + deltaB)

            ChiSquare = (1/len(B_rkt))*sum(np.array(newB - B_rkt[:,wComp])**2)
            # print(deltaB, ChiSquare)
            if ChiSquare < bestChi[1]:
                bestChi = [deltaB,ChiSquare]
        print('Best ChiSquare B: ', bestChi)
    else:  # apply the temporal offset
        # offsetDC = np.array([[5.714285714285714, -0.5306122448979592, 0.346938775510204], [175.78947368421052, 3.183673469387755, -183.6326530612245]])
        offsetDC = np.array([[0.10010010010010717,-1.2901290129012892,4.590459045904595], [0,0,0]])

    B_CHAOS_rkt_magTime = B_CHAOS_rkt_magTime + offsetDC[wRocket-4]

    #####################################
    # --- interpolate DCM to mag time ---
    #####################################

    # --- GET THE DCM FOR THE RKT with No time adjustments ---
    data_dict_attitude_temp = {}

    for key in ['a11', 'a12', 'a13', 'a21', 'a22', 'a23', 'a31', 'a32', 'a33']:

        # --- cubic interpolation ---
        # splCub = CubicSpline(Epoch_attitude_tt2000, data_dict_attitude[key][0])
        splCub = CubicSpline(Epoch_attitude_tt2000, data_dict_attitude[key][0])

        # --- evaluate the interpolation at all the epoch_mag points ---
        data_dict_attitude_temp = {**data_dict_attitude_temp, **{key:np.array([splCub(timeVal) for timeVal in Epoch_mag_tt2000])}}

    DCM_magEpoch_forRKT = np.array(
        [[ [data_dict_attitude_temp['a11'][i], data_dict_attitude_temp['a12'][i], data_dict_attitude_temp['a13'][i]],
           [data_dict_attitude_temp['a21'][i], data_dict_attitude_temp['a22'][i], data_dict_attitude_temp['a23'][i]],
           [data_dict_attitude_temp['a31'][i], data_dict_attitude_temp['a32'][i], data_dict_attitude_temp['a33'][i]]]
            for i in range(len(data_dict_mag['Epoch'][0]))
        ]
    )

    # apply DCM to modified CHAOS data
    B_rkt_ENU = np.array([np.matmul(DCM_magEpoch_forRKT[i], B_rkt[i]) for i in range(len(Epoch_mag_tt2000))])

    # --- GET THE DCM FOR CHAOS WITH time adjustments ---
    data_dict_attitude_temp = {}

    for key in ['a11', 'a12', 'a13', 'a21', 'a22', 'a23', 'a31', 'a32', 'a33']:
        # --- cubic interpolation ---
        splCub = CubicSpline(Epoch_interpolateThis, data_dict_attitude[key][0])

        # --- evaluate the interpolation at all the epoch_mag points ---
        data_dict_attitude_temp = {**data_dict_attitude_temp, **{key: np.array([splCub(timeVal) for timeVal in Epoch_mag_tt2000])}}

    DCM_magEpoch_forCHAOS = np.array(
        [[[data_dict_attitude_temp['a11'][i], data_dict_attitude_temp['a12'][i], data_dict_attitude_temp['a13'][i]],
          [data_dict_attitude_temp['a21'][i], data_dict_attitude_temp['a22'][i], data_dict_attitude_temp['a23'][i]],
          [data_dict_attitude_temp['a31'][i], data_dict_attitude_temp['a32'][i], data_dict_attitude_temp['a33'][i]]]
         for i in range(len(data_dict_mag['Epoch'][0]))
         ]
    )

    # apply DCM to CHAOS data
    B_CHAOS_ENU_magTime = np.array([np.matmul(DCM_magEpoch_forCHAOS[i], B_CHAOS_rkt_magTime[i]) for i in range(len(Epoch_mag_tt2000))])

    #####################################
    # --- Subtract B_Model from B_Rkt ---
    #####################################
    DeltaB_ENU = B_rkt_ENU - B_CHAOS_ENU_magTime

    ###################################
    # --- Handle NAN values in data ---
    ###################################
    if replaceNANS:

        # # find the Nans
        # nanIndices = []
        # for i in range(len(DeltaB_ENU)):
        #     if np.isnan(DeltaB_ENU[i][0]):
        #         nanIndices.append(i)

        def nan_helper(y):
            return np.isnan(y), lambda  z: z.nonzero()[0]

        newB = []

        # linearily interpolate the Nans
        for i in range(3):
            Bcomp =DeltaB_ENU[:, i]
            nans, x = nan_helper(Bcomp)
            Bcomp[nans] = np.interp(x(nans),x(~nans),Bcomp[~nans])
            newB.append(Bcomp)


        DeltaB_ENU = np.array([ [newB[0][i],newB[1][i],newB[2][i]] for i in range(len(DeltaB_ENU))])

    ###########################################
    # --- Add DeltaB to Modified CHAOS data ---
    ###########################################
    if include_CHAOS_model:
        data_for_output = B_CHAOS_ENU_magTime + DeltaB_ENU
    else:
        data_for_output = DeltaB_ENU

    if outputData:

        # create the output data_dict

        comps = ['Bx', 'By', 'Bz']
        newComps = ['B_East', 'B_North', 'B_Up', 'Bmag']
        data_for_output_despin = np.array([[data_for_output[i][0], data_for_output[i][1], data_for_output[i][2]] for i in range(len(data_for_output))])

        # --- Magnetic Components ---
        # get the attributes of the old components and replace them
        for i, key in enumerate(comps):
            newAttrs = deepcopy(data_dict_output[key][1])
            newAttrs['LABLAXIS'] = newComps[i]

            # remove the old key
            del data_dict_output[key]

            # append the new key
            data_dict_output = {**data_dict_output, **{newComps[i]: [data_for_output_despin[:, i], newAttrs]}}

        #####################################################
        # --- APPLY T0 KENTON/ANTONIO CORRECTION TO B_ENU ---
        #####################################################
        # [1] Create a new Epoch variable through multiplication
        # [2] interpolate the B-Field data onto the old time-base?
        # [3] export the newly interpolated B-Field data
        if DespinToggles.KentonAntonio_T0_Correction:
            if wRocket == 5:
                stl.prgMsg('Applying Kenton/Antonio T0 correction')
                # Calculate old timebase
                T0_time = dt.datetime(2022, 11, 20, 17, 20, 00, 000000)
                oldEpoch_TSL = stl.EpochTo_T0_Rocket(InputEpoch=data_dict_output['Epoch'][0], T0=T0_time)
                slope = DespinToggles.EB_East[0]
                interscept = DespinToggles.EB_East[1]
                newEpoch = np.array([ pycdf.lib.tt2000_to_datetime(int((slope * val + interscept)*1E9 + pycdf.lib.datetime_to_tt2000(T0_time))) for val in oldEpoch_TSL])

                # interpolate onto old timebase
                data_dict_output_interp = stl.InterpolateDataDict(InputDataDict=data_dict_output,
                                                              InputEpochArray=newEpoch,
                                                              targetEpochArray=deepcopy(data_dict_output['Epoch'][0]),
                                                              wKeys=[])
                data_dict_output = deepcopy(data_dict_output_interp)

                stl.Done(start_time)

        stl.prgMsg('Creating output file')
        if DespinToggles.KentonAntonio_T0_Correction:
            data_dict_output = {**data_dict_output, **{'Kenton_time_correction_was_applied':[np.array([1]), {'VAR_TYPE':'data'}]}}

        if include_CHAOS_model:
            data_dict_output = {**data_dict_output,
                                **{'CHAOS_BField_Model_Included': [np.array([1]), {'VAR_TYPE': 'data'}]}}
        else:
            data_dict_output = {**data_dict_output,
                                **{'CHAOS_BField_Model_NOT_Included': [np.array([1]), {'VAR_TYPE': 'data'}]}}

        outputPath = f'{rocketFolderPath}{outputPath_modifier_despin}\{ACESII.fliers[wflyer]}\\{fileoutName_despin}.cdf'
        stl.outputCDFdata(outputPath, data_dict_output, instrNam='RingCore')
        stl.Done(start_time)











# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 4:  # ACES II High
    rocketFolderPath = DataPaths.ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Low
    rocketFolderPath = DataPaths.ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wflyer]}\*.cdf')) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    if justPrintFileNames:
        RingCore_L1_to_L2_Despin(wRocket, 0, rocketFolderPath, justPrintFileNames,wflyer)
    elif not wFiles:
        for fileNo in (range(len(glob(f'{rocketFolderPath}{inputPath_modifier}\{ACESII.fliers[wflyer]}\*.cdf')))):
            RingCore_L1_to_L2_Despin(wRocket, fileNo, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        for filesNo in wFiles:
            RingCore_L1_to_L2_Despin(wRocket, filesNo, rocketFolderPath, justPrintFileNames,wflyer)