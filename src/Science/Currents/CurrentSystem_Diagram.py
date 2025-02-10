# # --- CurrentSystem_Diagram.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: Create a plot the describes (roughly) how the current system in ACESII should look



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *

start_time = time.time()
# --- --- --- --- ---



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
targetTimes = [dt.datetime(2022, 11, 20, 17, 24, 30, 000000), dt.datetime(2022, 11, 20, 17, 27, 30, 000000)]

subtractCHAOS = True

# filter toggles
filterData = True
lowCut_toggle, highcut_toggle, filttype_toggle, order_toggle = 0.4, 0.02, 'LowPass', 4  # filter toggles LOW FLYER

# ESA toggles
wPitch = 2
Vmin, Vmax = 5E5, 3E8


# E-Field Adjust
AdjustEField = True
E_Field_adjust = [120, -15, 0]

# plot the detrended data
plotfilteredData = False

# remove edge - The filters introduce bad data. Reduce this % of data from either side
edgePercent = 2

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from matplotlib.gridspec import GridSpec
from ACESII_code.class_var_func import butter_filter
from scipy.interpolate import CubicSpline

def reduceData(targetTimes,data_dict):
    lowCutoff, highCutoff = np.abs(data_dict['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict['Epoch'][0] - targetTimes[1]).argmin()
    for key, val in data_dict.items():
        data_dict[key][0] = np.array(data_dict[key][0][lowCutoff:highCutoff])
    return data_dict
# none


def CurrentSystem_Diagram(rocketFolderPath):

    prgMsg('Loading and Reducing Attitude Data')

    # get the attitude input files HIGH
    inputFile_attitude_high = glob(r'C:\Data\ACESII\attitude\high\*.cdf')[0]
    data_dict_attitude_high = loadDictFromFile(inputFile_attitude_high, {})
    data_dict_attitude_high = reduceData(targetTimes,data_dict_attitude_high)
    data_dict_attitude_high['Alt'][0] = data_dict_attitude_high['Alt'][0]/1000

    # get the attitude input files LOW
    inputFile_attitude_low = glob(r'C:\Data\ACESII\attitude\low\*.cdf')[0]
    data_dict_attitude_low = loadDictFromFile(inputFile_attitude_low, {})
    data_dict_attitude_low = reduceData(targetTimes, data_dict_attitude_low)
    data_dict_attitude_low['Alt'][0] = data_dict_attitude_low['Alt'][0] / 1000

    Done(start_time)
    prgMsg('Loading and Reducing Magnetometer Data')

    # get the MAG input files LOW
    inputFile_mag_low = glob(r'C:\Data\ACESII\L2\low\*RingCore_DeSpun*')[0]
    data_dict_mag_low = loadDictFromFile(inputFile_mag_low, {})
    data_dict_mag_low = reduceData(targetTimes, data_dict_mag_low)
    B_Field_low = np.array([[data_dict_mag_low['B_east'][0][i], data_dict_mag_low['B_north'][0][i], data_dict_mag_low['B_up'][0][i]] for i in range(len(data_dict_mag_low['Epoch'][0]))])

    # get the MAG input files HIGH
    inputFile_mag_high = glob(r'C:\Data\ACESII\L2\high\*RingCore_DeSpun*')[0]
    data_dict_mag_high = loadDictFromFile(inputFile_mag_high, {})
    data_dict_mag_high = reduceData(targetTimes, data_dict_mag_high)
    B_Field_high = np.array([[data_dict_mag_high['B_east'][0][i], data_dict_mag_high['B_north'][0][i], data_dict_mag_high['B_up'][0][i]] for i in range(len(data_dict_mag_high['Epoch'][0]))])

    Done(start_time)
    prgMsg('Loading and Reducing EEPAA Data')

    inputFile_EEPAA_high = glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')[0]
    data_dict_eepaa_high = loadDictFromFile(inputFile_EEPAA_high,{})
    lowCutoff, highCutoff = np.abs(data_dict_eepaa_high['Epoch_esa'][0] - targetTimes[0]).argmin(), np.abs(data_dict_eepaa_high['Epoch_esa'][0] - targetTimes[1]).argmin()
    for key, val in data_dict_eepaa_high.items():
        if key.lower() in ['eepaa','epoch_esa','differential_energy_flux']:
            data_dict_eepaa_high[key][0] = np.array(data_dict_eepaa_high[key][0][lowCutoff:highCutoff])


    inputFile_EEPAA_low = glob('C:\Data\ACESII\L2\low\*eepaa_fullCal*')[0]
    data_dict_eepaa_low = loadDictFromFile(inputFile_EEPAA_low, {})
    lowCutoff, highCutoff = np.abs(data_dict_eepaa_low['Epoch_esa'][0] - targetTimes[0]).argmin(), np.abs(data_dict_eepaa_low['Epoch_esa'][0] - targetTimes[1]).argmin()
    for key, val in data_dict_eepaa_low.items():
        if key.lower() in ['eepaa', 'epoch_esa', 'differential_energy_flux']:
            data_dict_eepaa_low[key][0] = np.array(data_dict_eepaa_low[key][0][lowCutoff:highCutoff])

    Done(start_time)
    prgMsg('Loading and Reducing Electric Field Data')

    # get the EFI input files LOW
    inputFile_elec_low = glob(r'C:\Data\ACESII\L2\low\*E_Field*')[0]
    data_dict_elec_low = loadDictFromFile(inputFile_elec_low, {})
    data_dict_elec_low = reduceData(targetTimes, data_dict_elec_low)
    E_Field = 1000*np.array([ [data_dict_elec_low['E_East'][0][i],data_dict_elec_low['E_North'][0][i],data_dict_elec_low['E_Up'][0][i]] for i in range(len(data_dict_elec_low['Epoch'][0]))])

    if E_Field_adjust:
        E_Field = E_Field + E_Field_adjust

    Done(start_time)

    # --- --- --- --- --- ---
    # --- SUBTRACT CHAOS  ---
    # --- --- --- --- --- ---
    if subtractCHAOS:

        # interpolate Attitude data up to magnetometer data

        # LF
        prgMsg('Interpolating LF Attitude Data')
        offsetResults_intercept = [0, 0]
        Epoch_attitude_TT2000 = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_attitude_low['Epoch'][0]])
        Epoch_mag_TT2000 = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_mag_low['Epoch'][0]])
        Epoch_attitude_loop = np.array([int(tme + offsetResults_intercept[5 - 4]) for tme in Epoch_attitude_TT2000])
        dataInterp_dict_attitude_low = {'Epoch':deepcopy(data_dict_mag_low['Epoch'][0]),
                                    'Alt':deepcopy(data_dict_attitude_low['Alt'][0]),
                                    'Latgd':deepcopy(data_dict_attitude_low['Latgd'][0]),
                                    'Long':deepcopy(data_dict_attitude_low['Long'][0])}

        for key, newDataList in dataInterp_dict_attitude_low.items():
            if key != 'Epoch':
                # --- cubic interpolation ---
                splCub = CubicSpline(Epoch_attitude_loop, dataInterp_dict_attitude_low[key])

                # evaluate the interpolation at all the epoch_mag points
                dataInterp_dict_attitude_low[key] = np.array([splCub(timeVal) for timeVal in Epoch_mag_TT2000])

        Done(start_time)

        # HF
        prgMsg('Interpolating HF Attitude Data')
        offsetResults_intercept = [0, 0]
        Epoch_attitude_TT2000 = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_attitude_high['Epoch'][0]])
        Epoch_mag_TT2000 = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_mag_high['Epoch'][0]])
        Epoch_attitude_loop = np.array([int(tme + offsetResults_intercept[4 - 4]) for tme in Epoch_attitude_TT2000])
        dataInterp_dict_attitude_high = {'Epoch': deepcopy(data_dict_mag_high['Epoch'][0]),
                                    'Alt': deepcopy(data_dict_attitude_high['Alt'][0]),
                                    'Latgd': deepcopy(data_dict_attitude_high['Latgd'][0]),
                                    'Long': deepcopy(data_dict_attitude_high['Long'][0])}

        for key, newDataList in dataInterp_dict_attitude_high.items():
            if key != 'Epoch':
                # --- cubic interpolation ---
                splCub = CubicSpline(Epoch_attitude_loop, dataInterp_dict_attitude_high[key])

                # evaluate the interpolation at all the epoch_mag points
                dataInterp_dict_attitude_high[key] = np.array([splCub(timeVal) for timeVal in Epoch_mag_TT2000])

        Done(start_time)


        prgMsg('Calculating CHAOS magnetic Field and Detrending')
        from ACESII_code.class_var_func import CHAOS
        B_CHAOS_low = np.array( CHAOS(dataInterp_dict_attitude_low['Latgd'],
                                      dataInterp_dict_attitude_low['Long'],
                                      dataInterp_dict_attitude_low['Alt'],
                                      dataInterp_dict_attitude_low['Epoch']))

        B_CHAOS_high = np.array(CHAOS(dataInterp_dict_attitude_high['Latgd'],
                                     dataInterp_dict_attitude_high['Long'],
                                     dataInterp_dict_attitude_high['Alt'],
                                     dataInterp_dict_attitude_high['Epoch']))

        # Subtracting the Chaos Model
        B_Field_low = B_Field_low - B_CHAOS_low
        B_Field_high = B_Field_high - B_CHAOS_high

        # --- detrend the data ---
        B_Field_detrend_low = np.array([scipy.signal.detrend(B_Field_low[:, 0]),
                                scipy.signal.detrend(B_Field_low[:, 1]),
                                scipy.signal.detrend(B_Field_low[:, 2])])

        B_Field_detrend_high = np.array([scipy.signal.detrend(B_Field_high[:, 0]),
                                        scipy.signal.detrend(B_Field_high[:, 1]),
                                        scipy.signal.detrend(B_Field_high[:, 2])])

        # Store the outputs
        B_Field_low = np.array([[B_Field_detrend_low[0][i], B_Field_detrend_low[1][i], B_Field_detrend_low[2][i]] for i in range(len(B_Field_low))])
        B_Field_high = np.array([[B_Field_detrend_high[0][i], B_Field_detrend_high[1][i], B_Field_detrend_high[2][i]] for i in range(len(B_Field_high))])

        Done(start_time)


    # --- --- --- --- --- ---
    # --- FILTER THE DATA ---
    # --- --- --- --- --- ---
    if filterData:
        prgMsg('Filtering Data')

        # filter the ULTRA DC component from the magnetometer data

        # B - LOW
        B_filtered_low = []
        for i in range(3):
            B_filtered_low.append(butter_filter(B_Field_low[:, i], lowcutoff=lowCut_toggle, highcutoff=highcut_toggle, filtertype=filttype_toggle, order=order_toggle, fs=128))
        B_Field_low = np.array([[B_filtered_low[0][i],B_filtered_low[1][i],B_filtered_low[2][i]] for i in range(len(B_Field_low))])


        # B - HIGH
        B_filtered_high = []
        for i in range(3):
            B_filtered_high.append(butter_filter(B_Field_high[:, i], lowcutoff=lowCut_toggle, highcutoff=highcut_toggle, filtertype=filttype_toggle, order=order_toggle, fs=128))
        B_Field_high = np.array([[B_filtered_high[0][i], B_filtered_high[1][i], B_filtered_high[2][i]] for i in range(len(B_Field_high))])

        # E - LOW
        # E_filtered = []
        # for i in range(3):
        #     E_filtered.append(butter_filter(E_Field[:, i], lowcutoff=lowCut_toggle, highcutoff=highcut_toggle, filtertype=filttype_toggle, order=order_toggle, fs=4000))
        # E_Field = np.array([[E_filtered[0][i], E_filtered[1][i], E_filtered[2][i]] for i in range(len(E_Field))])

        Done(start_time)

    # --- --- --- --- --- ---
    # --- REDUCE THE DATA ---
    # --- --- --- --- --- ---
    prgMsg('Reducing Edge Percent of Data')
    percent = edgePercent/100

    # E-Field
    dictonary = data_dict_elec_low
    N = len(dictonary['Epoch'][0])
    lowEdge,highEdge = int(percent*N),int((1- percent)*N)
    for key, val in dictonary.items():
        dictonary[key][0] = dictonary[key][0][lowEdge:highEdge]
    E_Field = E_Field[lowEdge:highEdge]

    # B-Field High
    dictonary = data_dict_mag_high
    N = len(dictonary['Epoch'][0])
    lowEdge, highEdge = int(percent * N), int((1 - percent) * N)
    for key, val in dictonary.items():
        dictonary[key][0] = dictonary[key][0][lowEdge:highEdge]
    B_Field_high = B_Field_high[lowEdge:highEdge]

    # B-Field Low
    dictonary = data_dict_mag_low
    N = len(dictonary['Epoch'][0])
    lowEdge, highEdge = int(percent * N), int((1 - percent) * N)
    for key, val in dictonary.items():
        dictonary[key][0] = dictonary[key][0][lowEdge:highEdge]
    B_Field_low = B_Field_low[lowEdge:highEdge]

    # ESA High
    N = len(data_dict_eepaa_high['Epoch_esa'][0])
    lowEdge, highEdge = int(percent * N), int((1 - percent) * N)
    data_dict_eepaa_high['Differential_Energy_Flux'][0] = data_dict_eepaa_high['Differential_Energy_Flux'][0][lowEdge:highEdge]
    data_dict_eepaa_high['Epoch_esa'][0] = data_dict_eepaa_high['Epoch_esa'][0][lowEdge:highEdge]

    # ESA Low
    N = len(data_dict_eepaa_low['Epoch_esa'][0])
    lowEdge, highEdge = int(percent * N), int((1 - percent) * N)
    data_dict_eepaa_low['Differential_Energy_Flux'][0] = data_dict_eepaa_low['Differential_Energy_Flux'][0][lowEdge:highEdge]
    data_dict_eepaa_low['Epoch_esa'][0] = data_dict_eepaa_low['Epoch_esa'][0][lowEdge:highEdge]
    Done(start_time)


    # Plot the Filtered High Flyer Magnetometer Data
    if plotfilteredData:

        fig, ax = plt.subplots(3)
        ylabels = ['$\Delta$B_East', '$\Delta$B_North', '$\Delta$B_Up']
        for i in range(3):
            ax[i].plot(data_dict_mag_high['Epoch'][0], B_Field_high[:, i])
            ax[i].set_ylabel(ylabels[i])

        ax[2].set_xlabel('Epoch')
        plt.show()


    # --- --- --- --- --- -
    # --- MAKE THE PLOT ---
    # --- --- --- --- --- -

    prgMsg('Making the Plot')

    fig = plt.figure()
    plt.subplots_adjust(wspace=0.1, hspace=0)
    gs = GridSpec(6, 2,width_ratios=[1,0.01])


    # E_Field Components
    axE = fig.add_subplot(gs[4:6, 0])
    Epoch = data_dict_elec_low['Epoch'][0]
    compsEfield = ['E_East (adjusted)','E_North (adjusted)', 'E_Up (adjusted)'] if E_Field_adjust else ['E_East','E_North', 'E_Up']
    axE.plot(Epoch, E_Field[:, 0], label=compsEfield[0], color='tab:blue')
    axE.plot(Epoch, E_Field[:, 1], label=compsEfield[1], color='tab:red')
    axE.plot(Epoch, E_Field[:, 2], label=compsEfield[2], color='tab:orange')
    axE.grid(True, axis='both')
    axE.legend(loc='best')
    axE.set_ylabel('E-Field [mV/m]')
    axE.set_ylim(-40,110)

    # B-Plot LOW
    axB_low = fig.add_subplot(gs[2:4, 0], sharex=axE)
    axB_low.tick_params(labelbottom=False)
    Epoch = data_dict_mag_low['Epoch'][0]
    axB_low.plot(Epoch, B_Field_low[:, 0], label='B_East', color='tab:blue')
    axB_low.plot(Epoch, B_Field_low[:, 1], label='B_North', color='tab:red')
    axB_low.plot(Epoch, B_Field_low[:, 2], label='B_Up', color='tab:orange')
    axB_low.grid(True, axis='both')
    axB_low.set_ylabel('$\Delta$ B [nT]')
    axB_low.legend(loc='best')
    axB_low.set_ylim(-57,40)

    # High Flyer - EEPAA
    axESA_high = fig.add_subplot(gs[0, 0],sharex=axE)
    axESA_high.tick_params(labelbottom=False)
    ESAdata = data_dict_eepaa_high['Differential_Energy_Flux'][0][:,wPitch,:]
    Epoch = data_dict_eepaa_high['Epoch_esa'][0]
    X,Y = np.meshgrid(Epoch, data_dict_eepaa_high['Energy'][0])
    cmap = axESA_high.pcolormesh(X, Y, np.transpose(ESAdata), vmin=Vmin, vmax=Vmax,cmap='turbo')
    pitches = data_dict_eepaa_high['Pitch_Angle'][0]
    axESA_high.set_title(r'$\alpha$='+ f'{pitches[wPitch]}' +'$^{\circ}$')
    axESA_high.set_ylabel('Energy [eV]')
    axESA_high.set_ylim(20, 1000)

    # High Flyer - ALTITUDE on EEPAA PLOT
    axESA_high_ALT = axESA_high.twinx()
    interpEpoch = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_eepaa_high['Epoch_esa'][0]])
    Epoch_attitude_high = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in dataInterp_dict_attitude_high['Epoch']])
    splCub = CubicSpline(Epoch_attitude_high, dataInterp_dict_attitude_high['Alt'])  # interpolate altitude onto EEPAA Epoch
    Alt_high = np.array([splCub(timeVal) for timeVal in interpEpoch])
    axESA_high_ALT.plot(Epoch, Alt_high, color='white', label='High Flyer Trajectory')
    axESA_high_ALT.set_ylabel('Alt [km]')
    axESA_high_ALT.legend(loc='best')

    # Low Flyer - EEPAA
    axESA_low = fig.add_subplot(gs[1, 0], sharex=axE)
    axESA_low.tick_params(labelbottom=False)
    ESAdata = data_dict_eepaa_low['Differential_Energy_Flux'][0][:,wPitch, :]
    Epoch = data_dict_eepaa_low['Epoch_esa'][0]
    X, Y = np.meshgrid(Epoch, data_dict_eepaa_low['Energy'][0])
    cmap = axESA_low.pcolormesh(X, Y, np.transpose(ESAdata), vmin=Vmin, vmax=Vmax,cmap='turbo')
    axESA_low.set_ylabel('Energy [eV]')
    axESA_low.set_ylim(20, 1000)

    # Low Flyer - ALTITUDE on EEPAA PLOT
    axESA_low_ALT = axESA_low.twinx()
    interpEpoch = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_eepaa_low['Epoch_esa'][0]])
    Epoch_attitude_low = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in dataInterp_dict_attitude_low['Epoch']])
    splCub = CubicSpline(Epoch_attitude_low, dataInterp_dict_attitude_low['Alt']) # interpolate altitude onto EEPAA Epoch
    Alt_low = np.array([splCub(timeVal) for timeVal in interpEpoch])
    axESA_low_ALT.plot(Epoch,Alt_low,color='white',label='Low Flyer Trajectory')
    axESA_low_ALT.set_ylabel('Alt [km]')
    axESA_low_ALT.legend(loc='best')

    # colorbar for ESA data
    axColorbar = fig.add_subplot(gs[0:2,1])
    cbar = plt.colorbar(cmap, cax=plt.subplot(gs[0:2, 1]))
    cbar.set_label('Differential Energy Flux')



    # B-Plot HIGH
    # axB_high = fig.add_subplot(gs[1, :], sharex=axE)
    # Epoch = data_dict_mag_high['Epoch'][0]
    # axB_high.plot(Epoch, B_Field_high[:, 0], label='B_East', color='tab:blue')
    # axB_high.plot(Epoch, B_Field_high[:, 1], label='B_North', color='tab:red')
    # # axB_high.plot(Epoch, B_Field_high[:, 2], label='B_Up', color='tab:orange')
    # axB_high.grid(True, axis='both')
    # axB_high.set_ylabel('B_High [nT]')



    plt.show()

    Done(start_time)














# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
CurrentSystem_Diagram(rocketFolderPath)