# --- LP_postFlight_cal_collect_data.py ---
# Desciption: Determine the data which will be used to calibrate
# the Fixed Langmuir Probes:
# (1) EISCAT density profiles
# (2) DERPA Te data + EISCAT Tr = Ti/Te data
# (3) Floating Potential Data
# (4) Rocket Altitude
# Interpolate Everything on the Langmuir Probe timebase

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"

from src.ACESII.my_imports import *
import time
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False  # Just print the names of files

# --- Select the Rocket ---
# 0 -> Integration High Flier
# 1 -> Integration Low Flier
# 2 -> TRICE II High Flier
# 3 -> TRICE II Low Flier
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5
showEISCAT_profiles = False
altitude_cutoff_upleg = [100, 166]  # in [km]. Everything below this altitude uses the ul EISCAT profiles, above it uses md and below it on downleg uses dl
altitude_cutoff_downleg = [200, 183]  # in [km]. Everything below this altitude uses the ul EISCAT profiles, above it uses md and below it on downleg uses dl
titles = ['upleg', 'middle', 'downleg']

# --- OutputData ---
outputData = True

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from glob import glob
from scipy.signal import savgol_filter
import spaceToolsLib as stl
import datetime as dt


#######################
# --- MAIN FUNCTION ---
#######################
def LP_collect_postFlight_cal_data(wRocket):

    # --- FILE I/O ---
    stl.prgMsg('Loading Data')

    # load the rocket ion saturation current data
    data_path = glob(rf'C:\Data\ACESII\L2\{ACESII.fliers[wRocket-4]}\\*langmuir_fixed.cdf*')[0]
    data_dict_LP_current = stl.loadDictFromFile(data_path)

    # load the EISCAT data
    data_path = r'C:\Data\ACESII\science\EISCAT\tromso\UHF\MAD6400_2022-11-20_beata_ant@uhfa.cdf'
    data_dict_EISCAT_tromso = stl.loadDictFromFile(data_path)

    data_path = r'C:\Data\ACESII\science\EISCAT\svalbard\UKFI_radar\MAD6400_2022-11-20_ipy_60@42m.cdf'
    data_dict_EISCAT_svalbard = stl.loadDictFromFile(data_path)

    # load the attitude data
    data_path = glob(rf'C:\Data\ACESII\attitude\{ACESII.fliers[wRocket-4]}\\*.cdf*')[0]
    data_dict_attitude = stl.loadDictFromFile(data_path)

    # load the trajectory data
    data_path = glob(rf'C:\Data\ACESII\trajectories\{ACESII.fliers[wRocket - 4]}\\*_auroral.cdf*')[0]
    data_dict_traj = stl.loadDictFromFile(data_path)

    # load the L-Shell data
    data_path = glob(rf'C:\Data\ACESII\coordinates\Lshell\{ACESII.fliers[wRocket - 4]}\\*.cdf*')[0]
    data_dict_LShell = stl.loadDictFromFile(data_path)

    # load the simulation m_eff data
    data_path = glob(rf'C:\Data\physicsModels\ionosphere\plasma_environment\\plasma_environment.cdf')[0]
    data_dict_sim = stl.loadDictFromFile(data_path)

    # load the DERPA Te data
    data_path = glob(rf'C:\Data\ACESII\L2\{ACESII.fliers[wRocket - 4]}\\*_ERPA1.cdf*')[0]
    data_dict_DERPA1 = stl.loadDictFromFile(data_path)
    data_path = glob(rf'C:\Data\ACESII\L2\{ACESII.fliers[wRocket - 4]}\\*_ERPA2.cdf*')[0]
    data_dict_DERPA2 = stl.loadDictFromFile(data_path)

    # load the floating potential data
    data_path = glob(rf'C:\Data\ACESII\science\payload_potential\{ACESII.fliers[wRocket - 4]}\\*.cdf*')[0]
    data_dict_payload_potential = stl.loadDictFromFile(data_path)

    # --- prepare the output ---
    data_dict_output = {}
    stl.Done(start_time)

    # --- shorten the dataset ---
    # if wRocket-4 == 0:
    #     low_idx = np.abs(data_dict_LP_current['Epoch'][0] - dt.datetime(2022, 11, 20, 17, 21, 00)).argmin()
    #     high_idx = np.abs(data_dict_LP_current['Epoch'][0] - dt.datetime(2022, 11, 20, 17, 29,50)).argmin()
    # elif wRocket-4 ==1 :
    #     low_idx = np.abs(data_dict_LP_current['Epoch'][0] - dt.datetime(2022, 11, 20, 17, 23, 20)).argmin()
    #     high_idx = np.abs(data_dict_LP_current['Epoch'][0] - dt.datetime(2022, 11, 20, 17, 28)).argmin()
    # data_dict_LP_current['fixed_current'][0]= data_dict_LP_current['fixed_current'][0][low_idx:high_idx+1]
    # data_dict_LP_current['Epoch'][0] = data_dict_LP_current['Epoch'][0][low_idx:high_idx + 1]

    ###############################################
    # --- Collect the EISCAT Data and Smooth it ---
    ###############################################

    # --- UPLEG BACKGROUND PROFILES (Tromso) ---
    target_times_ul = [
                    dt.datetime(2022,11,20,17,21,25),
                    dt.datetime(2022,11,20,17,22,10),
                    # dt.datetime(2022,11,20,17,25,25)
                    ]


    # --- MIDDLE BACKGROUND PROFILES (Tromso) ---
    target_times_md = [
        dt.datetime(2022, 11, 20, 17, 22, 10),
        dt.datetime(2022, 11, 20, 17, 23, 25),
        dt.datetime(2022, 11, 20, 17, 26, 10),
    ]

    # --- DOWNLEG BACKGROUND PROFILES (Svalbard) ---
    # target_times_dl = [
    #                     # dt.datetime(2022, 11, 20, 17, 24, 55),
    #                    dt.datetime(2022, 11, 20, 17, 27, 25),
    #                    dt.datetime(2022, 11, 20, 17, 28, 10),
    #                    ]

    target_times_dl = [
        dt.datetime(2022, 11, 20, 17, 21, 25),
        dt.datetime(2022, 11, 20, 17, 22, 10),
    ]

    target_times_collection = [target_times_ul, target_times_md, target_times_dl]
    # data_dict_EISCAT = [deepcopy(data_dict_EISCAT_tromso),deepcopy(data_dict_EISCAT_tromso),deepcopy(data_dict_EISCAT_svalbard)]
    data_dict_EISCAT = [deepcopy(data_dict_EISCAT_tromso), deepcopy(data_dict_EISCAT_tromso), deepcopy(data_dict_EISCAT_tromso)]


    Ti_profiles_val = []
    Ti_profiles_alt = []
    Tr_profiles_val = []
    Tr_profiles_alt = []
    ne_profiles_val = []
    ne_profiles_alt = []

    for idx, target_times in enumerate(target_times_collection):
        target_idxs = [np.abs(tme - data_dict_EISCAT[idx]['Epoch'][0]).argmin() for tme in target_times]
        target_Ti = data_dict_EISCAT[idx]['ti'][0][target_idxs]
        target_ni = data_dict_EISCAT[idx]['ne'][0][target_idxs]
        target_Tr = data_dict_EISCAT[idx]['tr'][0][target_idxs]

        # Average the profiles
        # Ti
        Ti_profile = np.nanmean(target_Ti,axis=0) * (stl.kB/stl.q0)
        Ti_profile_idxs = np.isnan(Ti_profile)
        Ti_profile_alts = data_dict_EISCAT[idx]['range'][0][np.where(Ti_profile_idxs==False)[0]]
        Ti_profile_vals = Ti_profile[np.where(Ti_profile_idxs==False)[0]]

        # Tr
        Tr_profile = np.nanmean(target_Tr, axis=0)
        Tr_profile_idxs = np.isnan(Tr_profile)
        Tr_profile_alts = data_dict_EISCAT[idx]['range'][0][np.where(Tr_profile_idxs == False)[0]]
        Tr_profile_vals = Tr_profile[np.where(Tr_profile_idxs == False)[0]]

        # ni
        ne_profile = np.nanmean(target_ni, axis=0)
        ne_profile_idxs = np.isnan(ne_profile)
        ne_profile_alts = data_dict_EISCAT[idx]['range'][0][np.where(ne_profile_idxs==False)[0]]
        ne_profile_vals = ne_profile[np.where(ne_profile_idxs==False)[0]]

        # ---Smooth the profiles ---
        window_length = 20
        porder = 3
        Ti_profile_smooth = savgol_filter(x=Ti_profile_vals,window_length=window_length,polyorder=porder)
        Ti_profile_smoothed = savgol_filter(x=Ti_profile_smooth, window_length=window_length+20, polyorder=porder+2)

        Tr_profile_smooth = savgol_filter(x=Tr_profile_vals,window_length=window_length,polyorder=porder)
        Tr_profile_smoothed = savgol_filter(x=Tr_profile_smooth, window_length=window_length+20, polyorder=porder+2)

        ne_profile_smooth = savgol_filter(x=ne_profile_vals, window_length=window_length, polyorder=porder)
        ne_profile_smoothed = savgol_filter(x=ne_profile_smooth, window_length=window_length+20, polyorder=porder + 2)

        # store everything
        Ti_profiles_val.append(Ti_profile_smoothed)
        Ti_profiles_alt.append(Ti_profile_alts)

        Tr_profiles_val.append(Tr_profile_smoothed)
        Tr_profiles_alt.append(Tr_profile_alts)

        ne_profiles_val.append(ne_profile_smoothed)
        ne_profiles_alt.append(ne_profile_alts)


        if showEISCAT_profiles:
            import matplotlib.pyplot as plt
            fig, ax =plt.subplots(nrows=3,ncols=1)

            fig.suptitle(titles[idx])

            ax[0].plot(Ti_profile_alts, Ti_profile_vals, color='tab:blue')
            ax[0].plot(Ti_profile_alts, Ti_profile_smoothed, color='tab:red')
            ax[0].set_xlabel('Alt [km]')
            ax[0].set_ylabel('Ti [eV]')
            ax[0].set_xlim(0,400)

            ax[1].plot(ne_profile_alts, ne_profile_vals,color='tab:blue')
            ax[1].plot(ne_profile_alts, ne_profile_smoothed, color='tab:red')
            ax[1].set_yscale('log')
            ax[1].set_xlabel('Alt [km]')
            ax[1].set_ylabel('n [m^-3]')
            ax[1].set_xlim(0, 400)

            ax[2].plot(Tr_profile_alts, Tr_profile_vals, color='tab:blue')
            ax[2].plot(Tr_profile_alts, Tr_profile_smoothed, color='tab:red')
            ax[2].set_xlabel('Alt [km]')
            ax[2].set_ylabel('Tr')
            ax[2].set_xlim(0, 400)
            plt.show()

    ###############################################################################
    # --- Interpolate the Ti, Te, ni, mi data onto the langmuir current timebase ---
    ###############################################################################

    # --- Break the LP data into three sections ---

    # Interpolate the Altitude
    T0 = dt.datetime(2022, 11, 20, 17, 19)
    T0_LP_current = np.array(stl.EpochTo_T0_Rocket(data_dict_LP_current['Epoch'][0], T0=T0),dtype='float64')
    T0_attitude = stl.EpochTo_T0_Rocket(data_dict_attitude['Epoch'][0], T0=T0)
    alt_langmuir = np.interp(np.array(T0_LP_current,dtype='float64'), np.array(T0_attitude,dtype='float64'), np.array(data_dict_attitude['Alt'][0]/stl.m_to_km,dtype='float64'))

    # Interpolate the L-Shell
    T0_LShell  = stl.EpochTo_T0_Rocket(data_dict_LShell['Epoch'][0],T0=T0)
    LShell_langmuir = np.interp(T0_LP_current, T0_LShell,data_dict_LShell['L-Shell'][0])

    # Interpolate the LP_swept Floating Potential
    good_idxs = np.where(np.isnan(data_dict_payload_potential['floating_potential'][0])==False)
    T0_floating = stl.EpochTo_T0_Rocket(data_dict_payload_potential['Epoch'][0][good_idxs], T0=T0)
    floatingPotential_langmuir = np.interp(T0_LP_current, T0_floating, data_dict_payload_potential['floating_potential'][0][good_idxs])


    # Interpolate the DERPA Te Data

    # 1
    T0_DERPA1 = stl.EpochTo_T0_Rocket(data_dict_DERPA1['Epoch'][0],T0=T0)
    if wRocket-4 == 0:
        low_idx = np.abs(data_dict_DERPA1['Epoch'][0] - dt.datetime(2022, 11, 20, 17, 21, 37)).argmin()
        low_val = 0.09
    elif wRocket-4 == 1:
        low_idx = np.abs(data_dict_DERPA1['Epoch'][0] - dt.datetime(2022, 11, 20, 17, 23, 50)).argmin()
        low_val = 0.06

    Temp_DERPA1 = data_dict_DERPA1['temperature'][0]
    Temp_DERPA1[:low_idx] = low_val
    Te_DERPA1_langmuir = np.interp(T0_LP_current, T0_DERPA1, Temp_DERPA1)

    T0_DERPA2 = stl.EpochTo_T0_Rocket(data_dict_DERPA2['Epoch'][0], T0=T0)
    Temp_DERPA2 = data_dict_DERPA2['temperature'][0]
    Temp_DERPA2[:low_idx] = low_val
    Te_DERPA2_langmuir = np.interp(T0_LP_current, T0_DERPA2, Temp_DERPA2)


    # first find the peak altitude, beak datasets into upleg, middle and downleg
    max_alt_idx = np.abs(alt_langmuir- np.max(alt_langmuir)).argmin()
    upleg_idx = np.abs(alt_langmuir[0:max_alt_idx]-altitude_cutoff_upleg[wRocket-4]).argmin()
    downleg_idx = np.abs(alt_langmuir[max_alt_idx::]-altitude_cutoff_downleg[wRocket-4]).argmin()
    separation_idxs = [upleg_idx, len(alt_langmuir[0:max_alt_idx])+ downleg_idx] # the two points where the data is separated into three sections

    # Interpolate data onto the three EISCAT profiles
    Ti_interp = []
    Tr_interp = []
    ne_interp = []
    m_eff_i_interp = []

    for idx in range(3):
        if idx == 0:
            low_idx = 0
            high_idx = separation_idxs[0]
        elif idx == 1:
            low_idx = separation_idxs[0]
            high_idx = separation_idxs[1]
        elif idx == 2:
            low_idx = separation_idxs[1]
            high_idx = len(alt_langmuir)

        # get the specific EISCAT profiles

        # Create the EISCAT Ti profile for Langmuir
        Ti_langmuir = np.interp(alt_langmuir[low_idx:high_idx], Ti_profiles_alt[idx], Ti_profiles_val[idx])

        # Create the ni EISCAT profile for langmuir
        Tr_langmuir = np.interp(alt_langmuir[low_idx:high_idx], Tr_profiles_alt[idx], Tr_profiles_val[idx])

        # Create the ni EISCAT profile for langmuir
        ne_langmuir = np.interp(alt_langmuir[low_idx:high_idx], ne_profiles_alt[idx], ne_profiles_val[idx])

        # Create the m_eff_i profile for langmuir
        m_eff_i = np.interp(alt_langmuir[low_idx:high_idx],data_dict_sim['simAlt'][0]/stl.m_to_km,data_dict_sim['m_eff_i'][0][0,:])

        # store everything
        Ti_interp.append(Ti_langmuir)
        Tr_interp.append(Tr_langmuir)
        ne_interp.append(ne_langmuir)
        m_eff_i_interp.append(m_eff_i)

    # Flatten all the interpolations
    Ti_interp = [val for sublist in Ti_interp for val in sublist]
    Tr_interp = [val for sublist in Tr_interp for val in sublist]
    ne_interp = [val for sublist in ne_interp for val in sublist]
    m_eff_i_interp = [val for sublist in m_eff_i_interp for val in sublist]





    # store everything in the output data dict
    data_dict_output = {**data_dict_output,
                        **{
                            'Epoch': deepcopy(data_dict_LP_current['Epoch']),
                            'Ti': [np.array(Ti_interp), {'DEPEND_0': 'Epoch', 'UNITS': 'eV'}],
                            'Te': [np.array(Ti_interp) * np.array(Tr_interp), {'DEPEND_0': 'Epoch', 'UNITS': 'eV'}],
                            'Te_DERPA1': [np.array(Te_DERPA1_langmuir), {'DEPEND_0': 'Epoch', 'UNITS': 'eV'}],
                            'Te_DERPA2': [np.array(Te_DERPA2_langmuir), {'DEPEND_0': 'Epoch', 'UNITS': 'eV'}],
                            'Tr': [np.array(Tr_interp), {'DEPEND_0': 'Epoch', 'UNITS': 'Te/Ti'}],
                            'ne': [np.array(ne_interp), {'DEPEND_0': 'Epoch', 'UNITS': 'm^-3'}],
                            'Alt': [alt_langmuir,{'DEPEND_0': 'Epoch', 'UNITS': 'km'}],
                            'm_eff_i': [np.array(m_eff_i_interp),{'DEPEND_0': 'Epoch', 'UNITS': 'kg'}],
                            'floating_potential' : [np.array(floatingPotential_langmuir),{'DEPEND_0': 'Epoch', 'UNITS': 'Volts'}],
                            'L-Shell': [np.array(LShell_langmuir),deepcopy(data_dict_LShell['L-Shell'][1])]
                           }
                        }




    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---
    if outputData:
        stl.prgMsg('Creating output file')
        fileoutName = f'ACESII_{ACESII.payload_IDs[wRocket-4]}_postFlight_cal.cdf'
        outputPath = f'C:\Data\ACESII\calibration\LP_postFlight_calibration\\{ACESII.fliers[wRocket-4]}\\' + fileoutName
        stl.outputCDFdata(outputPath, data_dict_output)
        stl.Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
LP_collect_postFlight_cal_data(wRocket)

