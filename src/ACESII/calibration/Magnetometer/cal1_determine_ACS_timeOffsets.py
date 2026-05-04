# --- cal1_determine_ACS_timeOFfsets.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: adjust the timebase of the ACS DCM to see if a cleaner MAG despin can be derived. Also remove any NaNs from the data

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2026-03-30"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np

# --- --- --- --- ---


######################
# --- DATA TOGGLES ---
######################
just_print_file_names_bool = False
rocket_str = 'high'
wInstr = 'MAG'
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'{wInstr}':['L1', [[0],[0]]],
    'attitude':['',[[0],[0]]],
    '/':['science/ESA_currents/EEPAA',[[0],[0]]],
    '':['science/ESA_currents/IEPAA',[[0],[0]]]
}
outputData = False

# slope_init = {'high':0.0000185,
#               'low':0.000089}

slope_init = {'high':-0.000006,
              'low':0.000089}
# intercept_init = {'high':0.12756724137931032,
#                   'low':0.12078947368421052}
intercept_init = {'high':0.12756724137931032,
                  'low':0.12756724137931032}

plot_interactive_slider = True

#################
# --- IMPORTS ---
#################
from src.ACESII.data_tools.my_imports import *
import datetime as dt
import ppigrf
from matplotlib.widgets import Slider
from geopy import distance
start_time = time.time()

def cal1_determine_ACS_timeOFfsets(data_dicts):

    #######################
    # --- LOAD THE DATA ---
    #######################
    stl.prgMsg(f'Loading data')
    data_dict_MAG = deepcopy(data_dicts[0])
    data_dict_attitude = deepcopy(data_dicts[1])
    data_dict_Jpara_EEPAA = deepcopy(data_dicts[2])
    data_dict_Jpara_IEPAA = deepcopy(data_dicts[3])
    stl.Done(start_time)

    # --- prepare the output ---
    data_dict_output = {
                       }
    data_dict_raw = {
        f'a{i}{j}':[[],{}] for i in range(1,4) for j in range(1,4)
    }
    data_dict_raw = {**data_dict_raw,
                     **{'Lat':[[],{}],
                        'Alt':[[],{}],
                        'Long' :[[],{}],}
                     }

    # --- --- --- --- --- --- ---
    # --- ADJUST THE ACS TIME ---
    # --- --- --- --- --- --- ---

    # --- [0] Calculate the time lag due to speed of light ---
    # The IGNORES any curvature due to the Earth, which might skew results
    Alat = 69.294167
    Alon = 16.020833
    Aalt = 12.7  # in meters
    pt1 = np.array([[Alat, Alon, Aalt] for i in range(len(data_dict_attitude['Epoch'][0]))])
    pt2 = np.array([data_dict_attitude['Lat'][0], data_dict_attitude['Long'][0], data_dict_attitude['Alt'][0]]).T
    distance_2d = np.array([distance.distance(pt1[i][:2], pt2[i][:2]).m for i in range(len(data_dict_attitude['Epoch'][0]))])
    distance_3d = np.sqrt(np.square(distance_2d) + np.square(pt1[:, 2] - pt2[:, 2])) / stl.m_to_km
    light_time_delay = distance_3d * stl.m_to_km / stl.lightSpeed  # in seconds
    # fig, ax = plt.subplots(2)
    # ax[0].plot(data_dict_attitude['Epoch'][0], distance_3d)
    # ax[0].set_ylabel('Distance [km]')
    # ax[1].plot(data_dict_attitude['Epoch'][0], light_time_delay)
    # ax[1].set_ylabel('Light Time Delay [s]')
    # plt.show()


    # === [1] Interpolate the DCM directly (no adjustments) ===
    Bkeys = ['Bx', 'By', 'Bz']
    Bkeys_new = ['B_E', 'B_N', 'B_U']
    MAGdata = np.array([data_dict_MAG[key][0] for key in Bkeys]).T
    T0 = dt.datetime(2022,11,20,17,20)
    T0_MAG = stl.EpochTo_T0_Rocket(data_dict_MAG['Epoch'][0],T0=T0)
    T0_attitude = stl.EpochTo_T0_Rocket(data_dict_attitude['Epoch'][0], T0=T0) - light_time_delay

    for i in range(1,4):
        for j in range(1,4):
            DCMkey = f'a{i}{j}'
            data_dict_raw[DCMkey][0] = np.interp(T0_MAG, T0_attitude,data_dict_attitude[DCMkey][0])

    for akey in ['Lat','Long','Alt']:
        data_dict_raw[akey][0] = np.interp(T0_MAG,T0_attitude,data_dict_attitude[akey][0])
        # data_dict_raw[akey][0] = np.interp(T0_MAG, T0_attitude, data_dict_traj[akey][0])

    # --- Apply the DCM directly - no time adjustments ---
    DCM_raw = np.array([
        [
            [data_dict_raw['a11'][0][i],data_dict_raw['a12'][0][i],data_dict_raw['a13'][0][i]],
            [data_dict_raw['a21'][0][i],data_dict_raw['a22'][0][i],data_dict_raw['a23'][0][i]],
            [data_dict_raw['a31'][0][i],data_dict_raw['a32'][0][i],data_dict_raw['a33'][0][i]],
        ]
        for i in range(len(data_dict_MAG['Epoch'][0]))
    ])

    # apply the raw DCM
    B_ENU = np.array([np.matmul(DCM_raw[i],MAGdata[i]) for i in range(len(data_dict_MAG['Epoch'][0]))])


    # Calculate the deltaB
    dB_E_raw = np.diff(B_ENU[:, 0], prepend=B_ENU[:, 0][0])
    dB_N_raw = np.diff(B_ENU[:, 1], prepend=B_ENU[:, 1][0])
    dB_U_raw = np.diff(B_ENU[:, 2], prepend=B_ENU[:, 2][0])
    dB_ENU = np.array([dB_E_raw,dB_N_raw,dB_U_raw])

    # Filter the deltaB
    dB_filtered = [[],[],[]]
    for i in range(3):
        filtData = dB_ENU[i]
        filtData[np.isnan(filtData)] = 0
        dB_filtered[i] = stl.butterFilter().butter_filter(data=filtData,
                                                       lowcutoff=0.05,
                                                       highcutoff=0.05,
                                                       order=2,
                                                       fs=256,
                                                       filtertype='lowPass'
                                                       )


    # Calculate the CHAOS magnetic field for comparison
    stl.prgMsg('Calculating CHAOS')
    B_MODEL = stl.CHAOS(data_dict_raw['Lat'][0], data_dict_raw['Long'][0], data_dict_raw['Alt'][0] / stl.m_to_km, data_dict_MAG['Epoch'][0])
    stl.Done(start_time)

    # Calculate the IGRF magnetic field for comparison
    # stl.prgMsg('Calculating IGRF')
    # date = dt.datetime(2022,11,20,17,25)
    # Be,Bn,Bu = ppigrf.igrf(data_dict_raw['Long'][0],data_dict_raw['Lat'][0],data_dict_raw['Alt'][0]/stl.m_to_km,date)
    # B_MODEL = np.array([Be[0],Bn[0],Bu[0]]).T
    # stl.Done(start_time)

    if plot_interactive_slider:
        # === Plot the Despin (Raw) ===
        fig, ax = plt.subplots(nrows=4,ncols=3)

        if rocket_str == 'high':
            xLimits = [dt.datetime(2022, 11, 20, 17, 23),
                       dt.datetime(2022, 11, 20, 17, 28)]
        else:
            xLimits = [dt.datetime(2022, 11, 20, 17, 23,50),
                       dt.datetime(2022, 11, 20, 17, 28)]


        # --- Plot the J_parallel from ESAs ---
        for i in range(3):
            eData = data_dict_Jpara_EEPAA['j_para'][0] / (1E-6)
            iData = data_dict_Jpara_IEPAA['j_para'][0] / (1E-6)

            iDta_filtered = stl.butterFilter().butter_filter(data=iData,
                                                             lowcutoff=0.15,
                                                             highcutoff=0.15,
                                                             order=4,
                                                             fs=20,
                                                             filtertype='lowPass'
                                                             )

            ax[0,i].plot(data_dict_Jpara_EEPAA['Epoch'][0], eData, color='tab:blue', label=r'$J_{\parallel}^{e-}$')
            ax[0,i].plot(data_dict_Jpara_IEPAA['Epoch'][0], -1 * iDta_filtered, color='tab:red',
                          label=r'$J_{\parallel}^{i+}$')
            ax[0,i].set_ylim(-10, 10)
            ax[0,i].set_ylabel(r'Current [$\mu$A]')
            ax[0,i].set_xlim(*xLimits)
            ax[0,i].legend(loc='upper right')


        # === FIRST COLUMN ===
        line1= ax[1,0].plot(data_dict_MAG['Epoch'][0], B_ENU[:, 0], label=f'{Bkeys_new[0]}_raw')
        line2 = ax[2,0].plot(data_dict_MAG['Epoch'][0], B_ENU[:, 1], label=f'{Bkeys_new[1]}_raw')
        line3 = ax[3,0].plot(data_dict_MAG['Epoch'][0], B_ENU[:, 2], label=f'{Bkeys_new[2]}_raw')

        ax[1,0].plot(data_dict_MAG['Epoch'][0], B_MODEL[:, 0], label=f'{Bkeys_new[0]}_MODEL')
        ax[2,0].plot(data_dict_MAG['Epoch'][0], B_MODEL[:, 1], label=f'{Bkeys_new[1]}_MODEL')
        ax[3,0].plot(data_dict_MAG['Epoch'][0], B_MODEL[:, 2], label=f'{Bkeys_new[2]}_MODEL')

        line1_adjust, = ax[1,0].plot(data_dict_MAG['Epoch'][0], B_ENU[:, 0],color='tab:red')
        line2_adjust, = ax[2,0].plot(data_dict_MAG['Epoch'][0], B_ENU[:, 1],color='tab:red')
        line3_adjust, = ax[3,0].plot(data_dict_MAG['Epoch'][0], B_ENU[:, 2],color='tab:red')

        ax[1, 0].set_ylim(-2500, 6000)
        ax[2, 0].set_ylim(7000, 16000)
        ax[3, 0].set_ylim(-50000, -42000)

        # === SECOND COLUMN ===
        ax[1, 1].plot(data_dict_MAG['Epoch'][0], B_MODEL[:, 0] - B_ENU[:, 0], label='CHAOS - BMAG')
        ax[2, 1].plot(data_dict_MAG['Epoch'][0], B_MODEL[:, 1] - B_ENU[:, 1], label='CHAOS - BMAG')
        ax[3, 1].plot(data_dict_MAG['Epoch'][0], B_MODEL[:, 2] - B_ENU[:, 2], label='CHAOS - BMAG')
        ax[1,1].axhline(y=0,color='red',alpha=0.5,linestyle='--')
        ax[2, 1].axhline(y=0, color='red', alpha=0.5, linestyle='--')
        ax[3, 1].axhline(y=0, color='red', alpha=0.5, linestyle='--')

        line1_sub_adjust, = ax[1, 1].plot(data_dict_MAG['Epoch'][0], B_MODEL[:,0]-B_ENU[:, 0],color='tab:red',label='CHAOS - BMAG')
        line2_sub_adjust, = ax[2, 1].plot(data_dict_MAG['Epoch'][0], B_MODEL[:,1]-B_ENU[:, 1], color='tab:red',label='CHAOS - BMAG')
        line3_sub_adjust, = ax[3, 1].plot(data_dict_MAG['Epoch'][0], B_MODEL[:,2]-B_ENU[:, 2], color='tab:red',label='CHAOS - BMAG')

        ax[1, 1].set_ylim(-1000, 1000)
        ax[2, 1].set_ylim(-500, 500)
        ax[3, 1].set_ylim(-250, 250)

        # === THIRD COLUMN ===

        line1_db_E = ax[1, 2].plot(data_dict_MAG['Epoch'][0],dB_filtered[0],color='tab:blue',label='dB$_{E}$ (raw)')
        line2_db_N = ax[2, 2].plot(data_dict_MAG['Epoch'][0], dB_filtered[1], color='tab:blue', label='dB$_{N}$ (raw)')
        line3_db_U = ax[3, 2].plot(data_dict_MAG['Epoch'][0], dB_filtered[2], color='tab:blue', label='dB$_{U}$ (raw)')

        line1_db_E_adjust, = ax[1, 2].plot(data_dict_MAG['Epoch'][0], dB_filtered[0], color='tab:red', label='dB$_{E}$')
        line2_db_N_adjust, = ax[2, 2].plot(data_dict_MAG['Epoch'][0], dB_filtered[1], color='tab:red', label='dB$_{N}$')
        line3_db_U_adjust, = ax[3, 2].plot(data_dict_MAG['Epoch'][0], dB_filtered[2], color='tab:red', label='dB$_{U}$')


        # === [2] Adjust the DCM interpolation time (with slider) ===

        axSlope = plt.axes([0.25,0.05,0.65,0.03])
        slopeSlider = Slider(
            axSlope,'Slope',valmin=-0.0001,valmax=0.0001,valinit=slope_init[rocket_str]
        )

        axnonlinearSlope = plt.axes([0.25, 0.07, 0.65, 0.03])
        nonlinearSlopeSlider = Slider(
            axnonlinearSlope, 'Nonlinear-Slope', valmin=-0.000001, valmax=0.000001, valinit=0
        )

        axInter = plt.axes([0.25, 0.025, 0.65, 0.03])
        interSlider = Slider(
            axInter, 'Intercept', valmin=-0.15, valmax=0.15, valinit=intercept_init[rocket_str]
        )


        def update(val):
            nonLinearSlopeVal = nonlinearSlopeSlider.val
            slopeVal = slopeSlider.val
            interVal = interSlider.val

            # adjust the attitude timebase
            # T0_attitude_new = (nonLinearSlopeVal)*np.square(T0_attitude)+(1+slopeVal)*T0_attitude + interVal
            T0_attitude_new =  (1 + slopeVal) * T0_attitude + interVal

            # re-interpolate the DCM
            DCM_adjust = np.zeros(shape=(len(data_dict_MAG['Epoch'][0]), 3, 3))
            for i in range(3):
                for j in range(3):
                    DCM_adjust[:, i, j] = np.interp(T0_MAG,T0_attitude_new, data_dict_attitude[f'a{i + 1}{j + 1}'][0])

            # apply the time-adjusted DCM to the MAG data
            MAGdata_adjust = np.array([np.matmul(DCM_adjust[i], MAGdata[i]) for i in range(len(T0_MAG))])

            # Calculate the dB
            db_filt_adjust = [[],[],[]]
            for i in range(3):
                filtData = MAGdata_adjust[:,i]
                filtData[np.isnan(filtData)]=0
                db = np.diff(filtData,prepend=filtData[0])
                db_filt_adjust[i] = stl.butterFilter().butter_filter(data=db,
                                                                  lowcutoff=0.05,
                                                                  highcutoff=0.05,
                                                                  order=4,
                                                                  fs=256,
                                                                  filtertype='lowPass'
                                                                  )

            # Update the plot Data
            line1_adjust.set_ydata(MAGdata_adjust[:, 0])
            line2_adjust.set_ydata(MAGdata_adjust[:, 1])
            line3_adjust.set_ydata(MAGdata_adjust[:, 2])

            line1_sub_adjust.set_ydata(B_MODEL[:,0]-MAGdata_adjust[:, 0])
            line2_sub_adjust.set_ydata(B_MODEL[:,1]-MAGdata_adjust[:, 1])
            line3_sub_adjust.set_ydata(B_MODEL[:,2]-MAGdata_adjust[:, 2])

            line1_db_E_adjust.set_ydata(db_filt_adjust[0])
            line2_db_N_adjust.set_ydata(db_filt_adjust[1])
            line3_db_U_adjust.set_ydata(db_filt_adjust[2])

            fig.canvas.draw_idle()

        slopeSlider.on_changed(update)
        interSlider.on_changed(update)
        # nonlinearSlopeSlider.on_changed(update)

        # Make some adjustments
        for i in range(4):
            if i >0:
                ax[i,2].legend()
                ax[i,2].set_ylim(-0.3,0.3)
                ax[i,2].axhline(y=0,linestyle='--',color='tab:red',alpha=0.5)
        for i in range(1,4):
            ax[i, 0].legend()
            ax[i, 1].legend()

        for i in range(2):
            for j in range(3):
                ax[i,j].set_xticklabels([])

        for i in range(4):
            for j in range(3):
                ax[i,j].set_xlim(*xLimits)

        plt.show()

    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')


        # Apply the initial slope/offset to the DCM, then apply it
        T0_attitude_new = (1 + slope_init[rocket_str]) * T0_attitude + intercept_init[rocket_str]
        DCM_adjust = np.zeros(shape=(len(data_dict_MAG['Epoch'][0]), 3, 3))
        for i in range(3):
            for j in range(3):
                DCM_adjust[:, i, j] = np.interp(T0_MAG, T0_attitude_new, data_dict_attitude[f'a{i + 1}{j + 1}'][0])

        MAGdata_adjust = np.array([np.matmul(DCM_adjust[i], MAGdata[i]) for i in range(len(T0_MAG))])

        # === Construct the output data dict ===

        for key in ['Epoch','Lat','Long','Alt','L-Shell','mLat','mLong','ILat','ILong']:
            if key in data_dict_MAG.keys():

                data_dict_output = {
                    **data_dict_output,
                    **{key:deepcopy(data_dict_MAG[key])}
                }

        # remove nans in the data
        bad_idxs = []
        for i in range(3):
            bad_idxs.append(np.where(np.isnan(MAGdata_adjust[:,i]))[0])
        bad_idxs = np.array(list(set(np.array(bad_idxs).flatten()))) # get only the unique data in an array


        MAG_data_clean = [[],[],[]]
        B_MODEL_clean = [[],[],[]]
        if bad_idxs.tolist() != []:
            data_dict_output['Epoch'][0] = np.delete(data_dict_MAG['Epoch'][0],bad_idxs)
            for i in range(3):
                MAG_data_clean[i] = np.delete(MAGdata_adjust[:,i],bad_idxs)
                B_MODEL_clean[i] = np.delete(B_MODEL[:,i],bad_idxs)
        else:
            MAG_data_clean = MAGdata_adjust.T
            B_MODEL_clean = B_MODEL.T


        data_dict_output = {
            **data_dict_output,
            **{
                'B_model_E': [B_MODEL_clean[0],{'LABLAXIS':'B_model_East','UNITS':'nT','VAR_TYPE':'data'}],
                'B_model_N': [B_MODEL_clean[1],{'LABLAXIS':'B_model_North','UNITS':'nT','VAR_TYPE':'data'}],
                'B_model_U': [B_MODEL_clean[2],{'LABLAXIS':'B_model_Up','UNITS':'nT','VAR_TYPE':'data'}],
                'B_E': [MAG_data_clean[0], {'LABLAXIS':'B_East','UNITS':'nT','VAR_TYPE':'data'}],
                'B_N': [MAG_data_clean[1], {'LABLAXIS':'B_North','UNITS':'nT','VAR_TYPE':'data'}],
                'B_U': [MAG_data_clean[2], {'LABLAXIS':'B_Up','UNITS':'nT','VAR_TYPE':'data'}],
            }
        }

        # write out the data
        fileoutName = f'ACESII_{ACESII.fliers_dict[rocket_str]}_l1_RingCore_cal1_DCM_time_adjust.cdf'
        outputPath = f'{DataPaths.ACES_data_folder}/calibration/{wInstr}/{rocket_str}//{fileoutName}'
        stl.outputDataDict(outputPath, data_dict=data_dict_output)
        stl.Done(start_time)




#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(cal1_determine_ACS_timeOFfsets,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)

