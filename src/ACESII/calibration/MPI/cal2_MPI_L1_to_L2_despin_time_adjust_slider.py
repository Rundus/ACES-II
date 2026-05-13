# --- MPI_L1_to_L2_despin_time_adjust_slider.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Rotate the MPI rktFrm data into
# ENU coordinates
import matplotlib.pyplot as plt
# imports
import numpy as np
from copy import deepcopy
import spaceToolsLib as stl
from scipy.interpolate import CubicSpline
from src.ACESII.data_tools.data_paths import DataPaths
import datetime as dt
from matplotlib.widgets import Slider


plot_interactive_despin_MPI = True


def MPI_L1_to_L2_despin_time_adjust_slider():
    # 1. Load attitude solution (DCM)
    path_to_attitude = rf'{DataPaths.ACES_data_folder}/attitude/low/ACESII_36364_Attitude_Solution.cdf'
    path_to_MPI = rf'{DataPaths.ACES_data_folder}/L1/MPI/low/ACESII_36364_l1_MPI_rktFrm.cdf'
    data_dict_attitude = stl.loadDictFromFile(path_to_attitude)
    data_dict_MPI = stl.loadDictFromFile(path_to_MPI)
    data_dict_MPI_ENU = stl.loadDictFromFile(rf'{DataPaths.ACES_data_folder}/L2/MPI/low/ACESII_36364_L2_MPI_ENU.cdf')
    data_dict_ExB = stl.loadDictFromFile('/home/connor/Data/ROCKETS/ACESII/science/ExB/low/ACESII_36364_ExB.cdf')


    if plot_interactive_despin_MPI:

        ###########################
        # --- PLOT1 ExB Compare ---
        ###########################

        fig, ax = plt.subplots(nrows=2,ncols=1)

        # --- Exb ---
        T0_ExB = stl.EpochTo_T0_Rocket(data_dict_ExB['Epoch'][0],T0=dt.datetime(2022,11,20,17,20))
        ax[0].plot(T0_ExB, data_dict_ExB['ExB_E'][0], color='tab:gray', label='ExB_E',linewidth=0.25)
        ax[1].plot(T0_ExB, data_dict_ExB['ExB_N'][0], color='tab:gray', label='ExB_N',linewidth=0.25)

        # Plot the unperturbed data
        MPI_plot_Data_E = []
        MPI_plot_Data_N = []
        for i in range(4):

            # 8. Time adjust the data to better align. This does NOT affect the Despin at all
            spin_rate = 0.5474
            if i == 0:  # MPI1
                delta = (1 / spin_rate) * (2 + 0.75)
            elif i == 1:  # MPI2
                delta = (1 / spin_rate) * (1 + 0.5)
            elif i == 2:  # MPI3
                delta = (1 / spin_rate) * (1 + 0.25)
            elif i == 3:  # MPI4
                delta = 0
            time_sec_MPI = data_dict_MPI_ENU[f'time{i+1}'][0] - delta

            lineMPI_E, = ax[0].plot(time_sec_MPI,data_dict_MPI_ENU[f'MPI{i+1}_E'][0],label=f'MPI{i+1}_E')
            lineMPI_N, =ax[1].plot(time_sec_MPI, data_dict_MPI_ENU[f'MPI{i+1}_N'][0],label=f'MPI{i+1}_N')
            MPI_plot_Data_E.append(lineMPI_E)
            MPI_plot_Data_N.append(lineMPI_N)



        axSlope = plt.axes([0.25, 0.05, 0.65, 0.03])
        slopeSlider = Slider(
            axSlope, 'Slope', valmin=-0.01, valmax=0.01, valinit=-0.00073
        )

        axInter = plt.axes([0.25, 0.025, 0.65, 0.03])
        interSlider = Slider(
            axInter, 'Intercept', valmin=-0.5, valmax=0.5, valinit=-0.058
        )

        # Make some adjustments
        for i in range(2):
            ax[i].set_ylim(-500,500)
            ax[i].set_xlim(230,420)
            ax[i].legend()

        # --- UPDATE FUNCTION ---

        def update(val):
            slopeVal = slopeSlider.val
            interVal = interSlider.val


            T0 = dt.datetime(2022,11,20, 17,20)
            T0_attitude = stl.EpochTo_T0_Rocket(data_dict_attitude[f'Epoch'][0],T0=T0)
            T0_attitude_new = (1 + slopeVal) * T0_attitude + interVal

            # FORM the DCM Cubic Splines matrix
            DCM_splines = [
                    [[], [], []],
                    [[], [], []],
                    [[], [], []]
            ]

            for i in range(3):
                for j in range(3):
                    cs = CubicSpline(T0_attitude_new,data_dict_attitude[f'a{i+1}{j+1}'][0])
                    DCM_splines[i][j] = cs

            # Loop over each MPI
            for idx in range(4):
                time_epoch_MPI = data_dict_MPI[f'Epoch_MPI{idx+1}'][0]
                time_sec_MPI = stl.EpochTo_T0_Rocket(time_epoch_MPI, T0=T0)

                # 4. interp DCM as a fxn of time
                dcm_interp=np.zeros(shape=(len(time_sec_MPI),3,3)) #for MPI 1

                # Loop over each timepoint in the MPI data
                for i in range(3):
                    for j in range(3):
                        dcm_interp[:,i,j] = DCM_splines[i][j](time_sec_MPI)

                # Assume velocity is stored in fields ['Vx','Vy','Vz'] and time in 'Epoch'
                MPI_RKT = np.array([data_dict_MPI[f'Vx_rkt_MPI{idx+1}'][0],data_dict_MPI[f'Vy_rkt_MPI{idx+1}'][0],data_dict_MPI[f'Vz_rkt_MPI{idx+1}'][0]]).T

                # 6. Rotate MPI vectors to ENU, dot product each mpi_vel_vec
                MPI_ENU = np.array([np.matmul(dcm_interp[i],MPI_RKT[i]) for i in range(len(MPI_RKT))])

                # Set the new data
                MPI_plot_Data_E[idx].set_ydata(MPI_ENU[:, 0])
                MPI_plot_Data_N[idx].set_ydata(MPI_ENU[:, 1])

        slopeSlider.on_changed(update)
        interSlider.on_changed(update)
        plt.show()





###EXECUTE###
MPI_L1_to_L2_despin_time_adjust_slider()
