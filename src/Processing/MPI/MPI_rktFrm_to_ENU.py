# --- MPI_rktFrm_to_ENU.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Rotate the MPI rktFrm data into
# ENU coordinates


# imports
import numpy as np
from copy import deepcopy
import spaceToolsLib as stl
from scipy.interpolate import CubicSpline


def MPI_rktFrm_to_ENU():
    # 1. Load attitude solution (DCM)
    path_to_attitude = r'C:\Data\ACESII\attitude\low\ACESII_36364_Attitude_Solution.cdf'
    path_to_MPI = r'C:\Data\ACESII\L1\low\ACESII_36364_l1_MPI_rktFrm.cdf'
    data_dict_attitude = stl.loadDictFromFile(path_to_attitude)
    data_dict_MPI1 = stl.loadDictFromFile(path_to_MPI)

    # 2. Extract time (seconds) and DCM components
    time_epoch_attitude = data_dict_attitude[f'Epoch'][0]
    time_sec_attitude = stl.EpochTo_T0_Rocket(time_epoch_attitude, T0=time_epoch_attitude[0])
    time_epoch_MPI1 = data_dict_MPI1[f'Epoch_MPI1'][0]
    time_sec_MPI1 = stl.EpochTo_T0_Rocket(time_epoch_MPI1, T0=time_epoch_attitude[0])

    # [DCM elements]
    a11 = np.array(data_dict_attitude['a11'][0])
    a12 = np.array(data_dict_attitude['a12'][0])
    a13 = np.array(data_dict_attitude['a13'][0])
    a21 = np.array(data_dict_attitude['a21'][0])
    a22 = np.array(data_dict_attitude['a22'][0])
    a23 = np.array(data_dict_attitude['a23'][0])
    a31 = np.array(data_dict_attitude['a31'][0])
    a32 = np.array(data_dict_attitude['a32'][0])
    a33 = np.array(data_dict_attitude['a33'][0])

    # 3.  Interps for each DCM element
    dcm_a11 = CubicSpline(time_sec_attitude, a11)
    dcm_a12 = CubicSpline(time_sec_attitude, a12)
    dcm_a13 = CubicSpline(time_sec_attitude, a13)

    dcm_a21 = CubicSpline(time_sec_attitude, a21)
    dcm_a22 = CubicSpline(time_sec_attitude, a22)
    dcm_a23 = CubicSpline(time_sec_attitude, a23)

    dcm_a31 = CubicSpline(time_sec_attitude, a31)
    dcm_a32 = CubicSpline(time_sec_attitude, a32)
    dcm_a33 = CubicSpline(time_sec_attitude, a33)


    # 4. interp DCM as a fxn of time
    N = len(time_sec_MPI1)
    dcm_interp=np.zeros(shape=(N,3,3)) #for MPI 1

    for i in range(N):
        dcm_interp[i,0,0] = dcm_a11(time_sec_MPI1[i])
        dcm_interp[i,0,1] = dcm_a12(time_sec_MPI1[i])
        dcm_interp[i,0,2] = dcm_a13(time_sec_MPI1[i])
        dcm_interp[i,1,0] = dcm_a21(time_sec_MPI1[i])
        dcm_interp[i,1,1] = dcm_a22(time_sec_MPI1[i])
        dcm_interp[i,1,2] = dcm_a23(time_sec_MPI1[i])
        dcm_interp[i,2,0] = dcm_a31(time_sec_MPI1[i])
        dcm_interp[i,2,1] = dcm_a32(time_sec_MPI1[i])
        dcm_interp[i,2,2] = dcm_a33(time_sec_MPI1[i])



    # Assume velocity is stored in fields ['Vx','Vy','Vz'] and time in 'Epoch'
    mpi_vel_vec = np.array([data_dict_MPI1['Vx_rkt_MPI1'][0],data_dict_MPI1['Vy_rkt_MPI1'][0],data_dict_MPI1['Vz_rkt_MPI1'][0]]).T

    # 6. Rotate MPI vectors to ENU, dot product each mpi_vel_vec
    MPI_ENU = []
    for i in range(N):
        vec = mpi_vel_vec[i]
        matrix = dcm_interp[i]
        MPI_ENU.append(matrix @ vec)


    # 7.  output dictionary
    example_attrs = {'LABLAXIS': None, 'DEPEND_0': 'time', 'UNITS': None}

    data_dict_output = {
        'MPI_E': [np.array(MPI_ENU)[:,0], deepcopy(example_attrs)],
        'MPI_N': [np.array(MPI_ENU)[:,1], deepcopy(example_attrs)],
        'MPI_U': [np.array(MPI_ENU)[:,2], deepcopy(example_attrs)],
        'time': [np.array(time_sec_MPI1), deepcopy(example_attrs)],
        'DCM_INTERP': [np.array(dcm_interp), deepcopy(example_attrs)]
    }


    # --- 8. Output new CDF
    file_out_path = r'C:\Data\ACESII\L2\low\ACESII_35364_L2_MPI_ENU.cdf'
    stl.outputCDFdata(outputPath=file_out_path, data_dict=data_dict_output)


###EXECUTE###
MPI_rktFrm_to_ENU()
