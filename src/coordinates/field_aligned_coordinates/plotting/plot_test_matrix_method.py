# Description: Plot the old RingCore_Field_Aligned data against the new method
import numpy as np

from src.my_imports import *
import spaceToolsLib as stl
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

data_dict_mag_FAC_old =stl.loadDictFromFile('C:\Data\ACESII\L2\high\ACESII_36359_l2_RingCore_Field_Aligned.cdf')
data_dict_mag_ENU = stl.loadDictFromFile('C:\Data\ACESII\L2\high\ACESII_36359_l2_RingCore_ENU.cdf')
data_dict_transform_ENU = stl.loadDictFromFile('C:\Data\ACESII\coordinates\high\ACESII_36359_ECEF_to_ENU.cdf')
data_dict_transform_FAC = stl.loadDictFromFile('C:\Data\ACESII\coordinates\high\ACESII_36359_ECEF_to_FAC.cdf')

# Interpolate transformation matricies onto magnetometer timebase
Epoch_mag_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_mag_ENU['Epoch'][0]])
Epoch_trans_tt2000 = np.array([pycdf.lib.datetime_to_tt2000(val) for val in data_dict_transform_ENU['Epoch'][0]])

for key, val in data_dict_transform_ENU.items():
    if key != 'Epoch':
        cs = CubicSpline(Epoch_trans_tt2000,data_dict_transform_ENU[key][0])
        data_dict_transform_ENU[key][0] = cs(Epoch_mag_tt2000)

for key, val in data_dict_transform_FAC.items():
    if key != 'Epoch':
        cs = CubicSpline(Epoch_trans_tt2000,data_dict_transform_FAC[key][0])
        data_dict_transform_FAC[key][0] = cs(Epoch_mag_tt2000)


# form the vectors
B_ENU = np.array([data_dict_mag_ENU['B_East'][0], data_dict_mag_ENU['B_North'][0], data_dict_mag_ENU['B_Up'][0]]).T
B_FAC = np.array([data_dict_mag_FAC_old['B_r'][0], data_dict_mag_FAC_old['B_e'][0], data_dict_mag_FAC_old['B_p'][0]]).T

# apply the transform
B_FAC_new = np.array([np.matmul(data_dict_transform_FAC['ECEF_to_FAC'][0][i],np.matmul( data_dict_transform_ENU['ENU_to_ECEF'][0][i],vec)) for i,vec in enumerate(B_ENU)])

# plot the results to compare
fig, ax = plt.subplots(3)
# ax[0].scatter(data_dict_mag_FAC_old['Epoch'][0],data_dict_mag_FAC_old['B_r'][0],color='tab:blue',s=20)
# ax[0].scatter(data_dict_mag_ENU['Epoch'][0],B_FAC_new[:,0],color='tab:red',s=10)
ax[0].plot(data_dict_mag_ENU['Epoch'][0], data_dict_mag_FAC_old['B_r'][0]-B_FAC_new[:,0])

# ax[1].scatter(data_dict_mag_FAC_old['Epoch'][0],data_dict_mag_FAC_old['B_e'][0],color='tab:blue',s=20)
# ax[1].scatter(data_dict_mag_ENU['Epoch'][0],B_FAC_new[:,1],color='tab:red',s=10)
ax[1].plot(data_dict_mag_ENU['Epoch'][0], data_dict_mag_FAC_old['B_e'][0]-B_FAC_new[:,1])

# ax[2].scatter(data_dict_mag_FAC_old['Epoch'][0],data_dict_mag_FAC_old['B_p'][0],color='tab:blue',s=20)
# ax[2].scatter(data_dict_mag_ENU['Epoch'][0],B_FAC_new[:,2],color='tab:red',s=10)
ax[2].plot(data_dict_mag_ENU['Epoch'][0], data_dict_mag_FAC_old['B_p'][0]-B_FAC_new[:,2])

for i in range(3):
    ax[i].set_ylim(-1.5,1.5)

plt.show()

