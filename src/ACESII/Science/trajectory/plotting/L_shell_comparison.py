


import spaceToolsLib as stl
import matplotlib.pyplot as plt
import numpy as np

def compare_L_shell_data():

    data_dict_eepaa_high = stl.loadDictFromFile('C:\Data\ACESII\L2\high\ACESII_36359_l2_eepaa_fullCal.cdf')
    data_dict_eepaa_low = stl.loadDictFromFile('C:\Data\ACESII\L2\low\ACESII_36364_l2_eepaa_fullCal.cdf')

    data_dict_L_Shell_high =stl.loadDictFromFile('C:\Data\ACESII\science\L_shell\high\ACESII_36359_Lshell.cdf')
    data_dict_L_Shell_low = stl.loadDictFromFile('C:\Data\ACESII\science\L_shell\low\ACESII_36364_Lshell.cdf')

    # downsample the L-Shell data to match the eepaa data
    L_shell_high_indices = np.array([np.abs(data_dict_L_Shell_high['Epoch'][0] - val).argmin() for val in data_dict_eepaa_high['Epoch'][0]])
    L_shell_low_indices = np.array([np.abs(data_dict_L_Shell_low['Epoch'][0] - val).argmin() for val in data_dict_eepaa_low['Epoch'][0]])

    L_shells_high = data_dict_L_Shell_high['L-Shell'][0][L_shell_high_indices]
    L_shells_high = np.array([val[0] for val in L_shells_high])
    L_shells_low = data_dict_L_Shell_low['L-Shell'][0][L_shell_low_indices]
    L_shells_low = np.array([val[0] for val in L_shells_low])

    # plot everything

    fig, ax = plt.subplots(nrows=2)
    ax[0].pcolormesh(L_shells_high, data_dict_eepaa_high['Energy'][0], data_dict_eepaa_high['Differential_Number_Flux'][0][:,2,:].T, cmap='turbo', norm='log')
    ax[1].pcolormesh(L_shells_low, data_dict_eepaa_low['Energy'][0], data_dict_eepaa_low['Differential_Number_Flux'][0][:,2,:].T, cmap='turbo', norm='log')
    for i in range(2):
        ax[i].set_yscale('log')
        ax[i].set_xlim(7,10)

    plt.show()

compare_L_shell_data()