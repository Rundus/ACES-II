import matplotlib.pyplot as plt
from src.my_imports import *

# TOGGLES
wRocket = 5

# get the data
data_path = glob(rf'C:\Data\ACESII\science\payload_potential\{ACESII.fliers[wRocket - 4]}\\*.cdf*')[0]
data_dict_payload_potential = stl.loadDictFromFile(data_path)


for sweep_no in [81,82,83,84,85,86,87,88]:
    # get the sweep-specific data
    Epoch_val = data_dict_payload_potential['Epoch'][0][sweep_no]
    steps_data = data_dict_payload_potential['sweeps_voltage'][0][sweep_no]
    current_data = data_dict_payload_potential['sweeps_current'][0][sweep_no]

    # Plot everything
    fig, ax = plt.subplots()
    ax.scatter(steps_data,current_data)
    ax.set_ylabel('Current [nA]')
    ax.set_xlabel('Voltage [V]')
    fig.suptitle(f'Sweep No. {sweep_no}\n {Epoch_val} UTC')
    plt.show()

