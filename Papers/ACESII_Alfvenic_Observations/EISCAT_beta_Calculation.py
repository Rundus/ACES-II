# --- EISCAT_beta_Calculation.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Calculate the electron and ion beta at various altitudes using the EISCAT data
# and CHAOS model
import numpy as np

from ACESII_code.myImports import *
from ACESII_code.class_var_func import kB, q0,Bdip_mag

targetTime = [dt.datetime(2022, 11, 20, 17, 23, 0, 000000),
              dt.datetime(2022, 11, 20, 17, 29, 0, 000000)]
wRange = [400, 700] # in kilometers
fileName = r"C:\Data\ACESII\science\EISCAT\tromso\UHF\MAD6400_2022-11-20_beata_ant@uhfa.cdf"

data_dict = loadDictFromFile(inputFilePath=fileName)

low_idx, high_idx = np.abs(data_dict['Epoch'][0] - targetTime[0]).argmin(), np.abs(data_dict['Epoch'][0] - targetTime[1]).argmin()

low_idx_alt, high_idx_alt = np.abs(data_dict['range'][0] - wRange[0]).argmin(), np.abs(data_dict['range'][0] - wRange[1]).argmin()

tr = data_dict['tr'][0][low_idx:high_idx, low_idx_alt:high_idx_alt+1]
ti = data_dict['ti'][0][low_idx:high_idx, low_idx_alt:high_idx_alt+1]
ne = data_dict['ne'][0][low_idx:high_idx, low_idx_alt:high_idx_alt+1]
te = ti*tr

TempI = [[val for val in arr if not np.isnan(val)] for arr in ti]
TempE = [[val for val in arr if not np.isnan(val)] for arr in te]
Ne = [[val for val in arr if not np.isnan(val)] for arr in ne]
AvgNe = np.average([np.average(arr) for arr in Ne])
AvgT_ions = kB*np.average([np.average(arr) for arr in TempI])
AvgT_elec = kB*np.average([np.average(arr) for arr in TempE])
Bdip = np.average(Bdip_mag(Alt_km= [i for i in range(wRange[0], wRange[1], 10)], Lat_deg= [70 for i in range(10)]))

# --- Calculate Plasma Beta ---
beta_Ti = 2*u0*AvgNe*AvgT_ions / (Bdip**2)
beta_Te = 2*u0*AvgNe*AvgT_elec / (Bdip**2)

print(r'$\beta_{i}$ = ' + f'{beta_Ti}')
print(r'$\beta_{e}$ = ' + f'{beta_Te}')




