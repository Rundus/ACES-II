# --- Print_Payload_Coordinates.py ---
# description: plot the lat/long/alt of the rocket trajectories in order to be plotted in Google Earth

# imports
import spaceToolsLib as stl
from src.data_paths import DataPaths
import os


data_dict_attitude_high = stl.loadDictFromFile(r'C:\Data\ACESII\attitude\high\ACESII_36359_Attitude_Solution.cdf')
data_dict_attitude_low = stl.loadDictFromFile(r'C:\Data\ACESII\attitude\low\ACESII_36364_Attitude_Solution.cdf')
data_dicts = [data_dict_attitude_high, data_dict_attitude_low]

wFlyer = 1
dict = data_dicts[wFlyer]
N_downsample = 10

# downsample the data by N
Alt = dict['Alt'][0][::N_downsample]
Lat = dict['Lat'][0][::N_downsample]
Long = dict['Long'][0][::N_downsample]

# file path
file_path = rf"C:\Users\cfelt\PycharmProjects\ACESII\src\Science\GoogleEarth\KML_files\{DataPaths.fliers[wFlyer]}Flyer.txt"

if os.path.exists(file_path):
    os.remove(file_path)
    print('TEST')

with open(file_path, "w") as f:
    f.write('<coordinates>')
    for idx, val in enumerate(Alt):
        f.write(f'{Long[idx]},{Lat[idx]},{Alt[idx]}\n')
    f.write('</coordinates>')
