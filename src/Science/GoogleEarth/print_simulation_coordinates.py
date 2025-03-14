# --- Print_Payload_Coordinates.py ---
# description: plot the lat/long/alt of the rocket trajectories in order to be plotted in Google Earth

# imports
import spaceToolsLib as stl
from src.data_paths import DataPaths
import os


data_dict_spatial = stl.loadDictFromFile(r'C:\Data\physicsModels\ionosphere\spatial_environment\spatial_environment.cdf')


# downsample the data by N
simAlt = data_dict_spatial['simAlt'][0]
simLShell = data_dict_spatial['simLShell'][0]
grid_lat = data_dict_spatial['grid_lat'][0]
grid_alt = data_dict_spatial['grid_alt'][0]
grid_long = data_dict_spatial['grid_long'][0]


# file path
file_path = rf"C:\Users\cfelt\PycharmProjects\ACESII\src\Science\GoogleEarth\KML_files\simulation.txt"

if os.path.exists(file_path):
    os.remove(file_path)

with open(file_path, "w") as f:
    f.write('<coordinates>')

    for idx1, val1 in enumerate(simLShell):
        for idx2, val2 in enumerate(simAlt):
            f.write(f'{grid_long[idx1][idx2]},{grid_lat[idx1][idx2]},{grid_alt[idx1][idx2]}\n')

    f.write('</coordinates>')
