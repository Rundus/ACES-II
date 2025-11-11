# --- Plot3_simualtion_extent.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Plot the simulation extent + EISCAT Radar Overlap

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import numpy as np

from src.my_imports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import spaceToolsLib as stl
import matplotlib.pyplot as plt
from spacepy import coordinates as coord
from spacepy.time import Ticktock

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
timeTargetsUTC = [dt.datetime(2022,11,20,17,23,20,100000),
                      dt.datetime(2022,11,20,17,24,00,100000),
                      dt.datetime(2022,11,20,17,24,40,100000),
                      dt.datetime(2022,11,20,17,25,20,100000),
                      dt.datetime(2022,11,20,17,26,00,100000),
                      dt.datetime(2022,11,20,17,26,40,100000),
                      dt.datetime(2022,11,20,17,27,20,100000)] # find the UTC dates times of the specifically sampled labels
trajColors = ['tab:red', 'tab:orange']

# --- --- --- --- --- --- -
# --- LOAD ALL THE DATA ---
# --- --- --- --- --- --- -
print(stl.color.UNDERLINE + f'Plot1_AllSky' + stl.color.END)
stl.prgMsg(f'Loading ACESII traj data')

# trajectory
attitudeFolderPath = f'{DataPaths.ACES_data_folder}\\attitude'
inputFiles = [glob(f'{attitudeFolderPath}\\{ACESII.fliers[0]}\\*fullCal.cdf*')[0], glob(f'{attitudeFolderPath}\\{ACESII.fliers[1]}\\*fullCal.cdf*')[0]]
data_dicts_attitude = [stl.loadDictFromFile(inputFiles[0]),stl.loadDictFromFile(inputFiles[1])]
data_dict_EISCAT = stl.loadDictFromFile(r'C:\Data\ACESII\science\EISCAT\tromso\UHF\\MAD6400_2022-11-20_beata_ant@uhfa.cdf')
data_dict_simulation = stl.loadDictFromFile('C:\Data\physicsModels\ionosphere\spatial_environment\spatial_environment.cdf')

# define some variables
EpochRocket = [data_dicts_attitude[0]['Epoch'][0], data_dicts_attitude[1]['Epoch'][0]]
geoAlt = [data_dicts_attitude[0]['Alt'][0], data_dicts_attitude[1]['Alt'][0]]
geoLat = [data_dicts_attitude[0]['Lat'][0], data_dicts_attitude[1]['Lat'][0]]
geoLong = [data_dicts_attitude[0]['Long'][0], data_dicts_attitude[1]['Long'][0]]
geoMagLat = [data_dicts_attitude[0]['Lat_geom'][0], data_dicts_attitude[1]['Lat_geom'][0]]
stl.Done(start_time)


############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

# --- --- --- --- ---
# --- Plot plot ---
# --- --- --- --- ---
plt.style.use(r'C:\Users\cfelt\PycharmProjects\ACESII\src\Papers\ACESII_Joule_Heating\Plot3\Plot3_mpl_style_sheet.mplstyle')
fig, ax = plt.subplots()
ax.set_ylabel('Altitude [km]')
ax.set_xlabel('Distance from Launch (+N/-S) [km]')
ax.set_ylim(0, 430)
ax.set_xlim(-20, 600)
ax.minorticks_on()

# plot the pseudo geomagnetic field line
N = 20 # number of L-Shell lines to plot
LShellRange = np.linspace(6.5,11.5,N)
alt_range = np.linspace(0,430,N)*1000 # in meters
LongGeomRange = np.linspace(111.83,116.82,N)
for idx, Lval in tqdm(enumerate(LShellRange)):

    # get the geomagnetic coordinate of the P.O.I. based on L-Shell
    geomagAlts = [((alt + stl.Re * stl.m_to_km) / (stl.Re * stl.m_to_km)) for alt in alt_range]
    geomagLats = np.array([np.degrees(np.arccos(np.sqrt(radi / Lval))) for radi in geomagAlts])
    geomagLongs = np.array([LongGeomRange[idx] for i in range(len(alt_range))])
    times = [dt.datetime(2022,11,20,17,20) for i in range(len(alt_range))]

    # Convert to geographic coordinates
    Pos = np.array([geomagAlts, geomagLats, geomagLongs]).transpose()
    ISOtime = [times[i].isoformat() for i in range(len(times))]
    cvals_MAG = coord.Coords(Pos, 'MAG', 'sph')
    cvals_MAG.ticks = Ticktock(ISOtime, 'ISO')
    cvals_GDZ = cvals_MAG.convert('GEO', 'sph')

    # store the data for a particular B-Field line
    B_line_lat = (np.array(cvals_GDZ.lati) - geoLat[0][0]) * 111
    B_line_alt = alt_range / stl.m_to_km

    if idx == len(LShellRange)-1:
        ax.plot(B_line_lat, B_line_alt, color='tab:blue', linestyle='-.', alpha=0.2, label='B$_{Geo}$')
    else:
        ax.plot(B_line_lat, B_line_alt, color='tab:blue', linestyle='-.', alpha=0.2)


# plot the trajectory over everything
geoLat_km = [(geoLat[0]-geoLat[0][0])*111,(geoLat[1]-geoLat[1][0])*111]
ax.plot(geoLat_km[0], geoAlt[0]/1000, color=trajColors[0], label='ACES-II 36.359')  # High
ax.plot(geoLat_km[1], geoAlt[1]/1000, color=trajColors[1], label='ACES-II 36.364')  # Low

# --- Plot the EISCAT RADAR LINES ---
tromso_S0 = [(69.586297-geoLat[0][0])*111, 19.224617] # lat/long
start_idx = np.abs(data_dict_EISCAT['Epoch'][0] - dt.datetime(2022,11,20,17,20)).argmin()
end_idx = np.abs(data_dict_EISCAT['Epoch'][0] - dt.datetime(2022,11,20,17,30)).argmin()
azmuths = data_dict_EISCAT['azm'][0][start_idx:end_idx]-360
elevs = data_dict_EISCAT['elm'][0][start_idx:end_idx]
x_proj = np.sin(np.radians(90-elevs))*np.cos(np.radians(azmuths))
z_proj = np.cos(np.radians(90-elevs))
slopes = z_proj/x_proj
for i in range(len(azmuths)):
    if i ==0:
        ax.axline(xy1=(tromso_S0[0], 0), slope=slopes[i], color='black', linestyle='-.', alpha=1)
    else:
        ax.axline(xy1=(tromso_S0[0], 0), slope=slopes[i], color='black', linestyle='-.', alpha=1)

circle1 = plt.Circle((tromso_S0[0], 0),5,color='black')
ax.add_patch(circle1)
ax.text(tromso_S0[0]+10, 0+10,s='EISCAT')

# --- Plot the Simulation Grid ---
# Left Simulation Boundary
L_lat = (data_dict_simulation['grid_lat'][0][0,:] - geoLat[0][0])*111
L_alt = data_dict_simulation['grid_alt'][0][0,:]/stl.m_to_km
ax.scatter(L_lat, L_alt,marker='>', color='tab:cyan',label='Simulation Extent')

# Right Simulation Boundary
R_lat = (data_dict_simulation['grid_lat'][0][-1,:] - geoLat[0][0])*111
R_alt = data_dict_simulation['grid_alt'][0][-1,:]/stl.m_to_km
ax.scatter(R_lat, R_alt,marker='>', color='tab:cyan')

# Bottom Simulation Boundary
B_lat = (data_dict_simulation['grid_lat'][0][:,-1] - geoLat[0][0])*111
B_alt = data_dict_simulation['grid_alt'][0][:,-1]/stl.m_to_km
ax.scatter(B_lat, B_alt,marker='>', color='tab:cyan')

# Top Simulation Boundary
T_lat = (data_dict_simulation['grid_lat'][0][:,0] - geoLat[0][0])*111
T_alt = data_dict_simulation['grid_alt'][0][:,0]/stl.m_to_km
ax.scatter(T_lat, T_alt,marker='>', color='tab:cyan')

ax.legend()
fig.tight_layout()
plt.savefig(r'C:\Users\cfelt\OneDrive - University of Iowa\Research\ACESII\Feltman2025_ACESII_JouleHeating\PLOTS\Plot3\\Plot3_simulation_extent.png')


