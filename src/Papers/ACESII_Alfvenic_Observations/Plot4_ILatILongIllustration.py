# --- Plot4_ILatILongIllustration.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Create the zoomed-in ILatILong plot to aid the illustration in Figure 4

# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import spaceToolsLib as stl
from src.myImports import *
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
print(color.UNDERLINE + f'Plot4_ILatILongIllustration' + color.END)
from matplotlib.ticker import (AutoMinorLocator)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Physics Toggles
ILat_window = [235, 300]
ILongs_window = [-122, -72]
AndoyaCoordinates = [69.2950106, 16.0293464]

# Plot toggles
Figure_width = 10 # in inches
Figure_height =10# in inches

Text_FontSize = 20
Label_FontSize = 25
Tick_FontSize = 20
Tick_Length = 7
Tick_Width = 3
Plot_LineWidth = 5.5
Label_Padding = 15
Tick_Padding = 10
Legend_FontSize = 16
dpi = 200
Scatter_LineWidth_circle = 30
Scatter_LineWidth_X = 50

DistancePlot_LineWidth = 2

# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---
prgMsg('Loading Data')
rocketAttrs, b, c = ACES_mission_dicts()

# LP Particle Data
attitude_files = [glob(r'C:\Data\ACESII\attitude\high\*.cdf')[0], glob(r'C:\Data\ACESII\attitude\low\*.cdf')[0]]
data_dict_attitude_high = loadDictFromFile(inputFilePath=attitude_files[0])
data_dict_attitude_low = loadDictFromFile(inputFilePath=attitude_files[1])
Done(start_time)



# --- Find the target times for the UTC labels ---
target_Time_HF = dt.datetime(2022,11,20,17,24,56)
tTime_index_HF = np.abs(data_dict_attitude_high['Epoch'][0]-target_Time_HF).argmin()


target_Time_LF = dt.datetime(2022,11,20,17,24,54)
tTime_index_LF = np.abs(data_dict_attitude_low['Epoch'][0]-target_Time_LF).argmin()

# --- Get Andoya Coordiantes in kilometers ---
Origin_Lat_km = AndoyaCoordinates[0]*stl.lat_to_meter
Origin_Long_km = stl.long_to_meter(AndoyaCoordinates[1],AndoyaCoordinates[0])

# adjust data based on the origin
ILats_km = [data_dict_attitude_high['ILat_km'][0]- Origin_Lat_km, data_dict_attitude_low['ILat_km'][0]- Origin_Lat_km]
ILongs_km = [data_dict_attitude_high['ILong_km'][0]- Origin_Long_km, data_dict_attitude_low['ILong_km'][0]- Origin_Long_km]

############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

fig, ax = plt.subplots()
fig.set_figwidth(Figure_width)
fig.set_figheight(Figure_height)

#
# # High Flyer - Lat/Long
# rColor='red'
# ax.plot(data_dict_attitude_high['ILong'][0], data_dict_attitude_high['Lat'][0], color=rColor, label='High Flyer Trajectory')
#
#
# # Low Flyer - Lat/Long
# rColor='blue'
# ax.plot(data_dict_attitude_low['Long'][0], data_dict_attitude_low['Lat'][0], color=rColor,label='Low Flyer Trajectory')


# High Flyer - ILat/ILong
rColor='tab:orange'
ax.plot(ILongs_km[0], ILats_km[0], color=rColor, label='HF Magnetic Footpoint (150 km)',linewidth=Plot_LineWidth,linestyle='--',zorder=0)



# Low Flyer - ILat/ILong
rColor='tab:cyan'
ax.plot(ILongs_km[1], ILats_km[1], color=rColor, label='LF Magnetic Footpoint (150 km)',linewidth=Plot_LineWidth,linestyle='--',zorder=0)


# HF Point of Interest
pointOfInterest_HF = np.array([ILongs_km[0][tTime_index_HF], ILats_km[0][tTime_index_HF]])
ax.text(*pointOfInterest_HF + [9, -2], s=r'HF observes $\delta B_{\perp}$'+'\n17:24:56 UTC',weight='bold',fontsize=Text_FontSize,ha='center',va='bottom')
ax.scatter(*pointOfInterest_HF, color='red',marker='o',lw=Scatter_LineWidth_circle,zorder=1) # plot the circle
ax.scatter(*pointOfInterest_HF, color='red',marker='x',lw=Scatter_LineWidth_X,zorder=1) # plot the x


# LF Point of Interest
pointOfInterest_LF = np.array([ILongs_km[1][tTime_index_LF], ILats_km[1][tTime_index_LF]])
ax.text(*pointOfInterest_LF + [-5.5, -3.5], s=r'LF observes $\delta B_{\perp}$'+'\n17:24:54 UTC',weight='bold',fontsize=Text_FontSize,ha='center',va='top')
ax.scatter(*pointOfInterest_LF, color='blue',marker='o',lw=Scatter_LineWidth_circle,zorder=1) # plot the circle
ax.scatter(*pointOfInterest_LF, color='blue',marker='x',lw=Scatter_LineWidth_X,zorder=1) # plot the x


# Draw the black lines connecting the distances
ax.plot(*np.array([pointOfInterest_LF,pointOfInterest_HF]).T,color='black',linewidth=DistancePlot_LineWidth,zorder=0) # hypotenuse
ax.text(*[-96, 262.5],s=f'{ round(np.sqrt( (pointOfInterest_LF[0] - pointOfInterest_HF[0])**2 + (pointOfInterest_LF[1] - pointOfInterest_HF[1])**2),1)} km',fontsize=Text_FontSize-3,weight='bold',rotation=75)

ax.plot(*np.array([ pointOfInterest_LF,  [pointOfInterest_HF[0], pointOfInterest_LF[1]]   ]).T,color='black',linewidth=DistancePlot_LineWidth,zorder=0) # hypotenuse
ax.text(*[-94, 247.5],s=f'{ round(np.sqrt( (pointOfInterest_LF[0] - pointOfInterest_HF[0])**2 ),1)} km',fontsize=Text_FontSize-3,weight='bold')

ax.plot(*np.array([ [pointOfInterest_HF[0],pointOfInterest_LF[1]],  pointOfInterest_HF   ]).T,color='black',linewidth=DistancePlot_LineWidth,zorder=0) # hypotenuse
ax.text(*[-87.5, 263],s=f'{ round(np.sqrt( (pointOfInterest_LF[1] - pointOfInterest_HF[1])**2),1)} km',fontsize=Text_FontSize-3,weight='bold', rotation=90)

# Adjust Plot parameters
ax.set_ylim(*ILat_window)
ax.set_xlim(*ILongs_window)
ax.set_xlabel(f'(-W|+E) [km]', fontsize=Label_FontSize,labelpad=Label_Padding,weight='bold')
ax.set_ylabel(f'(-S|+N) [km]', fontsize=Label_FontSize,labelpad=Label_Padding,weight='bold')
ax.tick_params(axis='both', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax.tick_params(axis='both', which='minor', colors='black', labelsize=Tick_FontSize-5, length=Tick_Length-3, width=Tick_Width-1)
ax.grid(True,which='Major', alpha=0.5)
ax.grid(True,which='minor', alpha=0.2)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.legend(loc='lower right',fontsize=Legend_FontSize)
plt.tight_layout()

plt.savefig(r'C:\Users\cfelt\Desktop\Research\Feltman2024_ACESII_Alfven_Observations\PLOTS\Plot4\Plot4_ILatILongIllustration.png', dpi=dpi)







