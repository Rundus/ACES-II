# --- Plots9_PoyntingFlux.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

from src.my_imports import *

start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import spaceToolsLib as stl
import matplotlib.pyplot as plt
print(stl.color.UNDERLINE + f'Plot9_keyObservations' + stl.color.END)

#################
# --- TOGGLES ---
#################
figure_height = 2.3*(3.75)
figure_width = (15)

plot_LineWidth = 3
plot_textFontSize = 20
plot_MarkerSize = 20
plot_Colors = ['tab:blue', 'tab:red', 'tab:orange', 'tab:green', 'tab:purple', 'tab:olive', 'black']
text_FontSize = 25


title_FontSize = 25

labels_FontSize = 22
labels_subplot_fontsize = 19

Legend_FontSize = 22.5
legend_SubAxes_FontSize = 15

tick_LabelSize = 17
tick_SubplotLabelSize = 15
tick_Width = 2
tick_Length = 4

cbar_FontSize = 15
dpi = 100

dispersiveRegionTargetTime = [dt.datetime(2022,11,20,17,24,55,500000),
                              dt.datetime(2022,11,20,17,25,4,250000)]



# --- plot Poynting Flux Toggles ---
wSTEBtoPlot = [1, 2, 3, 4, 5] # STEB number NOT index. Don't -1
PoyntingScale = 1E3# convert from W/m^2 to ergs/cm^2



##########################
# --- --- --- --- --- ---
# --- LOADING THE DATA ---
# --- --- --- --- --- ---
##########################
stl.prgMsg('Loading Data')

#========================
# --- High Flyer Data ---
#========================

targetVarName = 'Epoch'
targetVar = dispersiveRegionTargetTime

# Magnetometer Data - High Flyer
inputFile_B_high = 'C:\Data\ACESII\L2\high\ACESII_36359_l2_RingCore_ENU.cdf'  # get the B data
data_dict_B_high = stl.loadDictFromFile(inputFile_B_high,targetVar=[targetVar, targetVarName])
inputFile_deltaB_high = glob('C:\Data\ACESII\L3\deltaB\high\ACESII_36359_RingCore_Field_Aligned_WL250_stitchedFlight.cdf')[0] # get the deltaB data
data_dict_deltaB_high = deepcopy(stl.loadDictFromFile(inputFile_deltaB_high,targetVar=[targetVar, targetVarName]))

# Langmuir Data - High Flyer
inputFile_Langmuir_high = 'C:\Data\ACESII\L3\Langmuir\high\ACESII_36359_l3_langmuir_fixed.cdf'
data_dict_langmuir_high = deepcopy(stl.loadDictFromFile(inputFile_Langmuir_high,targetVar=[targetVar, targetVarName], wKeys_Load=['ni', 'Epoch', 'ILat']))
indexVals = [np.abs(data_dict_langmuir_high['Epoch'][0] - tme).argmin() for i, tme in enumerate(data_dict_B_high['Epoch'][0])]
data_dict_langmuir_high['ni'][0] = deepcopy(data_dict_langmuir_high['ni'][0][indexVals])
data_dict_langmuir_high['Epoch'][0] = deepcopy(data_dict_langmuir_high['Epoch'][0][indexVals])

# EISCAT Data - High Flyer (# Up-sample the EISCAT data (it's fine since the EISCAT variables are very slowly varying))
inputFile_EISCAT_high = 'C:\Data\ACESII\science\EISCAT_ACESII_Slice\high\ACESII_36359_EISCAT_Tromso_rktSlice.cdf'
data_dict_EISCAT_high = deepcopy(stl.loadDictFromFile(inputFile_EISCAT_high, targetVar=[targetVar, targetVarName],wKeys_Load=['Ion_Comp', 'Op_Comp','ILat','Epoch']))
data_dict_EISCAT_interp_high = stl.InterpolateDataDict(InputDataDict=data_dict_EISCAT_high,InputEpochArray=data_dict_EISCAT_high['Epoch'][0],targetEpochArray=data_dict_B_high['Epoch'][0],wKeys=[])
data_dict_EISCAT_high = deepcopy(data_dict_EISCAT_interp_high)

# Electron Energy Flux - High Flyer
inputFile_EFlux_high = 'C:\Data\ACESII\L3\Energy_Flux\high\ACESII_36359_eepaa_Energy_Flux.cdf'
data_dict_EFlux_high = deepcopy(stl.loadDictFromFile(inputFile_EFlux_high, targetVar=[targetVar, targetVarName]))

#=======================
# --- Low Flyer Data ---
#=======================

# Magnetometer Data - Low Flyer
inputFile_B_low = 'C:\Data\ACESII\L2\low\ACESII_36364_l2_RingCore_ENU.cdf'  # get the B data
data_dict_B_low = stl.loadDictFromFile(inputFile_B_low,targetVar=[targetVar, targetVarName])
inputFile_deltaB_low = glob('C:\Data\ACESII\L3\deltaB\low\ACESII_36364_RingCore_Field_Aligned_WL601.cdf')[0] # get the deltaB data
data_dict_deltaB_low = deepcopy(stl.loadDictFromFile(inputFile_deltaB_low,targetVar=[targetVar, targetVarName]))

# Langmuir Data - Low Flyer
inputFile_Langmuir_low = 'C:\Data\ACESII\L3\Langmuir\low\ACESII_36364_langmuir_fixed.cdf'
data_dict_langmuir_low = deepcopy(stl.loadDictFromFile(inputFile_Langmuir_low,targetVar=[targetVar, targetVarName], wKeys_Load=['ni', 'Epoch', 'ILat']))
indexVals = [np.abs(data_dict_langmuir_low['Epoch'][0] - tme).argmin() for i, tme in enumerate(data_dict_B_low['Epoch'][0])]
data_dict_langmuir_low['ni'][0] = deepcopy(data_dict_langmuir_low['ni'][0][indexVals])
data_dict_langmuir_low['Epoch'][0] = deepcopy(data_dict_langmuir_low['Epoch'][0][indexVals])

# EISCAT Data - Low Flyer (# Up-sample the EISCAT data (it's fine since the EISCAT variables are very slowly varying))
inputFile_EISCAT_low = 'C:\Data\ACESII\science\EISCAT_ACESII_Slice\low\ACESII_36364_EISCAT_Tromso_rktSlice.cdf'
data_dict_EISCAT_low = deepcopy(stl.loadDictFromFile(inputFile_EISCAT_low, targetVar=[targetVar, targetVarName],wKeys_Load=['Ion_Comp', 'Op_Comp','ILat','Epoch']))
data_dict_EISCAT_interp_low = stl.InterpolateDataDict(InputDataDict=data_dict_EISCAT_low,InputEpochArray=data_dict_EISCAT_low['Epoch'][0],targetEpochArray=data_dict_B_low['Epoch'][0],wKeys=[])
data_dict_EISCAT_low = deepcopy(data_dict_EISCAT_interp_low)

# Poynting Flux from direct ExB
inputFile_PoyntingFlux_Low = 'C:/Data/ACESII/science/PoyntingFlux/low/ACESII_36364_PoyntingFlux_Field_Aligned.cdf'
data_dict_Poynting_low = deepcopy(stl.loadDictFromFile(inputFile_PoyntingFlux_Low, targetVar=[targetVar, targetVarName],wKeys_Load=['S_p','Epoch']))

stl.Done(start_time)

##########################
# --- --- --- --- --- ---
# --- PREPARE THE DATA ---
# --- --- --- --- --- ---
##########################

#===================
# --- HIGH FLYER ---
#===================

# Get the times for the High Flyer STEBs
from src.Science.AlfvenSingatureAnalysis.Particles.dispersionAttributes import dispersionAttributes
STEBtimes = [dispersionAttributes.keyDispersionDeltaT[idx-1] for idx in wSTEBtoPlot]
# --- calculate Flux ---
B0 = 1E-9 * data_dict_B_high['Bmag'][0]
dB_e = data_dict_deltaB_high['B_e'][0] # in nanotesla
dB_r = data_dict_deltaB_high['B_r'][0]
ni = (stl.cm_to_m ** 3) * data_dict_langmuir_high['ni'][0]
rhoCalc_HF = (ni*data_dict_EISCAT_high['Ion_Comp'][0]*((stl.IonMasses[4] + stl.IonMasses[5] + stl.IonMasses[7])/3) + stl.IonMasses[1]*ni*data_dict_EISCAT_high['Op_Comp'][0])
VA_t_HF = B0 / np.sqrt(stl.u0 * rhoCalc_HF)
VA_avg_HF = sum(VA_t_HF)/len(VA_t_HF)
dBperp_HF = np.array([np.sqrt((dB_e[i]*1E-9)**2 + (dB_r[i]*1E-9)**2) for i in range(len(dB_e))])
S_est_HF = VA_t_HF*(dBperp_HF**2)/ (stl.u0) # calculate Estimated Poynting Flux


#==================
# --- LOW FLYER ---
#==================

# Get the times for the High Flyer STEBs
LaunchDateTime = pycdf.lib.datetime_to_tt2000(dt.datetime(2022,11,20,17,20,00,000000))
rktTime_deltaB_LF = stl.EpochTo_T0_Rocket(InputEpoch=data_dict_deltaB_low['Epoch'][0], T0=LaunchDateTime)

# --- calculate Flux ---
B0 = 1E-9 * data_dict_B_low['Bmag'][0]
dB_e = data_dict_deltaB_low['B_e'][0] # in nanotesla
dB_r = data_dict_deltaB_low['B_r'][0]
ni = (stl.cm_to_m ** 3) * data_dict_langmuir_low['ni'][0]
rhoCalc_LF = (ni*data_dict_EISCAT_low['Ion_Comp'][0]*((stl.IonMasses[4] + stl.IonMasses[5] + stl.IonMasses[7])/3) +stl.IonMasses[1]*ni*data_dict_EISCAT_low['Op_Comp'][0])
VA_t_LF = B0 / np.sqrt(stl.u0 * rhoCalc_LF)
VA_avg_LF = sum(VA_t_LF)/len(VA_t_LF)
dBperp = np.array([np.sqrt((dB_e[i]*1E-9)**2 + (dB_r[i]*1E-9)**2) for i in range(len(dB_e))])
S_est_LF = VA_t_LF*(dBperp**2)/ (stl.u0) # calculate Estimated Poynting Flux


############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

stl.prgMsg('Beginning Plot')
fig, ax = plt.subplots(nrows=2,sharex=True)
fig.set_size_inches(figure_width, figure_height)

fig.suptitle('Dispersive Region Energy Flux', weight='bold',fontsize=title_FontSize)

# --- High Flyer - Electron Flux ---
ax[0].plot(data_dict_EFlux_high['Epoch'][0],data_dict_EFlux_high['Energy_Flux_Downward'][0], linewidth=plot_LineWidth,zorder=2,label=r'Parallel $e^{-}$ Flux',color='maroon')
for k, tme in enumerate(STEBtimes):
    lowIdx = np.abs(data_dict_EFlux_high['Epoch'][0] - tme[0]).argmin()
    highIdx = np.abs(data_dict_EFlux_high['Epoch'][0] - tme[1]).argmin()
    ax[0].axvspan(data_dict_EFlux_high['Epoch'][0][lowIdx], data_dict_EFlux_high['Epoch'][0][highIdx], color=plot_Colors[k], alpha=0.25,zorder=1)
    textPlacement = pycdf.lib.tt2000_to_datetime(
        int((pycdf.lib.datetime_to_tt2000(data_dict_EFlux_high['Epoch'][0][lowIdx])+ pycdf.lib.datetime_to_tt2000(data_dict_EFlux_high['Epoch'][0][highIdx]))/2)
    )
    props = dict(boxstyle='round', facecolor='white', alpha=1, lw=4)
    ax[0].text(textPlacement, 1, f'S{wSTEBtoPlot[k]}', ha='center', weight='bold',fontsize=plot_textFontSize,bbox=props,va='center')

# --- High Flyer - Poynting Flux ---
ax[0].plot(data_dict_deltaB_high['Epoch'][0], PoyntingScale*S_est_HF, linewidth=plot_LineWidth,zorder=2,color='tab:orange', label=r'$\delta S_{\parallel} \sim \frac{\delta B_{\perp}^{2}}{\mu_{0}} V_{A}$')
ax[0].set_ylabel('High Flyer\n'+'[erg/cm$^{2}$s]',fontsize=labels_FontSize, weight='bold')
ax[0].set_xmargin(0)
ax[0].tick_params(axis='both', which='major', labelsize=tick_LabelSize, width=tick_Width, length=tick_Length)
ax[0].tick_params(axis='both', which='minor', labelsize=tick_LabelSize-2, width=tick_Width, length=tick_Length)
ax[0].tick_params(axis='x', which='major', labelbottom=False)
ax[0].text(dt.datetime(2022,11,20,17,24,56,000000),0.4,'(a)',color='black',fontsize=text_FontSize+10,va='center',ha='center',weight='bold')
ax[0].set_ylim(0.001, 1)
ax[0].legend(loc='upper right',fontsize=Legend_FontSize)
ax[0].set_yscale('log')
ax[0].grid(alpha=0.5)

# --- Low Flyer - Poynting Flux ---
ax[1].plot(data_dict_deltaB_low['Epoch'][0], PoyntingScale*data_dict_Poynting_low['S_p'][0], linewidth=plot_LineWidth,zorder=2, color='tab:blue',label=r'$(\mathbf{\vec{E}} \times \mathbf{\vec{B}})_{\parallel}$')
ax[1].plot(data_dict_deltaB_low['Epoch'][0], PoyntingScale*S_est_LF, linewidth=plot_LineWidth,zorder=2, color='tab:orange', label=r'$\delta S_{\parallel} \sim \frac{\delta B_{\perp}^{2}}{\mu_{0}} V_{A}$')
ax[1].set_ylabel('Low Flyer\n'+'[erg/cm$^{2}$s]',fontsize=labels_FontSize, weight='bold')
ax[1].set_ylim(-0.3E-2, 5.4E-2)
ax[1].set_xmargin(0)
ax[1].tick_params(axis='both', which='major', labelsize=tick_LabelSize, width=tick_Width, length=tick_Length)
ax[1].tick_params(axis='both', which='minor', labelsize=tick_LabelSize-2, width=tick_Width, length=tick_Length)
ax[1].tick_params(axis='x', which='major')
ax[1].grid(alpha=0.5)
ax[1].set_xlabel('Time [UTC]',fontsize=labels_FontSize, weight = 'bold')
ax[1].legend(fontsize=Legend_FontSize,loc='upper right')
ax[1].text(dt.datetime(2022,11,20,17,24,56,000000),0.047,'(b)',color='black',fontsize=text_FontSize+10,va='center',ha='center',weight='bold')


plt.tight_layout()
plt.savefig(r'C:\Users\cfelt\Desktop\Research\ACESII\Feltman2025_ACESII_Alfven_Observations\PLOTS\Plot9\Plot9_PoyntingFlux_base.png',dpi=dpi)
stl.Done(start_time)