# --- Plot1_AllSky.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Stack Plot detailing field-aligned
# particle data along with electric and magnetic signatures


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"
import spaceToolsLib as stl
from myImports import *
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
print(color.UNDERLINE + f'Plot2_Conjugacy' + color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
dpi = 200
Escale = 1000 # what to scale the deltaE field by
Conjugacy_targetILat = [71.31, 72.75]
Conjugacy_EBlimits = 9.5


# --- Cbar ---
# cbarMin, cbarMax = 5E6, 3E9
cbarMin, cbarMax = 1E6, 1E9
cbar_TickLabelSize = 14
my_cmap = stl.apl_rainbow_black0_cmap()
my_cmap.set_bad(color=(0,0,0))


# Plot toggles
Figure_width = 10 # in inches
Figure_height =15# in inches

Text_FontSize = 20
Label_FontSize = 20
Tick_FontSize = 15
Tick_Length = 5
Tick_Width = 2
Plot_LineWidth = 0.8
Label_Padding = 15
Tick_Padding = 10
Legend_FontSize = 13



# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---
prgMsg('Loading Data')

targetILat = Conjugacy_targetILat
targetVar = [targetILat,'ILat']

rocketAttrs, b, c = ACES_mission_dicts()

# delta B
inputMagFiles_high = glob('C:\Data\ACESII\L3\deltaB\high\*Field_Aligned*')[0]
data_dict_mag_high = loadDictFromFile(inputFilePath=inputMagFiles_high, targetVar=targetVar, wKeys_Reduce=['B_e', 'B_r', 'B_p', 'ILat', 'Epoch', 'Alt'])
inputMagFiles_low = glob('C:\Data\ACESII\L3\deltaB\low\*Field_Aligned*')[0]
data_dict_mag_low = loadDictFromFile(inputFilePath=inputMagFiles_low, targetVar=targetVar, wKeys_Reduce=['B_e', 'B_r', 'B_p', 'ILat', 'Epoch', 'Alt'])

# delta E
inputEFIFiles_low = glob('C:\Data\ACESII\L3\deltaE\low\*Field_Aligned*')[0]
data_dict_Efield_low = loadDictFromFile(inputFilePath=inputEFIFiles_low, targetVar=targetVar, wKeys_Reduce=['E_e', 'E_r', 'E_p', 'ILat', 'Epoch', 'Alt'])

data_dict_Efield_low['E_e'][0] = Escale*data_dict_Efield_low['E_e'][0]
data_dict_Efield_low['E_p'][0] = Escale*data_dict_Efield_low['E_p'][0]
data_dict_Efield_low['E_r'][0] = Escale*data_dict_Efield_low['E_r'][0]

# EEPAA Particle Data
inputEEPAA_low = glob('C:\Data\ACESII\L2\low\*eepaa_fullCal*')[0]
data_dict_eepaa_low = loadDictFromFile(inputFilePath=inputEEPAA_low, targetVar=targetVar, wKeys_Reduce=['Differential_Energy_Flux', 'ILat', 'Epoch', 'Alt'])
inputEEPAA_high = glob('C:\Data\ACESII\L2\high\*eepaa_fullCal*')[0]
data_dict_eepaa_high = loadDictFromFile(inputFilePath=inputEEPAA_high, targetVar=targetVar, wKeys_Reduce=['Differential_Energy_Flux', 'ILat', 'Epoch', 'Alt'])


# LP Particle Data
inputLP_low = glob('C:\Data\ACESII\L3\Langmuir\low\*langmuir_fixed*')[1]
data_dict_LP_low = loadDictFromFile(inputFilePath=inputLP_low, targetVar=targetVar, wKeys_Reduce=['ni', 'ILat', 'Epoch'])
inputLP_high = glob('C:\Data\ACESII\L3\Langmuir\high\*langmuir_fixed*')[1]
data_dict_LP_high = loadDictFromFile(inputFilePath=inputLP_high, targetVar=targetVar, wKeys_Reduce=['ni', 'ILat', 'Epoch'])
Done(start_time)


############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################

# --- Calculate Omni-Directional Flux ---
prgMsg('Calculating OmniFlux')

omniDirFlux_low = np.zeros(shape=(len(data_dict_eepaa_low['Differential_Energy_Flux'][0]), len(data_dict_eepaa_low['Energy'][0])))
for tme in range(len(data_dict_eepaa_low['Epoch'][0])):
    for engy in range(len(data_dict_eepaa_low['Energy'][0])):

        sumVal = 0

        for ptch in range(2, 18+1):
            val = data_dict_eepaa_low['Differential_Energy_Flux'][0][tme, ptch, engy]
            if val > 0:
                sumVal += val

        # Average the Omni-flux by the number of bins. ONLY include bins 10deg - 170 since they have full coverage
        omniDirFlux_low[tme][engy] = sumVal/len(range(2, 18+1))

omniDirFlux_high = np.zeros(shape=(len(data_dict_eepaa_high['Differential_Energy_Flux'][0]), len(data_dict_eepaa_high['Energy'][0])))
for tme in range(len(data_dict_eepaa_high['Epoch'][0])):
    for engy in range(len(data_dict_eepaa_high['Energy'][0])):
        sumVal = 0

        for ptch in range(2, 18 + 1):
            val = data_dict_eepaa_high['Differential_Energy_Flux'][0][tme, ptch, engy]
            if val > 0:
                sumVal += val

        # Average the Omni-flux by the number of bins. ONLY include bins 10deg - 170 since they have full coverage
        omniDirFlux_high[tme][engy] = sumVal / len(range(2, 18 + 1))


Done(start_time)

fig, ax = plt.subplots(8, height_ratios=[1.5, 1, 0.5, 0.5, 1.5, 1, 1, 0.5],sharex=False)
fig.set_figwidth(Figure_width)
fig.set_figheight(Figure_height)

# ---HF EEPAA---
axNo = 0
cmap = ax[axNo].pcolormesh(data_dict_eepaa_high['ILat'][0],data_dict_eepaa_high['Energy'][0],omniDirFlux_high.T, cmap=my_cmap,vmin=cbarMin,vmax=cbarMax, norm='log')
ax[axNo].set_ylabel('Energy [eV]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
ax[axNo].set_yscale('log')
ax[axNo].set_ylim(28,1E4)

# --- delta B HF---
axNo +=1
ax[axNo].plot(data_dict_mag_high['ILat'][0],data_dict_mag_high['B_e'][0], color='darkviolet', linewidth=Plot_LineWidth, label='$\delta B_{e}$')
ax[axNo].plot(data_dict_mag_high['ILat'][0],data_dict_mag_high['B_r'][0], color='dodgerblue', linewidth=Plot_LineWidth, label='$\delta B_{r}$')
ax[axNo].set_ylabel('[nT]', fontsize=Label_FontSize,labelpad=Label_Padding)
ax[axNo].tick_params(axis='y',which='both', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].set_ylim(-Conjugacy_EBlimits, Conjugacy_EBlimits)
leg=ax[axNo].legend(loc='upper left',fontsize=Legend_FontSize)
for line in leg.get_lines():
    line.set_linewidth(4)

# --- LP High---
axNo +=1
colorChoice = 'black'
ax[axNo].plot(data_dict_LP_high['ILat'][0], data_dict_LP_high['ni'][0]/1E5, color=colorChoice, linewidth=Plot_LineWidth+1)
ax[axNo].set_ylabel('[10$^{5}$ cm$^{-3}$]', fontsize=Label_FontSize-2, color=colorChoice, labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize-3, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize - 6, length=Tick_Length-2, width=Tick_Width)
ax[axNo].tick_params(axis='x', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length + 4, width=Tick_Width, pad=Tick_Padding)
ax[axNo].tick_params(axis='x', which='minor', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width, pad=Tick_Padding)
ax[axNo].set_ylim(0, 2.55)
ax[axNo].minorticks_on()
ax[axNo].xaxis.set_tick_params(labelbottom=True)
ax[axNo].set_xlabel('ILat [deg] \n time [UTC]', fontsize=Tick_FontSize, weight='bold')
ax[axNo].xaxis.set_label_coords(-0.085, -0.26)



# --- BREAK AXIS ---
axNo +=1
ax[axNo].axis('off')
# ax[axNo].spines[['left', 'right']].set_visible(False)
# ax[axNo].set_yticks(ticks=[],labels=[])
# ax[axNo].set_axisbelow(False)

# --- LF EEPAA---
axNo +=1
cmap = ax[axNo].pcolormesh(data_dict_eepaa_low['ILat'][0], data_dict_eepaa_low['Energy'][0], omniDirFlux_low.T, cmap=my_cmap, vmin=cbarMin, vmax=cbarMax, norm='log')
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize-4, length=Tick_Length-2, width=Tick_Width-1)
ax[axNo].set_ylabel('Energy [eV]', fontsize=Label_FontSize, labelpad=Label_Padding)
ax[axNo].set_yscale('log')
ax[axNo].set_ylim(28,1E4)

# --- delta B LF---
axNo +=1
ax[axNo].plot(data_dict_mag_low['ILat'][0], data_dict_mag_low['B_e'][0], color='darkviolet', linewidth=Plot_LineWidth, label='$\delta B_{e}$')
ax[axNo].plot(data_dict_mag_low['ILat'][0], data_dict_mag_low['B_r'][0], color='dodgerblue', linewidth=Plot_LineWidth, label='$\delta B_{r}$')
ax[axNo].set_ylabel('[nT]', fontsize=Label_FontSize, color='black', labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=0, length=0, width=0)
ax[axNo].set_ylim(-Conjugacy_EBlimits, Conjugacy_EBlimits)
leg=ax[axNo].legend(loc='upper left',fontsize=Legend_FontSize)
for line in leg.get_lines():
    line.set_linewidth(4)

# --- delta E LF---
axNo +=1
ax[axNo].plot(data_dict_Efield_low['ILat'][0], data_dict_Efield_low['E_r'][0], linewidth=Plot_LineWidth, color='red',label=r'$\delta E_{r}$')
ax[axNo].plot(data_dict_Efield_low['ILat'][0], data_dict_Efield_low['E_e'][0], linewidth=Plot_LineWidth, color='blue',label=r'$\delta E_{e}$')
ax[axNo].set_ylabel('[mV/m]', fontsize=Label_FontSize, color='black', labelpad=Label_Padding)
ax[axNo].tick_params(axis='y', which='both', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=0, length=0, width=0)
ax[axNo].set_ylim(-Conjugacy_EBlimits, Conjugacy_EBlimits)
leg = ax[axNo].legend(loc='upper left',fontsize=Legend_FontSize)
# get the individual lines inside legend and set line width
for line in leg.get_lines():
    line.set_linewidth(4)

# --- LP Low ---
axNo +=1
ax[axNo].plot(data_dict_LP_low['ILat'][0], data_dict_LP_low['ni'][0]/1E5, color='black', linewidth=Plot_LineWidth+1)
ax[axNo].set_ylabel('[10$^{5}$ cm$^{-3}$]', fontsize=Label_FontSize-2, color='black', labelpad=Label_Padding )
ax[axNo].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize-3, length=Tick_Length, width=Tick_Width)
ax[axNo].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize - 6, length=Tick_Length-2, width=Tick_Width)
ax[axNo].tick_params(axis='x', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length + 4, width=Tick_Width, pad=Tick_Padding)
ax[axNo].tick_params(axis='x', which='minor', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width, pad=Tick_Padding)
ax[axNo].set_ylim(0, 2.55)
ax[axNo].set_xlabel('ILat [deg] \n time [UTC]', fontsize=Tick_FontSize, weight='bold')
ax[axNo].xaxis.set_label_coords(-0.085, -0.26)
ax[axNo].minorticks_on()

# --- get UTC labels and ILat Labels together ---
xTickLabels = ax[axNo].axes.get_xticklabels()
xTick_ILatLocations = [float(tickVal.get_text()) for tickVal in xTickLabels]
xTick_newLabels_high = [f'{iLat}\n{data_dict_eepaa_high["Epoch"][0][np.abs(data_dict_eepaa_high["ILat"][0] - iLat).argmin()].strftime("%H:%M:%S")}' for iLat in xTick_ILatLocations]
xTick_newLabels_low = [f'{iLat}\n{data_dict_eepaa_low["Epoch"][0][np.abs(data_dict_eepaa_low["ILat"][0] - iLat).argmin()].strftime("%H:%M:%S")}' for iLat in xTick_ILatLocations]
ax[2].set_xticks(xTick_ILatLocations,labels=xTick_newLabels_high)
ax[7].set_xticks(xTick_ILatLocations,labels=xTick_newLabels_low)

# -- Do some minor adjustments to labels/margins/limits ---
for i in range(8):
    ax[i].margins(0)
    ax[i].set_xlim(targetILat[0], targetILat[1])
fig.align_ylabels(ax[:])

# --- cbar 1---
cax = fig.add_axes([0.91, 0.796, 0.03, 0.184])
cbar = plt.colorbar(cmap, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbar_TickLabelSize + 5)

# --- cbar 2---
cax = fig.add_axes([0.91, 0.367, 0.03, 0.184])
cbar = plt.colorbar(cmap, cax=cax)
cbar.ax.minorticks_on()
cbar.ax.tick_params(labelsize=cbar_TickLabelSize + 5)

fig.subplots_adjust(left=0.12, bottom=0.06, right=0.9, top=0.98,hspace=0)  # remove the space between plots
plt.savefig(r'C:\Users\cfelt\Desktop\rockets\ACES-II\Papers\ACESII_Alfven_Observations\PLOTS\Plot2\Plot2_ConjugacyStack.png', dpi=dpi)







