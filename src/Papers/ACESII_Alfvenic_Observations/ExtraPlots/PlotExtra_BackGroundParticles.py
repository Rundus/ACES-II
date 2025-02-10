# --- Plots4_pitchAnglePlots.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Recreation of DiffEFlux pitch angle Plots, focusing on
# a few particle signatures


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import itertools

from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
from my_matplotlib_Assets.colorbars.matlab_parula import matlab_parula_cmap
from ACESII_code.class_var_func import Re

print(color.UNDERLINE + f'Plot4_BackGroundParticles' + color.END)


#################
# --- TOGGLES ---
#################

### DATA TOGGLES ###
targetTimes = [dt.datetime(2022, 11, 20, 17, 24, 59, 500000), dt.datetime(2022, 11, 20, 17, 25, 2, 000000)]
#             100, 110,120,130,140,150,160,170
wPitchBins = [10, 11, 12, 13, 14, 15, 16, 17, 18] # Index of the pitch bins to use in the figure
# wPitchBins = [13] # Index of the pitch bins to use in the figure
Engylimit = 17 # corresponds t the highest energy used
AltMin, AltMax = 400, 6000
N = 51



### PLOT TOGGLES ###

figure_height =8
figure_width = 12
# cmap = matlab_parula_cmap()
cmap = apl_rainbow_black0_cmap()
cbarLow_counts, cbarHigh_counts = 1, 100
titleFontSize = 18
labelsFontSize = 20
tickFontSize = 15
tickWidth = 1
tickLength = 3
lineWidth = 2
dpi = 100


# --- --- --- --- --- ---
# --- LOAD IN THE DATA ---
# --- --- --- --- --- ---

prgMsg('Loading Data')

dataKey = 'eepaa'
labelNam = 'Counts'

# Attitude data (for geomagnetic lat/long info)
inputFiles_Attitude = [glob(r'C:\Data\ACESII\attitude\high\*Attitude_Solution*')[0], glob(r'C:\Data\ACESII\attitude\low\*Attitude_Solution*')[0]]
data_dict_attitude_high = loadDictFromFile(inputFilePath=inputFiles_Attitude[0],targetVar=[targetTimes, 'Epoch'],wKeys_Reduce=['Alt','Epoch'])

# EEPAA Particle Data
cbarLow, cbarHigh = cbarLow_counts, cbarHigh_counts
inputFiles_eepaa_counts = [glob('C:\Data\ACESII\L1\high\*eepaa_fullCal*')[0], glob('C:\Data\ACESII\L1\low\*eepaa_fullCal*')[0]]
data_dict_counts_high = loadDictFromFile(inputFilePath=inputFiles_eepaa_counts[0],targetVar=[targetTimes, 'Epoch'],wKeys_Reduce=['eepaa','Epoch'])
Done(start_time)


# --- --- --- --- --- --- --- --- --
# --- Calc Initial Vpara & Vperp ---
# --- --- --- --- --- --- --- --- --
prgMsg('Calculating Vperp and Vparallel')
Energy = data_dict_counts_high['Energy'][0]
Pitch = data_dict_counts_high['Pitch_Angle'][0]
pitchI = np.zeros(shape=len(data_dict_counts_high['Epoch'][0])*len(wPitchBins)*len(Energy))
EmagI = np.zeros(shape=len(data_dict_counts_high['Epoch'][0])*len(wPitchBins)*len(Energy))
countsI = np.zeros(shape=len(data_dict_counts_high['Epoch'][0])*len(wPitchBins)*len(Energy))
AltI = np.zeros(shape=len(data_dict_counts_high['Epoch'][0])*len(wPitchBins)*len(Energy))
ranges = [range(len(data_dict_counts_high['Epoch'][0])), wPitchBins, range(Engylimit,len(Energy))]

counter = 0
for tme, ptch, engy in itertools.product(*ranges):
    pitchI[counter] = Pitch[ptch]
    EmagI[counter] = Energy[engy]
    countsI[counter] = data_dict_counts_high['eepaa'][0][tme][ptch][engy]
    AltI[counter] = int((1/1000)*sum(data_dict_attitude_high['Alt'][0])/len(data_dict_attitude_high['Alt'][0]))
    counter+=1

Done(start_time)


# --- --- --- --- --- --- --- --- --- --- --- --
# --- Calc Future Vparallel at all altitudes ---
# --- --- --- --- --- --- --- --- --- --- --- --
prgMsg('Calculating at all Altitudes')
simAlt = np.linspace(AltMin, AltMax, N)
simEpara = np.zeros(shape=(len(simAlt), len(countsI)))
simCounts = np.zeros(shape=(len(simAlt), len(countsI)))

def calcEpara(EmagI, alphaI, z):
    pitch = np.radians(alphaI)
    Val = EmagI * (np.cos(pitch)*np.cos(pitch)*(1 - ( (Re + AltI[0])/(Re + z) )**(3/2)  )  +  np.sin(pitch)*np.sin(pitch))
    return Val

# Calculate Vperp and Vpara for all time
for indx, alt in enumerate(simAlt):
    simEpara[indx] = calcEpara(EmagI, pitchI, np.array([alt for i in range(len(AltI))]))
    simCounts[indx] = countsI


Done(start_time)




# --- --- --- --- --- --- --- --- --- --
# --- Reformat the data for plotting ---
# --- --- --- --- --- --- --- --- --- --
prgMsg('Calculating Data Histogram')

# Epara_bins = sorted(list(set(list(itertools.chain(*[list(set(ar)) for ar in simEpara])))))
Epara = Energy[:Engylimit:-1]
Epara_bins = np.histogram_bin_edges(Epara, bins=len(Epara))
# Epara_bins = [-5 + 10*i for i in range(86)]
# Epara = [10*i for i in range(85)]
dataHistogram = np.zeros(shape=(len(simAlt), len(Epara)))

for indx, alt in tqdm(enumerate(simAlt)):

    # expand the Vpara variable for every count detected at a specific velocity
    # EparaExpanded = list(itertools.chain(*[[simEpara[indx][idx] for i in range(int(val))] for idx, val in enumerate(simCounts[indx])]))

    EparaExpanded = []

    # copy the energy for every count value
    for i, val in enumerate(simCounts[indx]):
        Engy = simEpara[indx][i]
        EparaExpanded.append([Engy for k in range(int(val))])

    # collapse the array
    EparaExpanded = list(itertools.chain(*EparaExpanded))

    # Make a histogram of the VparaExpanded
    hist, bins = np.histogram(EparaExpanded, Epara_bins)

    dataHistogram[indx] = 2*hist

Done(start_time)


# --- --- --- --- --- --- --- --- --- --
# --- Calculate STEB5 poly fit curve ---
# --- --- --- --- --- --- --- --- --- --
a1_STEB5 = 4.852E3
a2_STEB5 = -9.860E6
Emax_STEB5 = 828.27
Emin_STEB5 = 28.22
def polyfit(Emin, a1, a2, Emax):
    alpha = 10 * np.pi/180
    return a2*1000*np.sqrt(0.5*m_e/q0) * (Emin**(-0.5) - Emax**(-0.5)) / (np.cos(alpha)) + a1


polyFit_Altitudes = polyfit(np.array(Epara),
                            np.array([a1_STEB5 for i in range(len(Epara))]),
                            np.array([a2_STEB5 for i in range(len(Epara))]),
                            np.array([Emax_STEB5 for i in range(len(Epara))]))


############################
# --- --- --- --- --- --- --
# --- START THE PLOTTING ---
# --- --- --- --- --- --- --
############################
prgMsg('Beginning Plot')
fig, ax = plt.subplots()
fig.set_size_inches(figure_width, figure_height)


# the histogram plot
# cmap = ax.pcolormesh(Epara, simAlt, dataHistogram, cmap=cmap, norm='log', vmin=int(1E6),vmax=int(1E7))
cmap = ax.pcolormesh(Epara, simAlt, dataHistogram, cmap=cmap, norm='log')
ax.axvline(Emin_STEB5,alpha=1, linewidth=3)
ax.axvline(Emax_STEB5,alpha=1, linewidth=3)
ax.set_ylim(400, AltMax)
ax.set_ylabel('Altitude [km]', fontsize=labelsFontSize)
ax.set_xlabel('$E_{\parallel}$ [eV]', fontsize=labelsFontSize)
ax.set_xscale('log')
ax.set_xlim(1, 1000)
ax.set_xlim(min(Epara_bins), max(Epara_bins))
ax.tick_params(axis='both', which='major', labelsize=tickFontSize, width=tickWidth, length=tickLength)
ax.tick_params(axis='both', which='minor', labelsize=tickFontSize-5, width=tickWidth, length=tickLength/2)

# The polyfit plot
ax.plot(Epara, polyFit_Altitudes,color='black',label='STEB5 polyfit', linestyle='--',alpha=0.75, linewidth=3)
ax.legend()

# colorbar
cbar = plt.colorbar(cmap)
cbar.set_label('Counts (x2)', fontsize=labelsFontSize)

# output the figure
plt.tight_layout()
fileOutName = rf'C:\Users\cfelt\OneDrive\Desktop\Papers\ACESII_Alfven_Observations\Plot4\\Plot4_backgroundParticles.png'
plt.savefig(fileOutName, dpi=dpi)
Done(start_time)