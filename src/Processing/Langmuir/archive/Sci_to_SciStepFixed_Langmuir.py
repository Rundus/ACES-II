import spaceToolsLib as stl
from scipy import signal
from matplotlib import pyplot as plt
from src.my_imports import *
inputFile_Langmuir = 'C:\Data\ACESII\L3\Langmuir\low\ACESII_36364_langmuir_fixed.cdf'
data_dict_langmuir = deepcopy(stl.loadDictFromFile(inputFile_Langmuir,wKeys_Load=['ni','Epoch','ILat']))
d = data_dict_langmuir['ni'][0][1000:4000]

# Convolution part
dary = np.array(d)
dary -= np.average(dary)

step = np.hstack((np.ones(len(dary)), -1*np.ones(len(dary))))

dary_step = np.convolve(dary, step, mode='valid')

# Get the peaks of the convolution
peaks = signal.find_peaks(dary_step, width=20)[0]

# plots
plt.figure()
plt.plot(dary)
plt.plot(dary_step/10)

for ii in range(len(peaks)):
    plt.plot((peaks[ii], peaks[ii]), (-1500, 1500), 'r')

plt.show()