from scipy.signal import butter, lfilter

def butter_bandpass(lowcut, highcut, fs, order=5):
    return butter(order, [lowcut, highcut], fs=fs, btype='Bandpass')
    # return butter(order, lowcut, fs=fs, btype='Highpass')

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

from ACESII_code.myImports import *

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import freqz

# Sample rate and desired cutoff frequencies (in Hz).
fs = 128
lowcut = 1
highcut = 15
orderVal = 4

# Plot the frequency response for a few different orders.
plt.figure(1)
plt.clf()
for order in [orderVal]:
    b, a = butter_bandpass(lowcut=lowcut,highcut=highcut, fs=fs, order=order)
    w, h = freqz(b, a, fs=fs, worN=20000)
    plt.plot(w, abs(h), label="order = %d" % order)

plt.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)],
         '--', label='sqrt(0.5)')
plt.xlim(-0.1, 13)
plt.title('Butterworth Filter RollOff\n'
          f'Sample Freq: {fs}\n'
          f'Order: {orderVal}\n'
          r'$f_{L}$: '+f'{lowcut}\n'
          r'$f_{H}$: '+f'{highcut}\n')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain')
plt.grid(True)
plt.legend(loc='best')

plt.show()