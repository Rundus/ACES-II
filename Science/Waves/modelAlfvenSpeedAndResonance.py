# --- ModelAlfvenSpeed.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Plot a model of the Alfven Speed vs altitude in order to analyze
# Alfven dispersion features. Also, plot vs altitude a resonance condition
# useful for determining what particles will resonate with alfven waves.


# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *
from ACESII_code.class_var_func import Re


start_time = time.time()
# --- --- --- --- ---

# --- --- --- --- --- --- --- ---
# --- ALFVEN VELOCITY TOGGLES ---
# --- --- --- --- --- --- --- ---
SECTION_AlfvenVelocityPlot = True

# assume the kinetic term is just 2^-0.5
simplifyKinetic = True

# Altitude Range (in KILOMETERS) to plot alfven speed over
wkperp = (2*np.pi)/100 # This is ANGULAR wave number: k = 2 pi/lambda
AltLow = 1 # in km
AltHigh = 21*Re # in km
xscaling = Re # how to scale the xaxis of the plot
yscaling= 1E7 # how to scale the yaxis of the plot. Nominally: (1E4)*(1E3)
plotylabel = '$V_{A}$ [10,000 km/s]' # nominally: [10,000 km/s]
plotxlabel = 'Altitude [$R_{E}$]'
ylimits = (0,4)
xlimits = (21,-0.1) # the x-axis is inverted remember


# --- --- --- --- --- --- --- ---
# --- ALFVEN RESONANCE TOGGLES ---
# --- --- --- --- --- --- --- ---
SECTION_AlfvenResonancePlot = False
energyScaling = 1000 # converst eV to keV

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import pyIGRF
from ACESII_code.class_var_func import lightSpeed,u0, m_e, q0



# profile taken from Kletzing: Alfven Wave Electron Acceleration
def density(z): # returns density for altitude "z [km]" in m^-3
    h = 0.06*Re # in km from E's surface
    n0 = 6E4
    n1 = 1.34E7
    z0 = 0.05*Re # in km from E's surface
    n = n0*np.exp(-1*(z-z0)/h) + n1*(z**(-1.55)) # calculated density (in cm^-3)
    return (cm_to_m**3)*n

def kineticTerm(kperp, z, simplify): # represents the denominator of the Alfven velocity term: 1/(1 + (kperp*c/omega_pe)^2)^1/2
    if simplify:
        y = 1/np.sqrt(2)
    else:
        plasmaFreq = np.sqrt((density(z) * q0 * q0) / (m_e * ep0))
        y = 1 / np.sqrt(1 + (kperp * lightSpeed / plasmaFreq) ** 2)
    return y

def AlfvenSpeed(z,lat,long,year,kperp,simplify):
    # -- Output order forpyIGRF.igrf_value ---
    # [0] Declination (+ E | - W)
    # [1] Inclination (+ D | - U)
    # [2] Horizontal Intensity
    # [3] North Comp (+ N | - S)
    # [4] East Comp (+ E | - W)
    # [5] Vertical Comp (+ D | - U)
    # [6] Total Field

    B = pyIGRF.igrf_value(lat, long, z, year)
    V_A = (B[6]*1E-9)/np.sqrt(u0 * density(z) * IonMasses[0])

    if simplify:
        V = V_A*kineticTerm(1, z, simplify)
    else:
        V = V_A*kineticTerm(kperp, z, simplify)
    return V




def main(AltLow, AltHigh):

    # altitude range
    N = 10000
    altitudeAxis = np.linspace(AltLow, AltHigh, N)

    # fig,ax = plt.subplots()
    # ax.plot(altitudeAxis,density(altitudeAxis))
    # plt.show()

    # determine the alfven speed
    year = 2022 + 323 / 365  # Corresponds to 11/20/2022
    LaunchLat = 69.294167
    LaunchLong = 16.020833
    alfvenAxis = np.sqrt(2)*np.array([AlfvenSpeed(z, LaunchLat, LaunchLong, year, wkperp, simplifyKinetic) for z in altitudeAxis])

    # Resonance Bands
    potentials = [10, 100, 500]
    acceleration = []
    centerEnergy = []
    deceleration = []

    for i, potential in enumerate(potentials):

        atemp, ctemp, dtemp = [], [], []

        for j in range(len(alfvenAxis)):
            atemp.append(0.5 * (m_e / q0) * ((alfvenAxis[j] - np.sqrt(2 * q0 * potential / m_e))**2))
            dtemp.append(0.5 * (m_e / q0) * ((alfvenAxis[j] + np.sqrt(2 * q0 * potential / m_e))**2))
            ctemp.append(0.5*(m_e/q0) * (alfvenAxis[j]**2))

        acceleration.append(atemp)
        deceleration.append(dtemp)
        centerEnergy.append(ctemp)


    # normalize the data to whatever you wish
    xData = altitudeAxis/(xscaling) # to 1000 km
    yData = alfvenAxis/(yscaling) # to 10,000 km/s

    yData = (1/(u0*yscaling*alfvenAxis))

    acceleration = np.array(acceleration)/energyScaling
    deceleration = np.array(deceleration) / energyScaling
    centerEnergy = np.array(centerEnergy) / energyScaling

    if SECTION_AlfvenVelocityPlot:
        # --- --- --- ---
        # --- PLOTTING ---
        # --- --- --- ---
        fig, ax = plt.subplots()
        ax.set_ylabel(plotylabel)
        ax.set_xlabel(plotxlabel)
        fig.suptitle('Alfven Speed ($V_{A}$) vs Altitude')
        ax.plot(xData, yData)
        ax.invert_xaxis()
        # ax.set_ylim(*ylimits)
        ax.set_xlim(*xlimits)
        ax.set_yscale('log')
        plt.show()
        
    if SECTION_AlfvenResonancePlot:


        # --- --- --- ---
        # --- PLOTTING ---
        # --- --- --- ---
        fig, ax = plt.subplots(len(potentials) + 1, sharex=True)
        plt.subplots_adjust(wspace=0, hspace=0)

        # --- alfven Speed Plot ---
        ax[0].set_ylabel(plotylabel)
        fig.suptitle('Alfven Speed ($V_{A}$) vs Altitude')
        ax[0].plot(xData, yData)
        ax[0].invert_xaxis()
        # ax[0].set_ylim(*ylimits)
        ax[0].set_xlim(*xlimits)

        # --- Resonance Plot ---
        for i in range(len(potentials)):
            ax[i+1].plot(xData, acceleration[i], color='green',linestyle='--')
            ax[i+1].plot(xData, centerEnergy[i], label=f'$\phi =$ {potentials[i]}')
            ax[i+1].plot(xData, deceleration[i], color='red', linestyle='--')
            ax[i+1].set_xlim(*xlimits)
            ax[i+1].set_ylabel('Parallel Energy [keV]')
            ax[i+1].grid()
            ax[i+1].legend()

        ax[len(potentials)].set_xlabel(plotxlabel)
        plt.show()









# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
main(AltLow,AltHigh)