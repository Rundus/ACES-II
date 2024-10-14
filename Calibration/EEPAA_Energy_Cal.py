# --- EEPAA_Energy_Cal.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Using the 401 chamber calibration data, fit a linear curve to
# the data to determine a calibration curve for EEPAA steps/energy
import matplotlib.pyplot as plt
import numpy as np

# --- --- --- --- ---
from ACESII_code.myImports import *

start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
inputFile = r'C:\Users\cfelt\PycharmProjects\UIOWA_CDF_operator\ACESII_code\Calibration\EEPAA_response_comparison.xlsx'

wRow = 0 # row corresponding to where all the names of the variables are


plotCalData_vs_Theoretical = True

# Thereotical EEPAA Energies
# 12750.00, 10905.65, 9328.10, 7978.75, 6824.58, 5837.38, 4992.97, 4270.72,
#                           3652.94,  3124.52, 2672.55, 2285.95, 1955.28, 1672.44, 1430.51, 1223.58,
#                           1046.58,   895.19,  765.70,  654.94,  560.20,  479.16,  409.85,  350.56,
#                            299.85,   256.48,  219.38,  187.64,  160.50,  137.28,  117.42,  100.44,
#                             85.91,    73.48,   62.85,   53.76,   45.98,   39.33,   33.64,   28.78,
#                             24.61,    21.05,   18.01,   15.40,   13.17,   11.27,    9.64,    8.24, 7.05


# --- --- --- ---
# --- import ---
# --- --- --- ---

import pandas as pd


def EEPAA_Energy_Cal(inputFile):
    # collect the csv data
    df = pd.read_excel(inputFile)
    steps = np.array(df.STEP)
    Energy_Theory = np.array(df['Energy (Theoretical) [eV]'])
    Energy_cal = np.array(df['Calibration Estimate [eV]'])
    dataForCal = np.array([list(pair) for pair in zip(steps,Energy_cal) if not (np.isnan(pair[0]) or np.isnan(pair[1]))])

    xData = dataForCal[:, 0]
    yData = np.log(np.array(dataForCal[:, 1]))
    def fitFunc(x,A,B):
        return A*x + B

    params, cov = curve_fit(fitFunc, xData, yData)
    xDataFit = np.linspace(min(xData), max(xData), 200)
    yDataFit = np.array([fitFunc(x,*params) for x in xDataFit])

    if plotCalData_vs_Theoretical:
        fig, ax = plt.subplots()
        fig.suptitle('Measured E vs Step #')
        ax.scatter(xData, yData,label='Raw',color='blue')
        ax.plot(xDataFit, yDataFit, label=f'Fit Function: Ln(E_measured) = A*(Step #) + B\n A={params[0]} \n B = {params[1]}',color='red')
        ax.set_ylabel('$Ln(E_{measured})$ ')
        ax.set_xlabel('Step #')
        ax.legend()
        plt.show()

    np.set_printoptions(formatter={'float_kind':'{:25f}'.format})
    new_EEPAA_Energies = np.array([ round(np.exp(fitFunc(x,*params)),2) for x in range(1,50)],dtype='float64')[::-1]

    print('[',end='')
    for thing in new_EEPAA_Energies:
        print(thing,', ',end='')
    print(']',end='')




# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
EEPAA_Energy_Cal(inputFile)