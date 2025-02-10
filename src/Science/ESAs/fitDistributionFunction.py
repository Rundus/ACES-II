# --- fitDistributionFunction.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Combine the High Flyer EEPAA and LEESA data into a single, large electron dataset



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import math

from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
wPitches = [0,1,2,3,4,5,6] # indicies of pitch bins to use in plot. Nominally [0,1,2,3,4] --> -10deg to 30deg
plotDataRegion = True
useKappa = True
wTtime = 1
# ---------------------------
outputData = True
# ---------------------------


targetTimes = {
    'Quiet_Region' : [dt.datetime(2022,11,20,17,22,15,000000), dt.datetime(2022, 11, 20, 17, 22, 35, 000000)],
    'Before_Dispersive': [dt.datetime(2022,11,20,17,24,55,000000), dt.datetime(2022, 11, 20, 17, 24, 55, 500000)], # before dispersive region (First STEB occurs at dt.datetime(2022, 11, 20, 17, 24, 56, 000000))
    'Primary_Auroral': [dt.datetime(2022,11,20,17,25,30,000000), dt.datetime(2022, 11, 20, 17, 25, 35, 500000)]
}


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from numpy import sqrt,abs,exp, array,pi,linspace
from scipy.special import  gamma

def fitDistributionFunction(wRocket, rocketFolderPath):
    wflyer = wRocket-4

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    ModelData = L2_TRICE_Quick(wflyer)
    inputFile_superSetDist = glob(f'{rocketFolderPath}L3\distFunc\{fliers[wflyer]}\*eSuperSet*')[0]
    fileoutName_Dist = f'ACESII_{rocketID}_'

    # --- --- --- --- --- -
    # --- LOAD THE DATA ---
    # --- --- --- --- --- -

    # --- get the data from the file ---
    prgMsg(f'Loading Electron Files')
    # load data here
    data_dict_eSuperSet = loadDictFromFile(inputFile_superSetDist)
    Done(start_time)

    # --- Collect the Data ---
    prgMsg('Collecting Data')
    # in the region of interest, collect all the x,y points of the data [Energy, PSD]

    tTimeNames = [key for key,val in targetTimes.items()]
    timeMin, timeMax = abs(data_dict_eSuperSet['Epoch'][0] - targetTimes[tTimeNames[wTtime]][0]).argmin(),abs(data_dict_eSuperSet['Epoch'][0] - targetTimes[tTimeNames[wTtime]][1]).argmin()
    distFunc = data_dict_eSuperSet['Distribution_Function'][0][timeMin:timeMax+1]
    oneCountLevel = data_dict_eSuperSet['oneCountLevel'][0][timeMin:timeMax+1]
    Energy = data_dict_eSuperSet['Energy'][0]
    Pitch = data_dict_eSuperSet['Pitch_Angle'][0]
    xData = []
    yData = []
    xData_singleCount = []
    yData_singleCount = []
    ranges = [range(len(distFunc)), wPitches, range(len(distFunc[0][0]))]

    for tme, ptch, engy in itertools.product(*ranges):
        xData_singleCount.append(Energy[engy])
        yData_singleCount.append(oneCountLevel[tme][ptch][engy])
        if distFunc[tme][ptch][engy] != 0 and distFunc[tme][ptch][engy] != rocketAttrs.epoch_fillVal:
            xData.append(Energy[engy])
            yData.append(distFunc[tme][ptch][engy])


    # prepare the data
    xData, yData = array(xData), array(yData)
    xData_singleCount = array(xData_singleCount)
    yData_singleCount = array(yData_singleCount)
    Done(start_time)


    # --- Break Data into multi Regions for Multiple Fit ---
    R1_engyBound = 13678.4/1000
    xDataR1, yDataR1 = [], []
    xDataR2, yDataR2= [], []

    for i in range(len(xData)):
        if xData[i] <= R1_engyBound:
            xDataR1.append(xData[i])
            yDataR1.append(yData[i])
        else:
            xDataR2.append(xData[i])
            yDataR2.append(yData[i])

    # --- Define the Fitting Function ---
    def MaxwellianDist(x, n, T):  # Fits the NATURAL LOG of a Maxwellian for a Uniform/Thermalized Plasma
        return 4*pi*n*( (m_e/(2*pi*q0*T))**1.5 ) * (2*x/m_e) * exp(-x/(T))

    def KappaDist(x,kappa,n,T):

        return n * gamma(kappa+1) *( (1 + 2*q0*x/(m_e*kappa*((sqrt(q0*T*(2*kappa - 3)/(kappa*m_e)))**2)) ) **(-kappa - 1)) / ( ((sqrt(pi)*(sqrt(q0*T*(2*kappa - 3)/(kappa*m_e))))**3) * ( (kappa**1.5)  *gamma(kappa-0.5))   )

    xDataFitThis = xDataR2
    yDataFitThis = yDataR2


    if useKappa:
        # --- Low T Plasma ---
        # guess = [1.500001, 1E13, 13]
        # boundVals = [[1.5000000001, 2],[1E10, 1E14], [0.001, 1E3]]  # kappa, n, T

        # --- High T Plasma ---
        guess = [1.509, 1.3E6, 2000]
        boundVals = [[1.5000000001, 2], [1E4, 1E14], [0.001, 1E4]]  # kappa, n, T


        bounds = tuple([[boundVals[i][0] for i in range(len(boundVals))], [boundVals[i][1] for i in range(len(boundVals))]])
        params, cov = curve_fit(KappaDist, xDataFitThis, yDataFitThis, p0=guess, maxfev=int(1E9), bounds=bounds)
        # params, cov = curve_fit(KappaDist, xDataFitThis, yDataFitThis, maxfev=int(1E6))
        xDataFit = linspace(Energy.min(), Energy.max(), 1000)
        yDataFit = KappaDist(xDataFit, *params)
        print(f'kappa = {params[0]}')
        print(f'n= {params[1]} m^-3')
        print(f'T= {params[2]} eV')
    else:
        guess = [1E7, 10]
        boundVals = [[1E6, 1E15], [1E-1, 1E3]] # n, T
        bounds = tuple([[boundVals[i][0] for i in range(len(boundVals))], [boundVals[i][1] for i in range(len(boundVals))]])
        params, cov = curve_fit(MaxwellianDist, xDataFitThis, yDataFitThis, p0=guess, bounds=bounds, maxfev=int(1E6))
        xDataFit = linspace(Energy.min(),Energy.max(),1000)
        yDataFit = MaxwellianDist(xDataFit,*params)
        print(f'n= {params[0]} m^-3')
        print(f'T= {params[1]} eV')


    if plotDataRegion:
        prgMsg('Making Plot')

        fig, ax = plt.subplots()
        # ax.scatter(xDataFitThis, yDataFitThis, label='Region 1', color='blue')
        ax.scatter(xDataR1, yDataR1, label='Region 1', color='blue')
        ax.scatter(xDataR2, yDataR2, label='Region 2', color='red')
        ax.plot(xDataFit,yDataFit,label='Fit', color='black')
        ax.set_title(f'Electron PSD [Observation]\n {tTimeNames[wTtime]} \n ' + f'pa < {Pitch[wPitches[-1]]} deg   '+ f'{targetTimes[tTimeNames[wTtime]][0].strftime("%H:%M:%S")} - {targetTimes[tTimeNames[wTtime]][1].strftime("%H:%M:%S")} UTC \n' +
                     f'$\kappa =$ {params[0]}, n={params[1]/(cm_to_m**3)} $cm-3$ ,T= {params[2]}eV')
        # ax.text(x=40,
        #         y=0.7*yData.max(),
        #         s=f'pa < {Pitch[wPitches[-1]]} deg   '+ f'{targetTimes[tTimeNames[wTtime]][0].strftime("%H:%M:%S")} - {targetTimes[tTimeNames[wTtime]][1].strftime("%H:%M:%S")} UTC',
        #         ha='center',
        #         color='black',
        #         fontsize=13)
        ax.set_yscale('log')
        ax.set_ylabel('PSD [$s^{3}/m^{6}$]')
        ax.set_xscale('log')
        ax.set_ylim(1E-19, 1E-10)
        ax.grid(True)
        ax.set_xlim(1E0, 1.4E4)
        ax.set_xlabel('Energy [eV]')
        ax.legend()
        plt.show()
        Done(start_time)







# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
fitDistributionFunction(4, rocketFolderPath)
