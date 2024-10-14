# --- diffNFlux_fitting.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: using the method outline in Kaeppler's thesis, we can fit inverted-V distributions
# to get estimate the magnetospheric temperature, density and electrostatic potential that accelerated
# our particles
import numpy as np

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from ACESII_code.myImports import *
from my_matplotlib_Assets.colorbars.apl_rainbow_black0 import apl_rainbow_black0_cmap
from ACESII_code.Science.InvertedV.Evans_class_var_funcs import diffNFlux_for_mappedMaxwellian
from functools import partial
from scipy.integrate import simpson,trapezoid
plt.rcParams["font.family"] = "Arial"
start_time = time.time()
# --- --- --- --- ---


print(color.UNDERLINE + f'diffNFlux_fitting' + color.END)

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# --- Physics Toggles ---
invertedV_targetTimes = [dt.datetime(2022, 11, 20, 17, 24, 56, 000000), dt.datetime(2022,11,20,17,25,3,000000)]
invertedV_TargetTimes_data = [[dt.datetime(2022, 11, 20, 17, 25, 1, 160000), dt.datetime(2022,11,20,17,25,1,500000)]]

# ---  Density, Temperature and Potential ---
threshEngy = 140
invertedV_fitDensityTempPotential = True
wPitchsToFit = [1, 2, 3, 4]
countNoiseLevel = 4



# --- Plot toggles - General ---
figure_width = 10 # in inches
figure_height =8 # in inches
Title_FontSize = 20
Label_FontSize = 20
Label_Padding = 8
Tick_FontSize = 12
Tick_Length = 1
Tick_Width = 1
Tick_FontSize_minor = 10
Tick_Length_minor = 1
Tick_Width_minor = 1
Plot_LineWidth = 0.5
plot_MarkerSize = 14
Legend_fontSize = 14
dpi = 200

# --- Cbar ---
mycmap = apl_rainbow_black0_cmap()
cbarMin, cbarMax = 1E-18, 1E-14
cbarTickLabelSize = 14




##########################
# --- --- --- --- --- ---
# --- LOADING THE DATA ---
# --- --- --- --- --- ---
##########################
prgMsg('Loading Data')
from ACESII_code.Science.InvertedV.Evans_class_var_funcs import loadDiffNFluxData
diffNFlux,Epoch,Energy,Pitch = loadDiffNFluxData()
Done(start_time)

for pitchVal in wPitchsToFit:
    ###############################
    # --- --- --- --- --- --- --- -
    # --- FIT THE DISTRIBUTIONS ---
    # --- --- --- --- --- --- --- -
    ###############################

    # define my function at the specific pitch angle of interest
    fitFuncAtPitch = partial(diffNFlux_for_mappedMaxwellian, alpha=Pitch[pitchVal])

    ################################
    # --- DEFINE THE NOISE LEVEL ---
    ################################
    rocketAttrs, b, c = ACES_mission_dicts()
    diffNFlux_NoiseCount = np.zeros(shape=(len(Energy)))
    geo_factor = rocketAttrs.geometric_factor[0]
    count_interval = 0.8992E-3
    for engy in range(len(Energy)):
        deltaT = (count_interval) - (countNoiseLevel * rocketAttrs.deadtime[0])
        diffNFlux_NoiseCount[engy] = (countNoiseLevel) / (geo_factor[0] * deltaT * Energy[engy])

    ##############################
    # --- COLLECT THE FIT DATA ---
    ##############################

    paramTime = []
    modeled_T = []
    modeled_V = []
    modeled_n = []
    modeled_ChiVal = []


    for timeset in invertedV_TargetTimes_data:
        # collect the data to fit one single pitch
        low_idx, high_idx = np.abs(Epoch - timeset[0]).argmin(), np.abs(Epoch - timeset[1]).argmin()
        EpochFitData = Epoch[low_idx:high_idx + 1]
        fitData = diffNFlux[low_idx:high_idx+1, pitchVal, :]


        # for each slice in time, loop over the data and identify the peak differentialNumberFlux (This corresponds to the
        # peak energy of the inverted-V since the location of the maximum number flux tells you what energy the low-energy BULk got accelerated to)
        # Note: The peak in the number flux is very likely the maximum value AFTER 100 eV, just find this point

        for tmeIdx in range(len(EpochFitData)):
            try:

                # --- Determine the peak point based on a treshold limit ---
                EngyIdx = np.abs(Energy - threshEngy).argmin()
                peakDiffNVal = fitData[tmeIdx][:EngyIdx].max()
                peakDiffNVal_index = np.argmax(fitData[tmeIdx][:EngyIdx])
                parallelPotential = Energy[peakDiffNVal_index]

                # ---  get the subset of data to fit to and fit it. Only include data with non-zero points ---
                xData_fit = np.array(Energy[:peakDiffNVal_index+1])
                yData_fit = np.array(fitData[tmeIdx][:peakDiffNVal_index+1])
                nonZeroIndicies = np.where(yData_fit!=0)[0]
                xData_fit = xData_fit[nonZeroIndicies]
                yData_fit = yData_fit[nonZeroIndicies]

                # --- --- --- --- --- --- ---
                # --- FIRST ITERATION FIT ---
                # --- --- --- --- --- --- ---
                deviation = 0.18
                guess = [1, 120, 250]  # observed plasma at dispersive region is 0.5E5 cm^-3 BUT this doesn't make sense to use as the kappa fit since the kappa fit comes from MUCH less dense populations above
                boundVals = [[0.001, 5], # n [cm^-3]
                             [10, 300], # T [eV]
                             [(1-deviation)*parallelPotential, (1+deviation)*parallelPotential]]  # V [eV]

                bounds = tuple([[boundVals[i][0] for i in range(len(boundVals))], [boundVals[i][1] for i in range(len(boundVals))]])
                params, cov = curve_fit(fitFuncAtPitch,xData_fit,yData_fit,maxfev=int(1E9), bounds=bounds)

                # --- --- --- --- --- --- ---
                # --- SECOND ITERATION FIT ---
                # --- --- --- --- --- --- ---
                # --- Determine the number flux at this time for pitch angles 0 to 90 deg and energies above the parallel potential ---
                PitchIndicies = [1, 2, 3, 4, 5, 6, 7, 8, 10]
                allPitchFitData = diffNFlux[low_idx + tmeIdx, 1:PitchIndicies[-1]+1, :peakDiffNVal_index]
                numberFlux_perPitch = []

                for idx, ptch in enumerate(PitchIndicies):  # integrate over energy
                    dataValues = allPitchFitData[idx]
                    dataValues_cleaned = np.array([val if val > 0 else 0 for val in dataValues])  # Clean up the data by setting fillvalues to == 0
                    EnergyValues = Energy[:peakDiffNVal_index]
                    numberFlux_perPitch.append(simpson(y=np.cos(np.radians(Pitch[ptch])) * np.sin(np.radians(Pitch[ptch])) * dataValues_cleaned, x=EnergyValues))

                # integrate over pitch
                numberFlux = (2*2*np.pi)*simpson(y=-1*np.array(numberFlux_perPitch), x=Pitch[PitchIndicies])

                # define the new fit parameters
                # n0Guess = (1E-2) * np.sqrt((2*np.pi*m_e*params[1])/q0) * numberFlux / params[2] # simple verison
                # print(numberFlux, n0Guess)
                betaChoice = 6
                n0Guess = (-1E-6)*numberFlux*np.power(q0*params[1]/(2*np.pi*m_e),-0.5)*np.power(betaChoice*(1 - (1 - 1/betaChoice)*np.exp(-(0-parallelPotential)/(params[1]*(betaChoice-1)))),-1)
                print(n0Guess)

                # Second iterative fit
                deviation = 0.18
                guess = [n0Guess, 120, 250]  # observed plasma at dispersive region is 0.5E5 cm^-3 BUT this doesn't make sense to use as the kappa fit since the kappa fit comes from MUCH less dense populations above
                boundVals = [[0.0001, 5],  # n [cm^-3]
                             [10, 300],  # T [eV]
                             [100,400]]  # V [eV]

                bounds = tuple([[boundVals[i][0] for i in range(len(boundVals))],
                                [boundVals[i][1] for i in range(len(boundVals))]])
                params, cov = curve_fit(fitFuncAtPitch, xData_fit, yData_fit, maxfev=int(1E9), bounds=bounds)



                # --- Fit the second iteration fit ---
                fittedX = np.linspace(xData_fit.min(), xData_fit.max(), 100)
                fittedY = fitFuncAtPitch(fittedX,*params)



                # --- Calculate ChiSquare ---
                ChiSquare = (1/(3-1))*sum([(fitFuncAtPitch(xData_fit[i],*params) - yData_fit[i])**2 / (yData_fit[i]**2) for i in range(len(xData_fit))])


                # Plot the result to see
                fig, ax = plt.subplots(2)
                fig.set_size_inches(figure_width, figure_height)
                fig.suptitle(f'Pitch Angle = {Pitch[pitchVal]} \n {EpochFitData[tmeIdx]} UTC', fontsize=Title_FontSize)
                cmapObj =ax[0].pcolormesh(EpochFitData,Energy,fitData.T, vmin=9E4, vmax=1E7,cmap='turbo', norm='log')
                ax[0].set_yscale('log')
                ax[0].set_ylabel('Energy [eV]', fontsize=Label_FontSize)
                ax[0].set_xlabel('Time', fontsize=Label_FontSize)
                ax[0].axvline(EpochFitData[tmeIdx],color='black', linestyle='--')
                ax[0].set_ylim(28,1000)
                plt.colorbar(cmapObj)


                ax[1].plot(Energy, fitData[tmeIdx][:],'-o')
                ax[1].plot(fittedX, fittedY, color='red', label=f'n = {round(params[0],1)}' +' cm$^{-3}$' +f'\n T = {round(params[1],1)} eV\n' + f'V = {round(params[2],1)} eV\n' + r'$\chi^{2}_{\nu}= $' + f'{round(ChiSquare,3)}')

                for i in range(2):
                    ax[i].tick_params(axis='y', which='major', colors='black', labelsize=Tick_FontSize, length=Tick_Length, width=Tick_Width)
                    ax[i].tick_params(axis='y', which='minor', colors='black', labelsize=Tick_FontSize_minor, length=Tick_Length_minor, width=Tick_Width_minor)

                ax[1].axvline(Energy[peakDiffNVal_index],color='red')
                ax[1].set_yscale('log')
                ax[1].set_xscale('log')
                ax[1].set_xlabel('Energy [eV]', fontsize=Label_FontSize)
                ax[1].set_ylabel('diffNFlux [cm$^{-2}$s$^{-1}$str$^{-1}$ eV/eV]', fontsize=Label_FontSize-4)
                ax[1].set_xlim(28,1E4)
                ax[1].set_ylim(1E4, 1E7)

                # plot the noise
                ax[1].plot(Energy, diffNFlux_NoiseCount, color='black', label=f'{countNoiseLevel}-count noise')

                ax[1].legend(fontsize=Legend_fontSize)

                plt.savefig(rf'C:\Data\ACESII\science\invertedV\TempDensityPotential_Fitting\FitData_Pitch{Pitch[pitchVal]}_{tmeIdx}.png')
                plt.close()

                # Store the data to be plotted later
                paramTime.append(EpochFitData[tmeIdx])
                modeled_n.append(params[0])
                modeled_T.append(params[1])
                modeled_V.append(params[2])
                modeled_ChiVal.append(ChiSquare)
            except:
                print('no Pitch Data')


    # --- output Plot of the time series of modeled data ---
    fig, ax = plt.subplots(4, sharex=True)
    fig.set_size_inches(figure_width, figure_height)
    fig.suptitle(f'Pitch = {Pitch[pitchVal]}',fontsize=Label_FontSize)

    # Density
    ax[0].plot(paramTime,modeled_n,marker='o')
    avg_n = round(sum(modeled_n)/len(modeled_n),1)
    ax[0].axhline(avg_n,color='red',label=rf'n (Avg) = {avg_n}')# plot the average value
    ax[0].legend(fontsize= Legend_fontSize)
    ax[0].set_ylim(1,10)
    ax[0].set_ylabel('n  [cm$^{-3}$]')

    # Temperature
    ax[1].plot(paramTime,modeled_T,marker='o')
    avg_T = round(sum(modeled_T)/len(modeled_T),1)
    ax[1].axhline(avg_T,color='red',label=rf'T (Avg) = {avg_T} ')# plot the average value
    ax[1].legend(fontsize= Legend_fontSize)
    ax[1].set_ylim(50,200)
    ax[1].set_ylabel('T [eV]')

    # Potential
    ax[2].plot(paramTime,modeled_V,marker='o')
    avg_V = round(sum(modeled_V)/len(modeled_V),1)
    ax[2].axhline(avg_V,color='red',label=rf'V (Avg) = {avg_V}')# plot the average value
    ax[2].legend(fontsize= Legend_fontSize)
    ax[2].set_ylim(100,350)
    ax[2].set_ylabel('V [eV]')


    #ChiSquare
    ax[3].plot(paramTime, modeled_ChiVal, marker='o')
    avg_Chi = round(sum(modeled_ChiVal) / len(modeled_ChiVal), 1)
    ax[3].axhline(avg_Chi, color='red', label=rf'$\chi$ (Avg) = {avg_Chi}')  # plot the average value
    ax[3].legend(fontsize=Legend_fontSize)
    ax[3].set_ylim(0, 5)
    ax[3].set_ylabel(r'$\chi ^{2}_{\nu}$')

    plt.savefig(rf'C:\Data\ACESII\science\invertedV\TempDensityPotential_Fitting\Parameters_Pitch{Pitch[pitchVal]}.png')

