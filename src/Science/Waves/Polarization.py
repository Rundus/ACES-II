# --- PoyntingFlux.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Determine the PoyntingFLux of the data using E-Field and B-Field Measurements.
# For the low flyer, it ONLY accepts despun data, high flyer has its own case



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

import matplotlib.pyplot as plt
import numpy as np

from ACESII_code.myImports import *

start_time = time.time()
# --- --- --- --- ---



# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---

# Just print the names of files
justPrintFileNames = False
printMagFiles = False
printElecFiles = False

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 5

modifier = ''
inputPath_modifier_elec = 'science/deltaE'
wMagFile = 1
inputPath_modifier_mag = 'science/deltaB' # e.g. 'L1' or 'L1'. It's the name of the broader input folder
wEFIFile = 1
outputPath_modifier = 'science/Polarization' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder


# --- --- --- Which Data --- -- ---
useDelta_E_B = True # use the deltaB, deltaE data
# --- --- --- reduce data --- --- ---
reduceData = True
targetTimes = [pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 24, 30, 000)),
               pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 25, 15, 000))]
# --- --- --- PLOT --- --- ---
# --- --- --- ANIMATION --- --- ---
animatePolarization = False
plotSpecificLocations = True # only plot specific inidices
frame_skips = 1
fps = 5
# --- --- --- OUTPUT --- --- ---
outputData = True


timeOffset_from_Roglands_Analysis = 0.1157 # seconds

# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
from matplotlib import animation
from ACESII_code.class_var_func import u0, IonMasses,InterpolateDataDict
from matplotlib.gridspec import GridSpec

def Polarization(wRocket, rocketFolderPath, justPrintFileNames, wflyer):

    # --- ACES II Flight/Integration Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wflyer]
    globalAttrsMod = rocketAttrs.globalAttributes[wflyer]
    globalAttrsMod['Logical_source'] = globalAttrsMod['Logical_source'] + 'L2'
    outputModelData = L2_TRICE_Quick(wflyer)

    inputFiles_elec = glob(f'{rocketFolderPath}{inputPath_modifier_elec}\{fliers[wflyer]}{modifier}\*.cdf*')
    inputFiles_mag = glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wflyer]}{modifier}\*RingCore*')

    # determine which coordinates system its in based on magnetic field data
    FileName = inputFiles_mag[wMagFile].replace('.cdf','')

    if 'ENU' in FileName:
        wCoordinate = 'ENU'
        perpComponents = [0,1]
    elif 'Field_Aligned' in FileName:
        wCoordinate = 'Field_Aligned'
        perpComponents = [0,2]

    fileoutName = f'ACESII_{rocketID}_PolarizationAngle_delta_{wCoordinate}.cdf' if useDelta_E_B else f'ACESII_PolarizationAngle_flight_{wCoordinate}.cdf'

    if justPrintFileNames:
        if printMagFiles:
            print('--- B-FIELD FILES ---')
            for i, file in enumerate(inputFiles_mag):
                print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles_mag[i], round(getsize(file) / (10 ** 6), 1)))
            print('\n')

        if printElecFiles:
            print('--- E-FIELD FILES ---')
            for i, file in enumerate(inputFiles_elec):
                print('[{:.0f}] {:80s}{:5.1f} MB'.format(i, inputFiles_elec[i], round(getsize(file) / (10 ** 6), 1)))
            print('\n')
    else:
        print('\n')
        print(color.UNDERLINE + f'Calculating Poynting flux for ACESII {rocketID}' + color.END)

        # --- get the data from the mag file ---
        prgMsg(f'Loading data from mag Files')
        data_dict_mag = loadDictFromFile(inputFiles_mag[wMagFile], {})
        data_dict_mag['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_mag['Epoch'][0]])

        # component names for the magnetic field
        compNames_mag = [ key for key, val in data_dict_mag.items() if key.lower() not in ['epoch','db_mag']]

        if reduceData:
            indicies = [np.abs(data_dict_mag['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_mag['Epoch'][0] - targetTimes[1]).argmin()]

            for key, val in data_dict_mag.items():
                data_dict_mag[key][0] = data_dict_mag[key][0][indicies[0]:indicies[1]]

        # create vector variable and convert to tesla

        Epoch_mag_tt2000 = data_dict_mag['Epoch'][0]
        Epoch_mag_dt = [pycdf.lib.tt2000_to_datetime(data_dict_mag['Epoch'][0][i]) for i in range(len(data_dict_mag['Epoch'][0]))]
        B_Field = (1E-9)*np.array([np.array([data_dict_mag[compNames_mag[0]][0][i], data_dict_mag[compNames_mag[1]][0][i], data_dict_mag[compNames_mag[2]][0][i]]) for i in range(len(data_dict_mag['Epoch'][0]))])
        Done(start_time)

        if wRocket == 4:

            # collect the Magnitude of B from L1 spun data
            prgMsg('Getting Bmag')
            inputFileBmag = glob('C:\Data\ACESII\L1\high\*RingCore_rktFrm*')
            data_dict_Bmag = loadDictFromFile(inputFileBmag[0], {})
            data_dict_Bmag['Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_Bmag['Epoch'][0]])
            Done(start_time)

            prgMsg('Getting Plasma Density')
            inputFileBmag = glob('C:\Data\ACESII\science\Langmuir\high\*Temp&Density*')
            data_dict_density = loadDictFromFile(inputFileBmag[0], {})
            data_dict_density['fixed_Epoch'][0] = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dict_density['fixed_Epoch'][0]])
            Done(start_time)

            # reduce the datasets:
            prgMsg('Reducing Data')

            # Bmag
            indicies = [np.abs(data_dict_Bmag['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_Bmag['Epoch'][0] - targetTimes[1]).argmin()]
            for key, val in data_dict_Bmag.items():
                data_dict_Bmag[key][0] = deepcopy(data_dict_Bmag[key][0][indicies[0]:indicies[1]])

            # Density
            indicies = [np.abs(data_dict_density['fixed_Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_density['fixed_Epoch'][0] - targetTimes[1]).argmin()]
            for key, val in data_dict_density.items():
                data_dict_density[key][0] = deepcopy(data_dict_density[key][0][indicies[0]:indicies[1]])

            Done(start_time)

            #################################
            # --- CALCULATE POYNTING FLUX ---
            #################################
            prgMsg('Calculating Poynting Flux using EigenFunction')

            # down sample the density data onto the magnetometer data
            indiciesDownsampled = [np.abs(data_dict_density['fixed_Epoch'][0] - data_dict_mag['Epoch'][0][i]).argmin() for i in range(len(data_dict_mag['Epoch'][0]))]
            plasmaDensity = np.array([data_dict_density['fixed_ni_density'][0][index] for index in indiciesDownsampled])

            # calculate Alfven Velocity
            AlfvenVelocity = np.array([((1E-9)*data_dict_Bmag['Bmag'][0][i])/np.sqrt(u0*plasmaDensity[i]*IonMasses[0]) for i in range(len(plasmaDensity))])


            # calculate Alfven eigenfunction E


            # calculate Poyning Flux
            S = np.array([(AlfvenVelocity[i]/(2*u0))*np.array([B_Field[i][0]**2, B_Field[i][1]**2, B_Field[i][2]**2]) for i in range(len(data_dict_mag['Epoch'][0])) ])

            Done(start_time)

            if plotSPoynting:
                prgMsg('Plotting HF Poynting Flux')
                Epoch = [ pycdf.lib.tt2000_to_datetime(tme) for tme in data_dict_mag['Epoch'][0]]
                fig, ax = plt.subplots(3)
                fig.suptitle('Alfven EigenFunction')
                for i in range(3):
                    ax[i].plot(Epoch, S[:, i])
                    ax[i].set_ylabel(compNamesS[i] + '[$W/m^{2}$]')
                    ax[i].set_ylim(-0.01, 0.15)

                plt.show()
                Done(start_time)

        elif wRocket == 5:

            # --- get the data from the electric file ---
            prgMsg(f'Loading data from Electric Field Files')
            data_dict_elec = loadDictFromFile(inputFiles_elec[wEFIFile], {})
            data_dict_elec['Epoch'][0] = np.array([int(pycdf.lib.datetime_to_tt2000(tme) + timeOffset_from_Roglands_Analysis*1E9) for tme in data_dict_elec['Epoch'][0]])
            compNames_elec = [key for key, val in data_dict_elec.items() if key.lower() not in ['epoch', 'de_mag']]

            # --- reduce each dataset ---
            if reduceData:
                indicies = [np.abs(data_dict_elec['Epoch'][0] - targetTimes[0]).argmin(), np.abs(data_dict_elec['Epoch'][0] - targetTimes[1]).argmin()]
                for key, val in data_dict_elec.items():
                    data_dict_elec[key][0] = deepcopy(data_dict_elec[key][0][indicies[0]:indicies[1]])

            # get the electric field and convert it to V/m
            E_Field = (1/1000)*np.array([ [data_dict_elec[compNames_elec[0]][0][i], data_dict_elec[compNames_elec[1]][0][i],data_dict_elec[compNames_elec[2]][0][i]] for i in range(len(data_dict_elec['Epoch'][0]))])
            Done(start_time)

            ##########################################
            # --- DOWNSAMPLE EFI ONTO MAG TIMEBASE ---
            ##########################################

            prgMsg('Downsampling EFI Data via Interpolation onto Mag Time')

            data_dict_elec_interp = InterpolateDataDict(InputDataDict=data_dict_elec,
                                                        InputEpochArray=data_dict_elec['Epoch'][0],
                                                        wKeys=[],
                                                        targetEpochArray=data_dict_mag['Epoch'][0])

            E_y = np.array([data_dict_elec_interp[compNames_elec[perpComponents[1]]][0][i] for i in range(len(Epoch_mag_tt2000))])
            E_x = np.array([data_dict_elec_interp[compNames_elec[perpComponents[0]]][0][i] for i in range(len(Epoch_mag_tt2000))])

            B_y = np.array(data_dict_mag[compNames_mag[perpComponents[1]]][0])
            B_x = np.array(data_dict_mag[compNames_mag[perpComponents[0]]][0])
            Done(start_time)

            #################################################
            # --- Determine the Polarization of the Waves ---
            #################################################

            # NOTE: In general the polarization is determined from the E-Field of the Waves using the jones vector:
            # Jones Vector (in the plane perpendicular to wave propogation) = [Ex, Ey] / E_xy_mag
            #

            # --- get the angle between the vectors ---
            E_vec = np.array([[E_x[i], E_y[i]] for i in range(len(E_x))])
            B_vec = np.array([[B_x[i], B_y[i]] for i in range(len(B_x))])
            alpha = [round(np.degrees(np.arccos(np.dot(E_vec[i], B_vec[i])/(np.linalg.norm(E_vec[i]) * np.linalg.norm(B_vec[i])))),1) for i in range(len(B_vec))]


            if animatePolarization:

                prgMsg('Creating Polarization Plot')
                ###############################################
                # --- Animate the polarization of the Waves ---
                ###############################################

                normsB = np.array([np.linalg.norm(B_vec[i]) for i in range(len(B_vec))])
                normsE = np.array([np.linalg.norm(E_vec[i]) for i in range(len(B_vec))])

                # normalize the x-y data
                B_unit = np.array([B_vec[i]/normsB.max() for i in range(len(B_vec))])
                E_unit = np.array([E_vec[i]/normsE.max() for i in range(len(E_vec))])


                # get the Epoch information
                # convert all Epoch data to datetimes to display them
                datetimesEpoch = np.array([pycdf.lib.tt2000_to_datetime(data_dict_mag['Epoch'][0][i]).strftime("%H:%M:%S:%f") for i in range(len(data_dict_mag['Epoch'][0]))])

                # --- INITIALIZE PLOT ---

                # gridspec the data

                fig = plt.figure()
                fig.set_size_inches(18.5, 10.5)
                title = fig.suptitle(f'ACESII {rocketID}\n'
                                     f'{datetimesEpoch[0]}\n'
                                     r'$\alpha$ = ' + f'{alpha[0]} deg')

                ######################
                # the primary Gridspec
                ######################
                gs0 = GridSpec(ncols=1, nrows=2, height_ratios=[0.5,0.5])

                # --- the E-B data sub-Gridspec ---
                gs_top = gs0[0].subgridspec(nrows=4,ncols=1,hspace=0)

                # plt.subplots_adjust(wspace=0.1, hspace=0)
                ax_Bx = fig.add_subplot(gs_top[0, :])
                ax_By = fig.add_subplot(gs_top[1, :], sharex=ax_Bx)
                ax_Ex = fig.add_subplot(gs_top[2, :], sharex=ax_Bx)
                ax_Ey = fig.add_subplot(gs_top[3, :], sharex=ax_Bx)

                # plot the E/B data
                ax_Bx.plot(Epoch_mag_dt, B_x)
                ax_Bx.set_ylabel(compNames_mag[perpComponents[0]])
                BxVline = ax_Bx.axvline(Epoch_mag_dt[0],ymin=-20,ymax=20,color='red',alpha=0.5)
                ax_Bx.set_ylim(-5,5)

                ax_By.plot(Epoch_mag_dt, B_y)
                ax_By.set_ylabel(compNames_mag[perpComponents[1]])
                ByVline = ax_By.axvline(Epoch_mag_dt[0], ymin=-20, ymax=20, color='red', alpha=0.5)
                ax_By.set_ylim(-5, 5)

                ax_Ex.plot(Epoch_mag_dt, E_x)
                ax_Ex.set_ylabel(compNames_elec[perpComponents[0]])
                ExVline = ax_Ex.axvline(Epoch_mag_dt[0], ymin=-20, ymax=20, color='red', alpha=0.5)
                ax_Ex.set_ylim(-5, 5)

                ax_Ey.plot(Epoch_mag_dt, E_y)
                ax_Ey.set_ylabel(compNames_elec[perpComponents[1]])
                EyVline = ax_Ey.axvline(Epoch_mag_dt[0], ymin=-20, ymax=20, color='red', alpha=0.5)
                ax_Ey.set_ylim(-5, 5)

                # --- polarization gridspec ---
                gs_bottom = gs0[1].subgridspec(nrows=1,ncols=1)

                ax_Polar = fig.add_subplot(gs_bottom[:, :])

                B_quiver = ax_Polar.quiver(0, 0, B_x[0], B_y[0], color='blue',label='B-Field')
                E_quiver = ax_Polar.quiver(0, 0, E_x[0], E_y[0], color='red', label='E-Field')
                ax_Polar.set_xlabel(compNames_mag[perpComponents[0]].replace('B_', '') + ' [East-like]')
                ax_Polar.set_ylabel(compNames_mag[perpComponents[1]].replace('B_', '') + ' [North-like]')
                ax_Polar.legend(loc='best')
                ax_Polar.set_ylim(-0.05, 0.05)
                ax_Polar.set_xlim(-0.05, 0.05)
                plt.show()


                # --- Create Animation Function ---
                def animatePlot(i):

                    # update Epoch title
                    title.set_text(f'ACESII {rocketID}\n'
                             f'{datetimesEpoch[i]}\n'
                             r'$\alpha$ = ' + f'{alpha[i]} deg')

                    # update the quiver values
                    B_quiver.set_UVC(B_unit[i][0], B_unit[i][1])
                    E_quiver.set_UVC(E_unit[i][0], E_unit[i][1])

                    # update the vline locations
                    BxVline.set_xdata([Epoch_mag_dt[i],Epoch_mag_dt[i]])
                    ByVline.set_xdata([Epoch_mag_dt[i],Epoch_mag_dt[i]])
                    ExVline.set_xdata([Epoch_mag_dt[i],Epoch_mag_dt[i]])
                    EyVline.set_xdata([Epoch_mag_dt[i],Epoch_mag_dt[i]])



                if plotSpecificLocations:
                    locations = [i for i in range(0,int(len(datetimesEpoch)), 1)]
                else:
                    locations = [i for i in range(0, int(len(datetimesEpoch)), frame_skips)]  # NEEDS TO BE THE HIGH FLYER LENGTH


                # --- Create the Animation ---
                anim = animation.FuncAnimation(fig=fig, func=animatePlot, interval=1000 / fps, frames=locations)
                anim.save(f'C:\Data\ACESII\science\polarization\\{fliers[wflyer]}\\ACESII_{rocketID}_polarization_direction.mp4')
                Done(start_time)


            # --- --- --- --- --- --- ---
            # --- WRITE OUT THE DATA ---
            # --- --- --- --- --- --- ---

            if outputData:
                prgMsg('Creating output file')

                example_attrs = {'LABLAXIS': 'deltaEdeltaB_Angle', 'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None,
                                 'FILLVAL': rocketAttrs.epoch_fillVal, 'FORMAT': 'E12.2', 'UNITS': 'deg',
                                 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'data', 'SCALETYP': 'linear'}

                # clean up Alpha data so it only displays data within +/- 20deg of 90deg
                for i,val in enumerate(alpha):
                    if np.abs(90-val) > 20:
                        alpha[i]= rocketAttrs.epoch_fillVal


                data_dict = {'E_B_angle':[np.array(alpha),example_attrs], 'Epoch': [np.array(data_dict_mag['Epoch'][0]),data_dict_mag['Epoch'][1]]}

                outputPath = f'{rocketFolderPath}{outputPath_modifier}\{fliers[wflyer]}\\{fileoutName}'

                outputCDFdata(outputPath, data_dict, L2_TRICE_Quick(1), globalAttrsMod, 'Angle_Between_deltaE_deltaB')

                Done(start_time)










# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
if wRocket == 4:  # ACES II High
    rocketFolderPath = ACES_data_folder
    wflyer = 0
elif wRocket == 5: # ACES II Low
    rocketFolderPath = ACES_data_folder
    wflyer = 1

if len(glob(f'{rocketFolderPath}{inputPath_modifier_elec}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no electric field .cdf files in the specified directory' + color.END)
elif len(glob(f'{rocketFolderPath}{inputPath_modifier_mag}\{fliers[wflyer]}\*.cdf')) == 0:
    print(color.RED + 'There are no B-field .cdf files in the specified directory' + color.END)
else:
    if justPrintFileNames:
        Polarization(wRocket, rocketFolderPath, justPrintFileNames,wflyer)
    else:
        Polarization(wRocket, rocketFolderPath, justPrintFileNames, wflyer)
