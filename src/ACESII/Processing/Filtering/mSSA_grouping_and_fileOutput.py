# --- mSSA_grouping_and_fileOutput.py ---
# --- Author: C. Feltman ---
# DESCRIPTION:
# The output are SSA component files found in \L3\SSAcomponents_B
# or deltaB/deltaE files.
import matplotlib.pyplot as plt
import numpy as np

# it this code must be able to do the following SIMPLY
# [3] Plot the component grouping step
# [4] produce the deltaB/deltaE datafiles through only a couple toggles

# --- bookkeeping ---
from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---

# Just print the names of files
justPrintFileNames = False

# --- Select the Rocket ---
# 4 -> ACES II High Flier
# 5 -> ACES II Low Flier
wRocket = 4

# --- Select the DataSet ---
inputPath_modifier = 'L3\SSAcomponents_B'

#################
# --- TOGGLES ---
#################
figure_width = 12 # in inches
figure_height = 8 # in inches
Title_FontSize = 15
Label_FontSize = 15
Label_Padding = 8
Line_LineWidth = 2.5
Text_Fontsize = 25
Tick_FontSize = 25
Tick_FontSize_minor = 20
Tick_Length = 10
Tick_Width = 2
Tick_Length_minor = 5
Tick_Width_minor = 1
Plot_LineWidth = 0.5
Plot_MarkerSize = 14
Legend_fontSize = 20
dpi = 400

# --- Plotting ---
plot_GroupingSSA = True
plot_wSSAComponetFile = 2
plot_wAxesSSA = 0

# --- OUTPUT DATA ---
outputDataIntoOneMasterFile = False
# -------------------

#################
# --- IMPORTS ---
#################
from numpy.fft import rfft, fftfreq
from myspaceToolsLib.filter import generateGrouping
from myspaceToolsLib.CDF_load import getInputFiles

def mSSA_grouping_and_fileOutput(wRocket, rocketFolderPath, justPrintFileNames):

    # --- ACES-II Attributes Data ---
    rocketAttrs, b, c = ACES_mission_dicts()
    rocketID = rocketAttrs.rocketID[wRocket-4]
    inputFiles, input_names, input_names_searchable = getInputFiles(rocketFolderPath=rocketFolderPath,wRocket=wRocket,inputPath_modifier=inputPath_modifier)

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:95s}{:5.1f} MB'.format(i, file.rsplit('\\')[-1], round(getsize(file) / (10 ** 6), 1)))
        return

    # --- Get the Input Data for one File to read its contents---
    data_dict, globalAttrs = loadDictFromFile(inputFiles[plot_wSSAComponetFile], getGlobalAttrs=True)

    # --- Determine which coordinate system the data uses ---
    compNames, coordSys, coordSet = getCoordinateKeys(data_dict)
    varName = compNames[0].replace(coordSet[0], '')

    wInstr = ''
    for vari in ['RingCore', 'PoyntingFlux', 'mBeEr', 'BrEe', 'EFI']:
        if vari in inputFiles[plot_wSSAComponetFile]:
            wInstr = vari
            break

    # DETERMINE SSA Window Length
    firstIndex = input_names_searchable[plot_wSSAComponetFile].find('_WL')
    lastIndex = input_names_searchable[plot_wSSAComponetFile].find('_subset')
    SSA_window_Size = int(input_names_searchable[plot_wSSAComponetFile][firstIndex+3:lastIndex])

    # DETERMINE if data is mirror
    MirrorData = True if 'mirrored' in input_names_searchable[plot_wSSAComponetFile] else False


    #############################
    # --- DEFINE THE GROUPING ---
    #############################
    groupings_dict = {
        # 'RingCore_high': [[i for i in range(16)] +
        #                                 [20, 21, 22, 23, 24, 25, 26, 28, 29, 30, 31, 32, 37, 38, 39, 40, 41, 43, 44, 45,
        #                                  51, 52, 53, 54, 56, 57, 58, 59] +
        #                                 [60 + i for i in range(10)] +
        #                                 [71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 85, 86, 88, 89] +
        #                                 [90, 91, 92, 93, 94, 95, 96, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107,
        #                                  108, 109, 110] +
        #                                 [],  # identified harmonics (10,11,48,51,81,82 maybe?)
        #                                 [[i] for i in range(5)],  # indicies to investigation
        #                                 200],
                      'RingCore_high': [[0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,18,19,25,41,42,43,44, 58], # identified harmonics (2,3, 13 maybe not?)
                                        [[i] for i in range(60, 64)],  # indicies to investigation
                                        200],
                      # limit of the noise, between this value and 3*SSA_window size is discarded as noise
                      'E_Field_low': [[0,1,2,3,4,5,12,13],
                                      [[i] for i in range(0,5)],
                                      3*SSA_window_Size],
                      # 'RingCore_low': [[0,1,4,5],
                      #                  [[i] for i in range(0,5)],
                      #                  3*SSA_window_Size],
                    'RingCore_low': [[0,1,4,5],
                                         [[2],[4],[6]],
                                         SSA_window_Size],
                      'PoyntingFlux_low': [[4, 5, 6, 7, 8, 9, 11, 12, 13],
                                           [[i] for i in range(18, 21)],
                                           3*SSA_window_Size]
                      }


    # generate the groupings variables

    groupings = generateGrouping(SSA_window_Size = SSA_window_Size,
                                 badCompIndicies = groupings_dict[wInstr+f'_{fliers[wRocket-4]}'][0],
                                 InvestigateIndicies = groupings_dict[wInstr+f'_{fliers[wRocket-4]}'][1],
                                 noiseLimit = groupings_dict[wInstr+f'_{fliers[wRocket-4]}'][2])


    ###############################
    # --- PLOTTING THE GROUPING ---
    ###############################


    if plot_GroupingSSA:

        prgMsg('Grouping mSSA elements')

        # some preamble
        spinFreq = 0.6442441031179716 if wRocket == 4 else 0.55
        Epoch = data_dict['Epoch'][0]
        SSA_vec = np.array([data_dict[compNames[i]][0] for i in range(len(compNames))])

        # --- --- --- --- --- --- --- --- --- --- --- --- --- -
        # --- Plot the FFT and groupings for one wAxes axes ---
        # --- --- --- --- --- --- --- --- --- --- --- --- --- -
        fig, ax = plt.subplots(nrows=len(groupings), ncols=2)
        fig.set_size_inches(figure_width, figure_height)

        rktName = 'ACESII 36359\n' if wRocket == 4 else 'ACESII 36364\n'
        fig.suptitle(rktName +
                     str(compNames[plot_wAxesSSA]) +
                     f'\n Window Length: {SSA_window_Size}, Time: ' +
                     f'{Epoch[0].strftime("%H:%M:%S")} to {Epoch[-1].strftime("%H:%M:%S")}',fontsize=Title_FontSize)

        origMax = 1
        origData = [[], []]
        physData = [[], []]

        # loop over all the groups in grouping
        for i in range(len(groupings)):
            data = np.array(SSA_vec[plot_wAxesSSA])

            # combine the groupings
            plotThisData = np.array([0 for i in range(len(data[:, 0]))], dtype='float64')

            for j, componentIndex in enumerate(groupings[i]):
                plotThisData += np.array(data[:, componentIndex])

            # plot the component
            ax[i, 0].plot([i for i in range(len(plotThisData))], plotThisData)

            # if the data is mirrored, plot vertical bars on the LHS plots showing the cutoff
            if MirrorData:
                # determine where the true data starts using the length of the Epoch Variable
                startPoint = int(0.2 * len(Epoch))
                endPoint = len(SSA_vec[0]) - int(0.2 * len(Epoch))
                ax[i, 0].axvline(x=startPoint, color='black')
                ax[i, 0].axvline(x=endPoint, color='black')
            if i == 0:
                ax[i, 0].set_ylabel('Orig.', fontsize=Label_FontSize)
                ax[i, 1].set_ylabel('Orig.', fontsize=Label_FontSize)
                ax[i, 1].set_xticks([])
            elif i == 1:
                ax[i, 0].set_ylabel('$f_{n}$', fontsize=Label_FontSize)
                ax[i, 1].set_ylabel('$f_{n}$', fontsize=Label_FontSize)
                ax[i, 1].set_xticks([])
            elif i == len(groupings) - 2:
                ax[i, 0].set_ylabel('Noise', fontsize=Label_FontSize)
                ax[i, 1].set_ylabel('Noise', fontsize=Label_FontSize)
                ax[i, 1].set_xticks([])
            elif i == len(groupings) - 1:
                ax[i, 0].set_ylabel('Phys.\n Signal', fontsize=Label_FontSize)
                ax[i, 1].set_xlabel('Frequency [Hz]', fontsize=Label_FontSize)
                ax[i, 1].set_ylabel('FFT', fontsize=Label_FontSize)

                # set the ticks
                if MirrorData:
                    tickLocs = np.linspace(startPoint,endPoint,4,dtype='int64',endpoint=False)
                    tickLabels = [Epoch[indx-startPoint].strftime("%H:%M:%S") for indx in tickLocs]
                    ax[i, 0].set_xticks(ticks=tickLocs,labels=tickLabels)

            else:
                ax[i, 0].set_ylabel('F{}'.format(i), fontsize=Label_FontSize)
                ax[i, 1].set_ylabel(groupings[i])


            # calculate the FFT and plot it BUT ONLY for the real data, NOT the mirrored data
            N, T = len(plotThisData[startPoint:endPoint]), 1 / 128
            yf, xf = rfft(plotThisData[startPoint:endPoint]), fftfreq(N, T)[:N // 2]
            FFT = (2.0 / N * np.abs(yf[0:N // 2]))
            if i == 0:
                origMax = FFT.max()
                origData = [plotThisData,FFT]
            elif i == len(groupings) - 1:
                physData = [plotThisData, FFT]

            # ax[i, 1].plot(xf, FFT/origMax)
            ax[i, 1].plot(xf, FFT)
            ax[i, 1].vlines([spinFreq * (i + 1) for i in range(50)], ymin=0, ymax=1, alpha=0.25, color='red')
            ax[i, 1].set_xlim(-0.1, 15)
            ax[i, 1].set_ylim(0, 1)

        plt.savefig(r'C:\Users\cfelt\OneDrive\Desktop\THESIS\ThesisPhotos\LF_mSSA_example')
        plt.show()

        fig, ax = plt.subplots(ncols=2, nrows=3)
        # ax[0, 0].plot([i for i in range(len(plotThisData))], origData[0]/np.abs(max(origData[0])))
        ax[0, 0].plot([i for i in range(len(plotThisData))], origData[0] )
        ax[0, 0].set_ylabel('orignal')
        # ax[0, 1].plot(xf, origData[1]/max(origData[1]))
        ax[0, 1].plot(xf, origData[1])
        ax[0, 1].set_ylabel('$f_{n}$')

        # ax[1, 0].plot([i for i in range(len(plotThisData))], physData[0]/np.abs(max(origData[0])))
        ax[1, 0].plot([i for i in range(len(plotThisData))], physData[0] )
        ax[1, 0].set_ylabel('Phys. \nSignal')
        # ax[1, 1].plot(xf, physData[1]/max(origData[1]))
        ax[1, 1].plot(xf, physData[1])
        ax[1, 1].set_ylabel('$f_{n}$')

        # ax[2, 0].plot([i for i in range(len(plotThisData))], (origData[0]-physData[0])/np.abs(max(origData[0])))
        ax[2, 0].plot([i for i in range(len(plotThisData))], (origData[0] - physData[0]))
        ax[2, 0].set_ylabel('Difference')
        # ax[2, 1].plot(xf, (origData[1]-physData[1])/max(origData[1]))
        ax[2, 1].plot(xf, (origData[1] - physData[1]))
        ax[2, 1].set_ylabel('$f_{n}$')
        for k in range(3):
            ax[k, 0].set_xlim(startPoint,endPoint)
            ax[k, 0].set_ylim(-1, 1)
            ax[k, 1].set_xlim(0,15)
            ax[k, 1].set_ylim(0, 1)
        plt.show()



    ##############################
    # --- OUTPUT THE DATA DICT ---
    ##############################

    if outputDataIntoOneMasterFile:
        prgMsg('Outputting Data')

        # --- load all the data ---
        allData = []
        for file in inputFiles:
            dict =deepcopy(loadDictFromFile(file))
            allData.append(dict)

        # --- apply the mSSA grouping to the data ---
        newData = [[] for i in range(len(compNames))] # should look like [[],[],[]]

        # loop through all the dictionaries and apple the groupings
        for i, SSAdict in enumerate(allData):

            for k in range(len(compNames)):  # loop through all the components, but only take the last grouping, which is the physical signal

                data = np.array(SSAdict[compNames[k]][0])

                # combine the groupings of the last group set only
                formattedData = np.array([0 for i in range(len(data[:, 0]))], dtype='float64')
                for j, componentIndex in enumerate(groupings[-1]):
                    formattedData += np.array(data[:, componentIndex])

                newData[k].append(formattedData)

        # --- eliminate the 20% padded data on the sides ---
        if MirrorData:

            # loop through all the subset dictionaries
            for i, SSAdict in enumerate(allData):

                # for each dictionary, determine what the start/end point of the real data is
                dataTrueLength = len(SSAdict['Epoch'][0])
                startPoint = int(0.2*dataTrueLength)
                endPoint = startPoint + dataTrueLength

                # apply the reduction to the subsets
                for k in range(len(compNames)):
                    newData[k][i] = newData[k][i][startPoint:endPoint]

        # --- combine all the subsets into one dataset ---
        stitchedData = np.array([[item for sublist in dataset for item in sublist] for dataset in newData])
        stitchedEpoch = np.array([item for dict in allData for item in dict['Epoch'][0]])

        data_dict_output = {'Epoch': [stitchedEpoch, allData[0]['Epoch'][1]],
                            compNames[0]:[np.array(stitchedData[0]), allData[0][compNames[0]][1]],
                            compNames[1]:[np.array(stitchedData[1]), allData[0][compNames[1]][1]],
                            compNames[2]:[np.array(stitchedData[2]), allData[0][compNames[2]][1]]}

        # --- output the data ---
        outputPath = f'{rocketFolderPath}\\L3\\delta{varName}\\{fliers[wRocket-4]}\\ACESII_{rocketID}_{wInstr}_{coordSys}_WL{SSA_window_Size}.cdf'
        outputCDFdata(outputPath, data_dict_output)
        Done(start_time)









# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
rocketFolderPath = ACES_data_folder
if len(glob(f'{rocketFolderPath}{inputPath_modifier}\{fliers[wRocket-4]}\*.cdf')) == 0:
    print(color.RED + 'There are no .cdf files in the specified directory' + color.END)
else:
    mSSA_grouping_and_fileOutput(wRocket, rocketFolderPath, justPrintFileNames)
