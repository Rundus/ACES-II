# --- newCDF4_to_cdf.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: converts EISCAT data from netCDF4 to .cdf



# --- bookkeeping ---
# !/usr/bin/env python
__author__ = "Connor Feltman"
__date__ = "2022-08-22"
__version__ = "1.0.0"

from src.my_imports import *
start_time = time.time()
# --- --- --- --- ---


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
justPrintFileNames = False
wFiles = [0]
inputPath_modifier = r'\science\EISCAT\tromso\UHF'
outputPath_modifier = r'\science\EISCAT\tromso\UHF' # e.g. 'L2' or 'Langmuir'. It's the name of the broader output folder
fixEpoch = True
# ---------------------------
outputData = True
# ---------------------------
# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import netCDF4


def netCDF4_to_cdf(rocketFolderPath, justPrintFileNames, wFile):

    inputFiles = glob(f'{rocketFolderPath}{inputPath_modifier}\*.nc')
    input_names = [ifile.replace(f'{rocketFolderPath}{inputPath_modifier}\\', '') for ifile in inputFiles]
    input_names_searchable = [ifile.replace(inputPath_modifier.lower() +'_', '') for ifile in input_names]

    if justPrintFileNames:
        for i, file in enumerate(inputFiles):
            print('[{:.0f}] {:80s}{:5.1f} MB '.format(i, input_names_searchable[i], round(getsize(file) / (10 ** 6), 1)))
        return

    print('\n')

    # --- --- --- --- --- -
    # --- LOAD THE DATA ---
    # --- --- --- --- --- -

    # --- get the data from the netCDF4file ---
    stl.prgMsg(f'Loading data from {inputPath_modifier} Files')
    cdf4File = netCDF4.Dataset(inputFiles[wFile],"r", format="NETCDF4")
    stl.Done(start_time)

    # --- Global Attributes ---
    globalAttrs = cdf4File.__dict__
    varAttrs =cdf4File.variables['power'].ncattrs()
    var = cdf4File.variables['power']


    # --- Variables ---
    ExampleVarAttrs = {'FIELDNAM': None,
                       'LABLAXIS': None,
                       'DEPEND_0': None,
                       'DEPEND_1': None,
                       'DEPEND_2': None,
                       'FILLVAL': None,
                       'FORMAT': None,
                       'UNITS': None,
                       'VALIDMIN': None,
                       'VALIDMAX': None,
                       'VAR_TYPE': 'data',
                       'SCALETYP': 'linear',
                       'description':None}
    fileVars = cdf4File.variables
    varKeys = [str(key) for key in fileVars.keys()]
    data_dict = {}
    for key in varKeys:
        var = fileVars[key]
        varAttrs = deepcopy(ExampleVarAttrs)
        if key not in ['timestamps','range']:
            varAttrs['DEPEND_0'] = 'timestamps'
            varAttrs['DEPEND_1'] = 'range'
        varAttrs['UNITS'] = str(var.units)
        varAttrs['FIELDNAM'] = str(var.units)
        varAttrs['LABLAXIS'] = str(var.description)
        varAttrs['description'] = str(var.description)
        varAttrs['FILLVAL'] = 9.969209968386869e+36
        data_dict = {**data_dict, **{key: [np.array(list(var[:])),varAttrs]} }

    if fixEpoch:
        # for EISCAT data, the timestamps are from Jan 1st, 1970. We convert them to tt2000 here
        secondsSince2000 =  24*3600*pycdf.lib.epoch_to_num(pycdf.lib.datetime_to_epoch(pycdf.lib.tt2000_to_datetime(0)))
        Epoch = np.array([  pycdf.lib.tt2000_to_datetime(int((data_dict['timestamps'][0][i] - secondsSince2000)*1E9)) for i in range(len(data_dict['timestamps'][0]))])

        # replace old data
        data_dict["Epoch"] = data_dict.pop("timestamps")
        data_dict['Epoch'][0] = Epoch

        for key in data_dict.keys():
            data_dict[key][1]['DEPEND_0'] = 'Epoch'
    # --- --- --- --- --- --- ---
    # --- WRITE OUT THE DATA ---
    # --- --- --- --- --- --- ---

    if outputData:
        stl.prgMsg('Creating output file')

        outputPath = f'{rocketFolderPath}{outputPath_modifier}\{input_names_searchable[wFile].replace(".nc",".cdf")}'
        stl.outputCDFdata(outputPath, data_dict,globalAttrsMod=globalAttrs, )

        stl.Done(start_time)





# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---

if len(glob(f'{DataPaths.ACES_data_folder}{inputPath_modifier}\*.nc')) == 0:
    print(stl.color.RED + 'There are no .cdf files in the specified directory' + stl.color.END)
else:
    if wFiles == []:
        for fileNo in  range(len(glob(f'{DataPaths.ACES_data_folder}{inputPath_modifier}\*.nc'))):
            netCDF4_to_cdf(DataPaths.ACES_data_folder, justPrintFileNames, fileNo)
    else:
        for fileNo in wFiles:
            netCDF4_to_cdf(DataPaths.ACES_data_folder, justPrintFileNames, fileNo)
