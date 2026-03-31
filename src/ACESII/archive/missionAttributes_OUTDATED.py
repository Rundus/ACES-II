# --- missionAttributes_OUTDATED.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Place to store the attributes specific to certain rockets


# --- --- --- ---
# --- IMPORTS ---
# --- --- --- ---
import numpy as np
import datetime as dt

# --- --- --- --- ---
# --- ATTRIBUTES ---
# --- --- --- --- ---

class makeRocketAttrs:
    def __init__(self, missionAttrs):
        self.missionNam = missionAttrs['missionNam']
        self.rocketID = missionAttrs['rocketID']
        self.globalAttributes = missionAttrs['globalAttributes']
        self.g_nWords = missionAttrs['g_nWords']
        self.InstrNames = missionAttrs['InstrNames']
        self.InstrNames_LC = missionAttrs['InstrNames_LC']
        self.InstrNames_Full = missionAttrs['InstrNames_Full']
        self.ESA_words = missionAttrs['ESA_words']
        self.nDataLen = missionAttrs['nDataLen']
        self.gSync = missionAttrs['gSync']
        self.ESA_CLK_INPUT = missionAttrs['ESA_CLK_INPUT']
        self.LP_words = missionAttrs['LP_words']
        self.words_per_mf = missionAttrs['words_per_mf']
        self.num_of_sector_counts = missionAttrs['num_of_sector_counts']
        self.g_nNomRecLen = missionAttrs['g_nNomRecLen']
        self.nNomDataLen = missionAttrs['nNomDataLen']
        self.NumOfInstr = missionAttrs['NumOfInstr']
        self.epoch_fillVal = missionAttrs['epoch_fillVal']
        self.mfPerFullSweep = missionAttrs['mfPerFullSweep']
        self.Instr_Energy = missionAttrs['Instr_Energy']
        self.geometric_factor = missionAttrs['geometric_factor']
        self.deadtime = missionAttrs['deadtime']
        self.Instr_sector_to_pitch = missionAttrs['Instr_sector_to_pitch']
        self.Instr_words_per_sweep = missionAttrs['Instr_words_per_sweep']
        self.Instr_Energy = missionAttrs['Instr_Energy']
        self.esaMaxCounts = missionAttrs['esaMaxCounts']
        self.LP_Variables = missionAttrs['LP_Variables']
        self.MinorFrameTime = missionAttrs['MinorFrameTime']
        self.LPSamplePeriod = missionAttrs['LPSamplePeriod']
        self.Launch_Times = missionAttrs['Launch_Times']
        self.Andoya_magnetic_inclination = missionAttrs['Andoya_magnetic_inclination']
        self.LPFixed_calResistances = missionAttrs['LPFixed_calResistances']
        self.LPswept_voltage_range = missionAttrs['LPswept_voltage_range']
        self.LPFixedProbeBias = missionAttrs['LPFixedProbeBias']
        self.Epoch_range_to_determine_stepDAC = missionAttrs['Epoch_range_to_determine_stepDAC']
        self.LP_probe_areas = missionAttrs['LP_probe_areas']
        self.startEndLangmuirBreakIntoCurves = missionAttrs['startEndLangmuirBreakIntoCurves']
        self.Andoya_Space_Lat_Long = missionAttrs['Andoya_Space_Lat_Long']
        self.LPswept_cal_epoch_ranges = missionAttrs['LPswept_cal_epoch_ranges']
        self.LP_probe_areas = missionAttrs['LP_probe_areas']
        self.LPswept_cal_resistances = missionAttrs['LPswept_cal_resistances']
        self.LPswept_cal_epoch_ranges_single_sweep = missionAttrs['LPswept_cal_epoch_ranges_single_sweep']
        self.timeBetweenSteps_in_ns = missionAttrs['timeBetweenSteps_in_ns']
        self.ringCoreScaleFactors = missionAttrs['ringCoreScaleFactors']
        self.ringCore5thOrderCorrections = missionAttrs['ringCore5thOrderCorrections']
        self.ringCoreCalMatrix = missionAttrs['ringCoreCalMatrix']
        self.AvgConingRate = missionAttrs['AvgConingRate']
        self.AvgRollRate = missionAttrs['AvgRollRate']



# --- --- --- --- --- ---
# --- ACES II Mission ---
# --- --- --- --- --- ---
def ACES_mission_dicts():

    from copy import deepcopy

    ACESII_attrs_dict = {
        'Andoya_Space_Lat_Long':[69.294167, 16.020833],
        'Andoya_magnetic_inclination':78.1300,
        'missionNam': 'ACESII',
        'rocketID': ['36359', '36364'],
        'globalAttributes': [
            {'Source_name': f'ACESII_36359>ACES II 36.359',
             'Data_type': 'K0>Key Parameter',
             'PI_name': 'Bounds',
             'Logical_source': f'aces_36359_',
             'Logical_file_id': f'aces_36359_00000000_v01',
             'Logical_source_description': 'Raw Data from the ACESII mission organized by minorframe.150 words per minor frame.40 minor frames to a major frame.',
             'TEXT': 'Raw Data from the ACESII mission organized by minorframe.150 words per minor frame.40 minor frames to a major frame.'
             }, {
                'Source_name': f'ACESII_36364>ACES II 36.364',
                'Data_type': 'K0>Key Parameter',
                'PI_name': 'Bounds',
                'Logical_source': f'aces_36364_',
                'Logical_file_id': f'aces_36364_00000000_v01',
                'Logical_source_description': 'Raw Data from the ACESII mission organized by minorframe.150 words per minor frame.40 minor frames to a major frame.',
                'TEXT': 'Raw Data from the ACESII mission organized by minorframe.150 words per minor frame.40 minor frames to a major frame.'
            }],
        'Launch_Times': [722236869184000000, 722236969184000000], # TT2000 values corresponding to 17:20:00 and 17:21:40 for high/low flyer, respectively.
        'AvgConingRate': [0.05397837461866082, 0.11071336310292652],
        'AvgRollRate': [0.671, 0.547],
        'g_nWords': 150,
        'g_nNomRecLen': 150 * 2 + 12,
        'nNomDataLen': 150 - 2,
        'gSync': b'\x40\x28\x6b\xfe',
        # g_sSync = "\xfe\x6b\x28\x40"  is the real frame sync but the data is little endian
        'nDataLen': "74I",
        'InstrNames': ['EEPAA', 'LEESA', 'IEPAA', 'LP'],
        'InstrNames_LC': ['eepaa', 'leesa', 'iepaa', 'lp'],
        'NumOfInstr': 4,
        'mfPerFullSweep': 4,
        'num_of_sector_counts': [21, 21, 7, 8],
        'words_per_mf': [9, 9, 7, 9],
        'ESA_words': [[10, 26, 43, 60, 85, 101, 114, 126, 142],
                      [14, 28, 45, 61, 89, 103, 116, 132, 144],
                      [11, 27, 44, 86, 102, 115, 143]],
        'epoch_fillVal': -9223372036854775808,
        'esaMaxCounts': 4095,
        'timeBetweenSteps_in_ns':0.001 * (10**(9)),
        'ESA_CLK_INPUT': 625,
        'LP_words': [15, 29, 46, 66, 90, 104, 117, 133, 145],
        'LP_Variables': ["deltaNdivN","step", "ne_swept", "ni_swept","ni"],
        'MinorFrameTime': 250000,
        'Count_Interval': 917, # measured in ns
        'LPSamplePeriod': 31250,
        'LPFixedProbeBias':[-5.05,-4.96],
        'LPFixed_calResistances':[ {'Open':165,500000000:2170,250000000:2308,100000000:2540,50000000:2710,10000000:3122,5000000:3299},
                                   {'Open':165,500000000:2063,250000000:2260,100000000:2500,50000000:2660,10000000:3060,5000000:3194}],
        # 'LPFixed_calResistances':[ {100000000:2540,50000000:2710,10000000:3122,5000000:3299},
        #                            {100000000:2500,50000000:2660,10000000:3060,5000000:3194}],
        'LPswept_voltage_range':[[-4.72, 2.12], [-4.68, 2.08]],
        'LPswept_cal_resistances':[10*10**(9),500*10**(6),250*10**(6),100*10**(6),50*10**(6),10*10**(6),5*10**(6)],
        'LPswept_cal_epoch_ranges_single_sweep':[
                [ # HIGH FLYER
             [dt.datetime(2022, 11, 3, 13, 56, 11, 910063), dt.datetime(2022, 11, 3, 13, 56, 14, 88450)], #open
             [dt.datetime(2022, 11, 3, 13, 53, 30, 000000), dt.datetime(2022, 11, 3, 13, 53, 32, 141247)], #500M
             [dt.datetime(2022, 11, 3, 13, 54, 9, 973502), dt.datetime(2022, 11, 3, 13, 54, 12, 198665)],# 250M
             [dt.datetime(2022, 11, 3, 13, 54, 42, 500000), dt.datetime(2022, 11, 3, 13, 54, 44, 000000)],# 100M
             [dt.datetime(2022, 11, 3, 13, 55, 1, 984859), dt.datetime(2022, 11, 3, 13, 55, 4, 164512)],# 50M
             [dt.datetime(2022, 11, 3, 13, 55, 31, 964699), dt.datetime(2022, 11, 3, 13, 55, 34, 207151)],# 10M
             [dt.datetime(2022, 11, 3, 13, 56, 3, 893412), dt.datetime(2022, 11, 3, 13, 56, 6, 218441)] # 5M
                ],
                [
             [dt.datetime(2022, 11, 3, 14, 25, 24, 823952), dt.datetime(2022, 11, 3, 14, 25, 27, 22060)],# open
             [dt.datetime(2022, 11, 3, 14, 23, 4, 846869), dt.datetime(2022, 11, 3, 14, 23, 7, 845)],# 500M
             [dt.datetime(2022, 11, 3, 14, 23, 56, 839534), dt.datetime(2022, 11, 3, 14, 23, 59, 1309)],# 250M
             [dt.datetime(2022, 11, 3, 14, 23, 56, 825072), dt.datetime(2022, 11, 3, 14, 23, 59, 18086)],# 100M
             [dt.datetime(2022, 11, 3, 14, 24, 22, 871113), dt.datetime(2022, 11, 3, 14, 24, 24, 994890)],# 50M
             [dt.datetime(2022, 11, 3, 14, 24, 54, 864118), dt.datetime(2022, 11, 3, 14, 24, 56, 974901)],# 10M
             [dt.datetime(2022, 11, 3, 14, 25, 16, 860398), dt.datetime(2022, 11, 3, 14, 25, 18, 981178)],# 5M

                ]

             ],
        'LPswept_cal_epoch_ranges': # open 500M 250M 100M 50M 10M 5M
            [
                [ # HIGH FLYER
             [dt.datetime(2022, 11, 3, 13, 56, 11, 264921), dt.datetime(2022, 11, 3, 13, 56, 13, 20998)], #open
             [dt.datetime(2022, 11, 3, 13, 53, 28, 500000), dt.datetime(2022, 11, 3, 13, 53, 32, 200000)], #500M
             [dt.datetime(2022, 11, 3, 13, 54, 8, 000000), dt.datetime(2022, 11, 3, 13, 54, 16, 000000)],# 250M
             [dt.datetime(2022, 11, 3, 13, 54, 38, 000000), dt.datetime(2022, 11, 3, 13, 54, 46, 000000)],# 100M
             [dt.datetime(2022, 11, 3, 13, 55, 2, 000000), dt.datetime(2022, 11, 3, 13, 55, 16, 000000)],# 50M
             [dt.datetime(2022, 11, 3, 13, 55, 38, 000000), dt.datetime(2022, 11, 3, 13, 55, 42, 000000)],# 10M
             [dt.datetime(2022, 11, 3, 13, 56, 0, 000000), dt.datetime(2022, 11, 3, 13, 56, 8, 000000)] # 5M
                ],
                [
             [dt.datetime(2022, 11, 3, 14, 25, 23, 000000), dt.datetime(2022, 11, 3, 14, 25, 26, 500000)],# open
             [dt.datetime(2022, 11, 3, 14, 23, 5, 000000), dt.datetime(2022, 11, 3, 14, 23, 9, 000000)],# 500M
             [dt.datetime(2022, 11, 3, 14, 23, 33, 445653), dt.datetime(2022, 11, 3, 14, 23, 41, 321627)],# 250M
             [dt.datetime(2022, 11, 3, 14, 23, 57, 127446), dt.datetime(2022, 11, 3, 14, 24, 4, 669378)],# 100M
             [dt.datetime(2022, 11, 3, 14, 24, 22, 943856), dt.datetime(2022, 11, 3, 14, 24, 28, 845243)],# 50M
             [dt.datetime(2022, 11, 3, 14, 24, 46, 753576), dt.datetime(2022, 11, 3, 14, 24, 57, 182238)],# 10M
             [dt.datetime(2022, 11, 3, 14, 25, 11, 122669), dt.datetime(2022, 11, 3, 14, 25, 19, 107600)],# 5M
                ]

             ],
        'LP_probe_areas': [[0.002014, 0.002014], [0.002014, 0.002014]], # square meters
        'Epoch_range_to_determine_stepDAC': [[dt.datetime(2022, 11, 20,17,24,30,300000),dt.datetime(2022, 11, 20,17,25,10,310000)],
                                             [dt.datetime(2022, 11, 20,17,24,30,440000),dt.datetime(2022, 11, 20,17,25,10,450000)]],
        'startEndLangmuirBreakIntoCurves':[[dt.datetime(2022, 11, 20,17,21,00,890000),dt.datetime(2022, 11, 20,17,29,57,700000) ],
                                           [dt.datetime(2022, 11, 20,17,23,5,  10000),dt.datetime(2022, 11, 20,17,28,9,900000)]],
        'InstrNames_Full': ['EEPAA>Electron Energy Pitch Angle Analyzer',
                            'LEESA>Low Energy Electrostatic Analyzer',
                            'IEPAA>Ion Energy Pitch Angle Analyzer',
                            'Langmuir_Probe>Langmuir Probe'],
        # 'geometric_factor': [
        #     [0.000174 for i in range(21)],
        #     [0.000174/100 for i in range(21)], # LEESA geofactor was ~EEPAA/100
        #     [0.000174 for i in range(7)]],
        'geometric_factor': [
            [8.63E-5 for i in range(21)], # CONFIRMED: in units of cm^2 str^1
            [8.63E-5 / 100 for i in range(21)],  # LEESA geofactor was ~EEPAA/100
            [8.63E-5 for i in range(7)]],
        # 'deadtime': [674E-9, 674E-9],
        'deadtime': [324E-9, 324E-9],
        'Instr_sector_to_pitch': [
            [-10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190],
            [-10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190],
            [180, 150, 120, 90, 60, 30, 0]
        ],
        'Instr_words_per_sweep': [49, 49, 49],
        'Instr_Energy': [ #eepaa
            np.array(
                [13678.4, 11719.21, 10040.64, 8602.5, 7370.34, 6314.67, 5410.2, 4635.29,
                 3971.37, 3402.54, 2915.18, 2497.64, 2139.89, 1833.39, 1570.79, 1345.8,
                 1153.04, 987.89, 846.39, 725.16, 621.29, 532.3, 456.06, 390.74,
                 334.77, 286.82, 245.74, 210.54, 180.39, 154.55, 132.41, 113.45,
                 97.2, 83.28, 71.35, 61.13, 52.37, 44.87, 38.44, 32.94,
                 28.22, 24.18, 20.71, 17.75, 15.21, 13.03, 11.16, 9.56, 8.19]),
                        # leesa
            np.array(
                [13678.4, 11719.21, 10040.64, 8602.5, 7370.34, 6314.67, 5410.2, 4635.29,
                 3971.37, 3402.54, 2915.18, 2497.64, 2139.89, 1833.39, 1570.79, 1345.8,
                 1153.04, 987.89, 846.39, 725.16, 621.29, 532.3, 456.06, 390.74,
                 334.77, 286.82, 245.74, 210.54, 180.39, 154.55, 132.41, 113.45,
                 97.2, 83.28, 71.35, 61.13, 52.37, 44.87, 38.44, 32.94,
                 28.22, 24.18, 20.71, 17.75, 15.21, 13.03, 11.16, 9.56, 8.19])/1000,
                        # iepaa
            np.array(
                [13678.4, 11719.21, 10040.64, 8602.5, 7370.34, 6314.67, 5410.2, 4635.29,
                 3971.37, 3402.54, 2915.18, 2497.64, 2139.89, 1833.39, 1570.79, 1345.8,
                 1153.04, 987.89, 846.39, 725.16, 621.29, 532.3, 456.06, 390.74,
                 334.77, 286.82, 245.74, 210.54, 180.39, 154.55, 132.41, 113.45,
                 97.2, 83.28, 71.35, 61.13, 52.37, 44.87, 38.44, 32.94,
                 28.22, 24.18, 20.71, 17.75, 15.21, 13.03, 11.16, 9.56, 8.19])
                        ],
        'ringCoreScaleFactors': [[-0.00759273866010453, -0.00983337688949566, -0.00910656152667623],
                                 [0.00887255961330481, 0.00857834057158005, 0.00795888490103171]],
        'ringCore5thOrderCorrections': [
            # high flyer
            [[-2.01426584458453e-22, 4.78339584524288e-18, 5.14890545473976e-12, - 2.35923633686087e-08, 0.989892268554107, 851.689485242279],
             [-3.94182035267970e-21, 3.78429607076847e-18, 3.25132527020627e-11, - 2.51430072058701e-08, 0.938231980320960, 521.911450224810],
             [-6.08971576025967e-21, 1.48213776699574e-17, 4.75798877553707e-11, - 5.78468752225862e-08, 0.909294860266228, - 40.0287112936637]],
            # low flyer
            [[1.66474802438242e-22, 1.69219967182715e-18, -2.80634927683823e-13, -3.70529011377136e-09, 1.00277320989464, 974.278755837680],
             [-3.47456898510402e-21, 4.22701226361097e-18, 3.01652581051724e-11, -5.86821891732926e-09, 0.943213559983693, 610.204104911557],
             [-2.90175291682798e-22, 3.73789350207806e-18, 5.48134539223753e-12, -5.53834631668198e-09, 0.985171491538179, 69.4056944418974]]
        ],
        'ringCoreCalMatrix':[
            # high flyer - APPLY TO ROCKET FRAME
            [[0.980811647255180,	-0.00505723015508628,	-0.0188341290296427, -704.129906360017],
             [0.00847115917661956,	 0.973814621560943,      0.0751864354553729,  811.489357644093],
             [0.0233608047714457,	-0.0792669223938717,	 0.943802171760235,	  59.4473792734371]],
            # low flyer - APPLY TO ROCKET FRAME
            [[0.980811647255180,	-0.00505723015508628,	-0.0188341290296427,	-704.129906360017],
             [0.00847115917661956,	0.973814621560943,	0.0751864354553729,	811.489357644093],
             [0.0233608047714457,	-0.0792669223938717,	0.943802171760235,	59.4473792734371]]
        ]
    }

    ACESII_attrs = makeRocketAttrs(ACESII_attrs_dict)


    # Initialization Information for L0 Variables
    EEPAA_data_dict = {
        'Epoch':                    [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': ACESII_attrs.epoch_fillVal,'FORMAT': 'I5','UNITS': 'ns','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','MONOTON':'INCREASE','TIME_BASE':'J2000','TIME_SCALE':'Terrestrial Time','REFERENCE_POSITION':'Rotating Earth Geoid','SCALETYP':'linear'}],
        'STEPPER_Voltage_ADC':      [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'MCP_Voltage_ADC':          [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'MCP_Current_ADC':          [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'STACK_Voltage_ADC':        [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'STACK_Current_ADC':        [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'minor_frame_counter':      [[],                                             {'DEPEND_0':  None,  'DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'Count_Interval':           [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'major_frame_counter':      [[],                                             {'DEPEND_0':  None,  'DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'Sector_Counts':            [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': 'Sector_Number','DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'data','SCALETYP':'linear'}],
        'Sector_Number':            [[[i + 1 for i in (range(21))]],                 {'DEPEND_0':  None,  'DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'Pitch_Angle':              [[[(-10 + i*10) for i in (range(21))]],          {'DEPEND_0':  None,  'DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': 'deg','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear','LABLAXIS': 'Pitch_Angle'}],
        '625kHz_Clock_Input':       [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        '3p3V_Voltage_Monitor_ADC': [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        '5V_Voltage_Monitor_ADC':   [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'Temperature_Monitor_ADC':  [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'sfid':                     [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'EXP_Current':              [[],                                             {'DEPEND_0': 'Epoch','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        '28V_Monitor':              [[],                                             {'DEPEND_0': 'Epoch_monitors','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'Boom_Monitor':             [[],                                             {'DEPEND_0': 'Epoch_monitors','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': -1,'FORMAT': 'I5','UNITS': '#','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','SCALETYP':'linear'}],
        'Epoch_monitors':           [[],                                             {'DEPEND_0': 'Epoch_monitors','DEPEND_1': None,'DEPEND_2': None,'FILLVAL': ACESII_attrs.epoch_fillVal,'FORMAT': 'I5','UNITS': 'ns','VALIDMIN':None,'VALIDMAX': None,'VAR_TYPE':'support_data','MONOTON':'INCREASE','TIME_BASE':'J2000','TIME_SCALE':'Terrestrial Time','REFERENCE_POSITION':'Rotating Earth Geoid','SCALETYP':'linear'}],
        'sync_word':                [[],                                             {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -1, 'FORMAT': 'I5', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}],
        'status_word_1':            [[],                                             {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -1, 'FORMAT': 'I5', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}],
        'status_word_3':            [[],                                             {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -1, 'FORMAT': 'I5', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}],
        'HV_div16':                 [[],                                             {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -1, 'FORMAT': 'I5', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}],
        'HV_enable':                [[],                                             {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -1, 'FORMAT': 'I5', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}],
        'sweep_step':               [[],                                             {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -1, 'FORMAT': 'I5', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}],
        'TP5_enable':               [[],                                             {'DEPEND_0': 'Epoch', 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': -1, 'FORMAT': 'I5', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]
    }
    for key, val in EEPAA_data_dict.items():
        EEPAA_data_dict[key][1]['LABLAXIS'] = key
        EEPAA_data_dict[key][1]['FIELDNAM'] = key
        EEPAA_data_dict[key][1]['CATDESC'] = key

    # -- LEESA ---
    LEESA_data_dict = deepcopy(EEPAA_data_dict)

    # -- IEPAA ---
    IEPAA_data_dict = deepcopy(EEPAA_data_dict)
    del IEPAA_data_dict['Sector_Counts'], IEPAA_data_dict['Sector_Number'], IEPAA_data_dict['Pitch_Angle']
    IEPAA_data_dict['Sector_Counts'] = [[], EEPAA_data_dict['Sector_Counts'][1]]  # IEPAA has 7 anodes with 30deg/anode totally 210deg coverage
    IEPAA_data_dict['Sector_Number'] = [[j + 1 for j in (range(ACESII_attrs.words_per_mf[2]))], EEPAA_data_dict['Sector_Number'][1]]  # IEPAA has 7 anodes with 30deg/anode totally 210deg coverage
    IEPAA_data_dict['Pitch_Angle'] = [[j*30 for j in range(ACESII_attrs.words_per_mf[2])], EEPAA_data_dict['Pitch_Angle'][1]]

    # -- LPs ---
    LP_data_dict = deepcopy(EEPAA_data_dict)
    LPremove = ['STEPPER_Voltage_ADC',    'MCP_Voltage_ADC',         'MCP_Current_ADC',
                'STACK_Voltage_ADC',      'STACK_Current_ADC',       'Sector_Counts',
                'Sector_Number',          'Pitch_Angle',             '3p3V_Voltage_Monitor_ADC',
                '5V_Voltage_Monitor_ADC', 'Temperature_Monitor_ADC',
                'Boom_Monitor',           'Count_Interval',          '625kHz_Clock_Input',
                'sync_word',              'status_word_1',
                'TP5_enable',             'status_word_3',           'HV_div16',
                'HV_enable',              'sweep_step']

    for thing in LPremove:
        del LP_data_dict[thing]

    LP_data_dict['Channel_Number'] =   [[[i for i in (range(len(ACESII_attrs_dict['LP_words'])))]], {'CATDESC': 'Channel_Number', 'DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FIELDNAM': 'Channel_Number', 'FILLVAL': -1, 'FORMAT': 'I5', 'LABLAXIS': 'Channel_Number', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]
    LP_data_dict['Channel_Counts'] =   [[],                            {'CATDESC': 'Channel_Counts', 'DEPEND_0': 'Epoch', 'DEPEND_1': 'Channel_Number', 'DEPEND_2': None, 'FIELDNAM': 'Channel_Counts', 'FILLVAL': -1, 'FORMAT': 'I5','LABLAXIS': 'Channel_Counts', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None,'VAR_TYPE': 'data', 'SCALETYP': 'linear'}]
    LP_data_dict['Boom_Monitor_1'] =   [[],                            {'CATDESC': 'Boom_Monitor_1', 'DEPEND_0': 'Epoch_monitor_1', 'DEPEND_1': None, 'DEPEND_2': None, 'FIELDNAM': 'Boom_Monitor_1', 'FILLVAL': -1, 'FORMAT': 'I5', 'LABLAXIS': 'Boom_Monitor_1', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]
    LP_data_dict['Boom_Monitor_2'] =   [[],                            {'CATDESC': 'Boom_Monitor_2', 'DEPEND_0': 'Epoch_monitor_2', 'DEPEND_1': None, 'DEPEND_2': None, 'FIELDNAM': 'Boom_Monitor_2', 'FILLVAL': -1, 'FORMAT': 'I5', 'LABLAXIS': 'Boom_Monitor_2', 'UNITS': '#', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'SCALETYP': 'linear'}]
    LP_data_dict['Epoch_monitor_1'] =  [[],                            {'CATDESC': 'Epoch_monitor_1','DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': ACESII_attrs.epoch_fillVal, 'FORMAT': 'I5', 'UNITS': 'ns', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000', 'TIME_SCALE': 'Terrestrial Time', 'REFERENCE_POSITION': 'Rotating Earth Geoid', 'SCALETYP': 'linear'}]
    LP_data_dict['Epoch_monitor_2'] =  [[],                            {'CATDESC': 'Epoch_monitor_2','DEPEND_0': None, 'DEPEND_1': None, 'DEPEND_2': None, 'FILLVAL': ACESII_attrs.epoch_fillVal, 'FORMAT': 'I5', 'UNITS': 'ns', 'VALIDMIN': None, 'VALIDMAX': None, 'VAR_TYPE': 'support_data', 'MONOTON': 'INCREASE', 'TIME_BASE': 'J2000', 'TIME_SCALE': 'Terrestrial Time', 'REFERENCE_POSITION': 'Rotating Earth Geoid', 'SCALETYP': 'linear'}]

    data_dicts = [EEPAA_data_dict, LEESA_data_dict, IEPAA_data_dict,LP_data_dict,LP_data_dict]



    # --- Deconvolution Keys for the ESAs ---
    deConvolveKey_EEPAA = {
        1: None,
        2: 'sync_word',
        3: 'status_word_1',
        4: 'status_word_2',
        5: 'status_word_3',
        6: 'STEPPER_Voltage_ADC',
        7: 'MCP_Voltage_ADC',
        8: 'MCP_Current_ADC',
        9: 'STACK_Voltage_ADC',
        10: 'STACK_Current_ADC',
        11: '3p3V_Voltage_Monitor_ADC',
        12: '5V_Voltage_Monitor_ADC',
        13: 'Temperature_Monitor_ADC',
        14: 'Sector_Counts',
        15: 'Sector_Counts',
        16: 'Sector_Counts',
        17: 'Sector_Counts',
        18: 'Sector_Counts',
        19: 'Sector_Counts',
        20: 'Sector_Counts',
        21: 'Sector_Counts',
        22: 'Sector_Counts',
        23: 'Sector_Counts',
        24: 'Sector_Counts',
        25: 'Sector_Counts',
        26: 'Sector_Counts',
        27: 'Sector_Counts',
        28: 'Sector_Counts',
        29: 'Sector_Counts',
        30: 'Sector_Counts',
        31: 'Sector_Counts',
        32: 'Sector_Counts',
        33: 'Sector_Counts',
        34: 'Sector_Counts',
        35: '625kHz_Clock_Input',
        36: None
    }
    deConvolveKey_LEESA = deepcopy(deConvolveKey_EEPAA)
    deConvolveKey_IEPAA = deepcopy(deConvolveKey_EEPAA)
    deConvoleKeys = [deConvolveKey_EEPAA,deConvolveKey_LEESA,deConvolveKey_IEPAA]
    for i in [21 + i for i in (range(16))]:
        del deConvolveKey_IEPAA[i]
    deConvolveKey_IEPAA[21] = '625kHz_Clock_Input'
    for i in [22,23,24,25,26,27,28]:
        deConvolveKey_IEPAA[i] = None

    return ACESII_attrs,data_dicts,deConvoleKeys



# --- --- --- --- --- ---
# --- TRICEII Mission ---
# --- --- --- --- --- ---
def TRICE_mission_dicts():
    TRICE_attrs_dict = {
        'missionNam' : 'TRICE',
        'rocketID': ['52003', '52004'],
        'g_nWords': 150,
        'g_nNomRecLen': 150 * 2 + 12,
        'nNomDataLen': 150 - 2,
        'gSync' : b'\x40\x28\x6b\xfe',
        'nDataLen': "59I",
        'InstrNames_LC': ['eepaa']
    }

    TRICEII_attrs = makeRocketAttrs(TRICE_attrs_dict)

    return TRICEII_attrs



