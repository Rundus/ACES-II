# # --- magCal_investigation.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: script for general tinkering

import glob
from ACESII_code.myImports import *

directory = r'C:\Users\cfeltman\PycharmProjects\UIOWA_CDF_operator\Data\Integration\tmCDF\high'
files = glob.glob(rf'{directory}\*.cdf')
usefilter = True

for i in [i for i in (range(len(files))) if i not in [5]]:
    file = files[i]

    # print(f'Collecting Data for {files[i].replace(directory,"")}')
    with pycdf.CDF(file) as cdfFile:
        sfidCDF = cdfFile['sfid'][...]
        # epochCDF= cdfFile['Epoch'][...]
        # mfCDF = cdfFile['minorframe'][...]

    startpt,endpt = np.where(sfidCDF == 0),np.where(sfidCDF == 39)
    start,end = startpt[0][1], endpt[0][-1] + 1
    data_sfid = sfidCDF[start:end]

    # data_epoch= epochCDF[start:end]
    # data_mf = mfCDF[start:end]

    counter = data_sfid[0]
    corruptData = []
    storefirst = []

    if usefilter:
        for index in range(len(data_sfid)):
            if data_sfid[index] != counter:

                if len(storefirst) == 0:
                    storefirst.append([i,index,data_sfid[index]])

                if data_sfid[index] >= 40:
                    data_sfid[index] = counter
                    # data_epoch[index] = rocketAttributes.ACESattrs.epoch_fillVal
                    # data_mf[index] = [-1 for i in (range(rocketAttributes.ACESattrs.g_nWords))]
                else:
                    data_sfid = np.insert(data_sfid, index, counter, axis=0)
                    # data_epoch=np.insert(data_epoch,index,rocketAttributes.ACESattrs.epoch_fillVal,axis=0)
                    # data_mf=np.insert(data_mf,index,[-1 for i in (range(rocketAttributes.ACESattrs.g_nWords))],axis=0)
            counter += 1
            if counter == 40:
                counter = 0

    print(storefirst)



    # counter = data_sfid[0]
    # for index, num in enumerate(tqdm(data_sfid)):
    #     if num != counter:
    #         corruptData.append([index, num, counter])
    #
    #     counter += 1
    #     if counter >= 40:
    #         counter = 0
    #
    # if corruptData == []:
    #     print(f'THERE IS NO CORRUPT DATA IN [{i}]')
    # else:
    #     print(f'THERE IS CORRUPT DATA IN [{i}]: {100*len(corruptData)/len(data_sfid)}%', len(corruptData))









corruptions = [[0, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_All_Fire_Sequence_Test_202208223.cdf'],
             [1, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_Forward_Deployments_20220829_1.cdf'],
             [2, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_GPS_Rollout_20220901_1.cdf'],
             [3, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_LEESA_Deploy_20220827_1.cdf'],
             [4, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_LEESA_Deploy_20220827_2.cdf'],
             [5, '38.5173031252387%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_Mag_Cal_20220906_1.cdf'],
             [6, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_Mag_check_20220825_3.cdf'],
             [7, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_NoseCone_Deploy_20220826_1.cdf'],
             [8, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_Operational_Spin_20220826_1.cdf'],
             [9, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_Post_Vibe_Sequence_20220831_1.cdf'],
             [10, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_Post_Vibe_Sequence_20220831_2.cdf'],
             [11, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_Post_Vibe_Sequence_20220831_3.cdf'],
             [12, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_Sequence_Test2_AllFire_202208224.cdf'],
             [13, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_Sequence_Test2_AllFire_202208224_1.cdf'],
             [14, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_Sequence_Test3_NoFire_202208224_2.cdf'],
             [15, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_Sequence_Test4_PowerBackup2_202208224_3.cdf'],
             [16, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_Vibe_RandomX_20220825_2.cdf'],
             [17, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_Vibe_RandomY_20220825_3.cdf'],
             [18, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_Vibe_RandomZ_20220825_1.cdf'],
             [19, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_359_Vibe_Sine_Thrust_20220825_2.cdf'],
             [20, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_Aft_Boom_Deplay_20220830_1.cdf'],
             [21, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_Forward_Deploy_20220831_1.cdf'],
             [22, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_GPS_RollOut_20220902_1.cdf'],
             [23, '53.21901896512871%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_Mag_Cal_20220906_1.cdf'],
             [25, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_NoseCone_Deploy_20220830_1.cdf'],
             [26, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_Operational_Spin_20220830_1.cdf'],
             [27, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_Post_Vibe_Sequence_20220902_1.cdf'],
             [28, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_Post_Vibe_Sequence_20220902_2.cdf'],
             [29, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_Post_Vibe_Sequence_20220902_3.cdf'],
             [30, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_Post_Vibe_Sequence_20220902_4.cdf'],
             [31, '68.2742960401532%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_Sequence_Test_Allfire_20220825_1.cdf'],
             [32, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_Sequence_Test_Nofire_20220825_2.cdf'],
             [33, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_Sequence_Test_PowerBackup2_20220825_3.cdf'],
             [34, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_Vibe_Random_Thrust_20220829_2.cdf'],
             [35, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_Vibe_Random_X_20220829_3.cdf'],
             [36, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_Vibe_Random_Y_20220829_4.cdf'],
             [37, '9.853066151022871e-05%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\36_364_Vibe_Sine_Thrust_20220829_1.cdf'],
             [38, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\IEPAA_Highflyer_Integration_Test_8172022_1.cdf'],
             [39, '0% (After Filter)', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\IEPAA_Highflyer_Integration_Test_8172022_2.cdf'],
             [40, '0% (After Filter)', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\IEPAA_Highflyer_Integration_Test_8172022_3.cdf'],
             [41, '0% (After Filter)', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\LEESA_Highflyer_Integration_Test_8172022_1.cdf'],
             [42, '0% (After Filter)', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\LP_Highflyer_Integration_Test_8172022_1.cdf'],
             [43, '0% (After Filter)', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\LP_Lowflyer_Integration_Test_8182022_1.cdf'],
             [44, '0% (After Filter)', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\MPI_LowFlyer_Integration_Test_8192022_2.cdf'],
             [45, '0.0%', 'C:\\Users\\cfeltman\\PycharmProjects\\UIOWA_CDF_operator\\Data\\Integration\\tmCDF\\MPI_LowFlyer_Integration_Test_8192022_3.cdf']]
