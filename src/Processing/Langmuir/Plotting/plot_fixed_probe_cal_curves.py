#--- plot_fixed_probe_cal_curves.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: plot the calibration curves used for the LP fixed probes


if plotFixedCalCurve:
    figWidth = 10
    figHeight = 10
    Label_Fontsize = 15
    Title_Fontsize = 15

    import matplotlib

    matplotlib.rc('figure', figsize=(5, 5))
    plt.rcParams["figure.figsize"] = (figHeight, figWidth)
    xDataFit = np.array([i for i in range(1, 4096)])
    yDataFit = [calFunction_fixed(val, *parameters) for val in xDataFit]
    plt.plot(xDataFit, yDataFit, color='red')
    plt.scatter(analog_vals, calibrationCurrents)
    plt.xlabel('ADC Value', fontsize=Label_Fontsize)
    plt.ylabel(r'Ln($I_{cal}$) [nA]', fontsize=Label_Fontsize)
    plt.suptitle(f'FIXED LP - {rocketAttrs.rocketID[wRocket - 4]}\n'
                 'Calculated calibration Current vs Analog Value', fontsize=Title_Fontsize)
    plt.legend(['ln(y) = mx + b\n'
                f'm: {parameters[0]}\n'
                f'b: {parameters[1]}'])
    plt.savefig(
        rf'C:\Users\cfelt\Desktop\rockets\ACES-II\Instruments\fixedLangmuirProbe\Calibration\ACESII_{fliers[wflyer]}_CalCurve.png')