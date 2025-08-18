

class DespinToggles:


    # --- --- --- Reduce Data --- --- ---
    import datetime as dt
    reduceTimes = [[dt.datetime(2022, 11, 20, 17, 22, 00, 000000), dt.datetime(2022, 11, 20, 17, 28, 00, 000000)],
                   [dt.datetime(2022, 11, 20, 17, 22, 00, 000000), dt.datetime(2022, 11, 20, 17, 28, 00, 000000)]]


    # --- --- --- Apply Kenton/Antonio T0 correction --- --- ---
    # Description: Antontio worked with kenton to determine how to best
    # determine the T0 for the internal integration_tad_files timer and the rocket T0. This matters
    # When trying to determine the deltat between E-Field and B-Field measurements.
    # For similar events observed in the E/B fields, two linear fits were stl.done to try to bring the epochs into alignment
    # based on the peaks of the two waveforms:
    KentonAntonio_T0_Correction = True
    EB_East = [0.999994, 0.03374]  # slope, intercept (in seconds)
    EB_North = [0.999986, 0.03294]
    # The fits are stl.done by finding 5 deltaT



    # --- Fit Results ---
    fitResults = {
        'Bx': {'Spin Amp': 25.42873940404161, 'Spin Freq': 0.6463295881639182, 'Spin Phase': 91.9759995936283,
               'Cone Amp': 625.8772357084948, 'Cone Freq': 0.05294818121871208, 'Cone Phase': -138.77308595997619,
               'Offset': -44919.748937299344},
        'By': {'Spin Amp': 7.378420193701481, 'Spin Freq': 0.6442248190622027, 'Spin Phase': 109.20255873087793,
               'Cone Amp': 1380.5616077430786, 'Cone Freq': 0.02700105226961604, 'Cone Phase': 109.87799606103452,
               'Offset': -139.74554466082876},
        'Bz': {'Spin Amp': 8.095746809541962, 'Spin Freq': 0.6442537451458561, 'Spin Phase': 19.11852573798773,
               'Cone Amp': 1257.0313161879794, 'Cone Freq': 0.026874206798816504, 'Cone Phase': -69.78175516947503,
               'Offset': 32.456720919269245}
    }

    coneFreq = sum([fitResults['By']['Cone Freq'], fitResults['Bz']['Cone Freq']]) / 2
    spinFreq = sum([fitResults['Bz']['Spin Freq'], fitResults['By']['Spin Freq'],
                    fitResults['Bz']['Spin Freq']]) / 3 if wRocket == 4 else 0.55