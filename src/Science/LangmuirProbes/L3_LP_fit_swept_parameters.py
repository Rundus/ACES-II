# --- L3_LP_fit_swept_parameters.py ---
# THIS CODE WAS WRITTEN BY J-HOMANN

# Imports
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from src.my_imports import *

# Picking High vs. Low
var = 0 # pick which rocket to select
rocket = ['high', 'low']
Id = ['36359', '36364']

# Setting the minimum data value to begin fitting
x_data_min = -1

# Sweep Toggles
sweep_min = 1
sweep_max = [533, 312]  # 313 for Low, 533 for High
sweep_step = 2

# Physical constants
e = 1.602e-19
k = 1.381e-23
m_e = 9.11e-31
A = 0.002014

# Parameter guesses
delta_phi_guess = 0.5

# Parameter tolerance
n_tol = 3
T_tol = 0.05 * (e / k)
delta_phi_tol = 0.75

# Getting the Swept Langmuir Probe Data
file = rf'{DataPaths.ACES_data_folder}/science/payload_potential/{rocket[var]}/ACESII_{Id[var]}_FloatingPotential.cdf'
data_dict = stl.loadDictFromFile(file)
I_sweeps = deepcopy(data_dict['sweeps_current'][0])
V_sweeps = deepcopy(data_dict['sweeps_voltage'][0])
V_floating = deepcopy(data_dict['floating_potential'][0])
Epoch = deepcopy(data_dict['Epoch'][0])

# Getting the Density Data
density_file = rf"{DataPaths.ACES_data_folder}/L3/Langmuir/{rocket[var]}/ACESII_{Id[var]}_l3_langmuir_fixed_fullCal.cdf"
density_dict = stl.loadDictFromFile(density_file)

# Getting DERPA1 Temperature Data
DERPA1_file = rf"{DataPaths.ACES_data_folder}/L3/DERPA/{rocket[var]}/ACESII_{Id[var]}_l3_ERPA1_fullCal.cdf"
DERPA1_dict = stl.loadDictFromFile(DERPA1_file)

# Getting DERPA2 Temperature Data
DERPA2_file = rf"{DataPaths.ACES_data_folder}/L3/DERPA/{rocket[var]}/ACESII_{Id[var]}_l3_ERPA2_fullCal.cdf"
DERPA2_dict = stl.loadDictFromFile(DERPA2_file)

n_e = []
T = []
delta_phi = []
sweep_Epoch = []
chi_square = []
plasma_potential = []
for sweep_no in tqdm(range(sweep_min, sweep_max[var], sweep_step)):

    try:
        # Getting the Fixed Langmuir Probe Data
        idx = np.abs(density_dict['Epoch'][0] - Epoch[sweep_no]).argmin()
        n_guess = deepcopy(density_dict['ni'][0][idx])

        # Getting the DERPA1 Data
        idx = np.abs(DERPA1_dict['Epoch'][0] - Epoch[sweep_no]).argmin()
        T_guess_1 = (e / k) * deepcopy(DERPA1_dict['temperature'][0][idx])

        # Getting the DERPA2 Data
        idx = np.abs(DERPA2_dict['Epoch'][0] - Epoch[sweep_no]).argmin()
        T_guess_2 = (e / k) * deepcopy(DERPA2_dict['temperature'][0][idx])

        # Picking Between DERPA1 and DERPA2
        T_guess = T_guess_2

        # I-V Curves for the Data
        V = V_sweeps[sweep_no] - np.abs(V_floating[sweep_no])
        I = I_sweeps[sweep_no] - np.min(I_sweeps[sweep_no])  # Adding the ion offset

        # Bin equal voltages and average the currents
        # (We only want one y-value for every x-value)
        V_unique, inv_idx = np.unique(V, return_inverse=True)
        I_avg = np.zeros_like(V_unique)
        for i in range(len(V_unique)):
            I_avg[i] = np.mean(I[inv_idx == i])
        V_reduced = V_unique
        I_reduced = I_avg

        # Only include data greater than or equal to -1 V
        idxs = np.where(V_reduced > x_data_min)[0]
        V = V_reduced[idxs]
        I = I_reduced[idxs]

        # Only include data greater than or equal to the floating potential
        idxs = np.where(V < 0)[0]
        V = V[idxs]
        I = I[idxs]

        idxs = np.where(I < 2500)[0]
        V = V[idxs]
        I = I[idxs]


        # Defining the Model Current function
        def fit_func(V_reduced, n, T_e, delta_phi):
            return (A / 4) * 1e9 * n * e * np.sqrt((8 * k * T_e) / (np.pi * m_e)) * np.exp(
                e * (V_reduced - delta_phi) / (k * T_e))


        # Curve Fitting
        p0 = [n_guess, T_guess, delta_phi_guess]
        bound_values = tuple(np.array(
            [[n_guess / n_tol, n_guess * n_tol],
             [T_guess - T_tol, T_guess + T_tol],
             [0.01, 3]]).T)

        # bound_values = tuple(np.array(
        #     [[1e9, 1e13],
        #      [T_guess - T_tol, T_guess + T_tol],
        #      [0.01, 3]]).T)

        params, cov = curve_fit(fit_func, V, I, p0=p0, bounds=bound_values, maxfev=10000000)
        x_data_fit = V
        y_data_fit = fit_func(x_data_fit, *params)
        chi_squared = np.sum((I - y_data_fit) ** 2 / y_data_fit)

        # Plotting Data against Model
        fig, ax = plt.subplots()
        ax.plot(V_reduced, I_reduced, 'o', label=f'Sweep {sweep_no}')
        ax.plot(V, I, 'o', color='r', label='Fit Data')
        ax.plot(x_data_fit, y_data_fit, '-', color='g',
                label=f'Fitted Data\n n = {round(params[0], 3)}\n T = {round(params[1], 3)} \n delta_phi_guess = {round(params[2], 5)} \n chi_squared = {round(chi_squared, 5)}')
        ax.set_xlabel('Voltage w.r.t phi_floating (V)')
        ax.set_ylabel('Electron Current (nA)')
        ax.set_ylim(-100, 10000)
        ax.set_xlim(-2, 1)
        fig.suptitle(f'Sweep {sweep_no} with Model Curve \n{Epoch[sweep_no]}')
        ax.legend()
        fig.savefig(rf"{DataPaths.ACES_data_folder}/L3/Langmuir/{rocket[var]}/plots/fit{sweep_no}.png")
        plt.close()

        # if chi_squared <= 100:
        n_e.append(params[0])
        T.append(params[1])
        delta_phi.append(params[2])
        sweep_Epoch.append(Epoch[sweep_no])
        chi_square.append(chi_squared)

    except:
        print(f'Error{sweep_no}')
        n_e.append(np.nan)
        T.append(np.nan)
        delta_phi.append(np.nan)
        sweep_Epoch.append(Epoch[sweep_no])
        chi_square.append(np.nan)

        continue

# Prepare data dict for CDF
# data_dict_output = {
#     'n_e': [np.array(n_e), {'UNITS': 'm^-3', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'Density'}],
#     'T': [(k/e)*np.array(T), {'UNITS': 'eV', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'Temperature'}],
#     'delta_phi': [np.array(delta_phi), {'UNITS': 'V', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'Delta Phi'}],
#     'plasma_potential': [np.array(delta_phi) + np.array(np.abs(V_floating[[i for i in range(sweep_min, sweep_max[var], sweep_step)]])), {'UNITS': 'V', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'Plasma Potential'}],
#     'Epoch': [np.array(sweep_Epoch), {}],
#     'chi_squared': [np.array(chi_square, dtype='float64'), {'DEPEND_0': 'Epoch', 'LABLAXIS': 'Chi Squared'}],
#     }

data_dict_output = {
    'n_e': [np.array(n_e), {'UNITS': 'm^-3', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'Density'}],
    'T': [(k / e) * np.array(T), {'UNITS': 'eV', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'Temperature'}],
    'delta_phi': [np.array(delta_phi), {'UNITS': 'V', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'Delta Phi'}],
    'plasma_potential': [np.array(delta_phi) + np.array(np.abs(V_floating[[i for i in range(len(delta_phi))]])),
                         {'UNITS': 'V', 'DEPEND_0': 'Epoch', 'LABLAXIS': 'Plasma Potential'}],
    'Epoch': [np.array(sweep_Epoch), {}],
    'chi_squared': [np.array(chi_square, dtype='float64'), {'DEPEND_0': 'Epoch', 'LABLAXIS': 'Chi Squared'}],
}

# Output path
fileOutPath = rf"{DataPaths.ACES_data_folder}/L3/Langmuir/{rocket[var]}/ACESII_{Id[var]}_LP_FitParameters.cdf"
stl.outputDataDict(outputPath=fileOutPath, data_dict=data_dict_output)