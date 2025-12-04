# # --- TOFdiscriminantPlot.py ---
# # --- Author: C. Feltman ---
# # DESCRIPTION: scipt to make a plot investigating the solutions of the TOF equation


# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
import numpy as np
import matplotlib.pyplot as plt
from ACESII_code.class_var_func import m_e,q0
# Make the Data
RE = 6371
Zmin = 0
Zmax = 10*RE # in km
N = 50000
ZaccRange = np.linspace(Zmin,Zmax,N)
pitchs = [0, 30, 60, 90]
# pitchs = [0,30, 60]
colatitude = 90 - 70
BobsMean = 5E-5
zobs = 150/RE
B0 = 3.14E-5


plotDeltaZBetaAndGamma = False
plotDiscriminant = False
plotTwoSolutions = False
plotFinalEnergySolution = True


def beta(z,alpha,colat,Bobs):
    y = ((RE**3)*(np.sin(np.radians(alpha))**2) * 3*B0 *np.sqrt(1 + 3*np.cos(np.radians(colat))*np.cos(np.radians(colat))))/(Bobs*((RE + z)**4))
    return y

def gamma(z,alpha, colat,Bobs):
    Bacc = B0*((RE/(RE + z))**3) * np.sqrt(1 + 3*np.cos(np.radians(colat))*np.cos(np.radians(colat)))
    BB = Bacc/Bobs
    y = 1 / (1 - np.sin(np.radians(alpha))*np.sin(np.radians(alpha)) * BB)
    return y

gammaVals,betaVals = [],[]

for ptch in pitchs:
    gammaVals.append(np.array([1 / (gamma(z, ptch, colatitude, BobsMean)) for z in ZaccRange]))
    betaVals.append(np.array([beta(z, ptch, colatitude, BobsMean) * (z - zobs) for z in ZaccRange]))


colors=['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:cyan']

if plotDeltaZBetaAndGamma:

    # Make the plots
    plt.figure()

    # --- PLOT 1 ---
    plt.subplot(1,2,1)
    xData = ZaccRange/RE
    plt.title('PLOT OF ' + r'$\beta \Delta z$'+'\n'+
                 r'$\theta = $' + f'{colatitude}' + r'$^{\circ}$, ' +
                 r'B0=' + f'{B0}, ' + '$B_{obs}=$' + f'{BobsMean}')


    for i in range(len(pitchs)):
        plt.plot(xData,betaVals[i],color=colors[i])


    plt.xlabel('z [$R_{E}$]')
    plt.ylabel(r'$\beta\Delta z = \Delta z\frac{\sin^{2}(\alpha_{obs})}{ B_{obs}} (3B_{0} \sqrt{1 + 3 \cos^{2}\theta} \frac{R_{E}^{3}}{(R_{E} + z)^{4}})$')
    plt.legend([r'$\alpha$= '+f'{ptch}' for ptch in pitchs])

    # --- PLOT 2 ---
    plt.subplot(1,2,2)
    plt.title('PLOT OF ' + r'$2 / \gamma$'+'\n'+
                 r'$\theta = $' + f'{colatitude}' + r'$^{\circ}$, ' +
                 r'B0=' + f'{B0}, ' + '$B_{obs}=$' + f'{BobsMean}')


    for i in range(len(pitchs)):
        plt.plot(xData,gammaVals[i],color=colors[i])

    plt.xlabel('z [$R_{E}$]')
    plt.ylabel(r'$2 / \gamma = 2(1 - \sin^{2}\alpha_{obs} (B_{acc}/B_{obs}))$')
    plt.legend([r'$\alpha$= '+f'{ptch}' for ptch in pitchs])
    # plt.ylim(0.5,3)

    plt.show()

if plotDiscriminant:

    # --- --- --- --- --- --- --- -
    # --- Plot the discriminant ---
    # --- --- --- --- --- --- --- -
    discriminant = []
    for i in range(len(pitchs)):
        discriminant.append(
            np.array(
                [ gammaVals[i][j] + betaVals[i][j] - (betaVals[i][j])**2 for j in range(len(gammaVals[i]))]
            )
        )

    fig,ax = plt.subplots()
    for i in range(len(discriminant)):
        ax.plot(xData,discriminant[i],color = colors[i])
    plt.title('PLOT OF ' + r'$2/\gamma + \Delta z \beta - (\Delta z \beta)^{2}$'+'\n'+
                 r'$\theta = $' + f'{colatitude}' + r'$^{\circ}$, ' +
                 r'B0=' + f'{B0}, ' + '$B_{obs}=$' + f'{BobsMean}')
    ax.set_ylabel(r'$2/\gamma + \Delta z \beta - (\Delta z \beta)^{2}$')
    ax.set_xlabel('z [$R_{E}$]')
    plt.legend([r'$\alpha$= '+f'{ptch}' for ptch in pitchs])
    ax.set_xlim(-0.1,4)
    plt.show()

if plotTwoSolutions:
    # --- --- --- --- --- --- --- --- --
    # --- Plot the solution function ---
    # --- --- --- --- --- --- --- --- --
    plt.figure()
    plt.subplot(1,2,1)
    solution = []
    for i in range(len(pitchs)):
        solution.append(
            np.array(
                [ gammaVals[i][j] + betaVals[i][j] - 1 for j in range(len(gammaVals[i]))]
            )
        )


    for i in range(len(solution)):
        plt.plot(xData,solution[i],color = colors[i])
    plt.title('PLOT OF ' + r'$2/\gamma + \Delta z \beta - 1$'+'\n'+
                 r'$\theta = $' + f'{colatitude}' + r'$^{\circ}$, ' +
                 r'B0=' + f'{B0}, ' + '$B_{obs}=$' + f'{BobsMean}')
    plt.ylabel(r'$2/\gamma + \Delta z \beta - 1$')
    plt.xlabel('z [$R_{E}$]')
    plt.legend([r'$\alpha$= '+f'{ptch}' for ptch in pitchs])

    # -----------------------
    plt.subplot(1,2,2)
    solution = []
    for i in range(len(pitchs)):
        solution.append(
            np.array(
                [ gammaVals[i][j] + betaVals[i][j] + 1 for j in range(len(gammaVals[i]))]
            )
        )


    for i in range(len(solution)):
        plt.plot(xData,solution[i],color = colors[i])
    plt.title('PLOT OF ' + r'$2/\gamma + \Delta z \beta + 1$'+'\n'+
                 r'$\theta = $' + f'{colatitude}' + r'$^{\circ}$, ' +
                 r'B0=' + f'{B0}, ' + '$B_{obs}=$' + f'{BobsMean}')
    plt.ylabel(r'$2/\gamma + \Delta z \beta + 1$')
    plt.xlabel('z [$R_{E}$]')
    plt.legend([r'$\alpha$= '+f'{ptch}' for ptch in pitchs])


    plt.show()

if plotFinalEnergySolution:
    plt.figure()

    ZaccChoice = [0.5*RE, RE, 1.5*RE,2*RE]
    zaccels = ['$0.5R_{E}$', '$R_{E}$', '$1.5R_{E}$','$2R_{E}$']
    wPitch = 10
    colat = 20
    t0 = 0
    zobs = 400
    N = 10000
    timeMin, timeMax = 0, 4

    def totalEnergyData(t,t0,wZacc,alpha,colat,Bobs):
        discriminant = np.sqrt(2/gamma(wZacc,alpha,colat,Bobs)+(wZacc-zobs)*beta(wZacc,alpha,colat,Bobs) - ((wZacc-zobs)*beta(wZacc,alpha,colat,Bobs))**2)
        y = (1/q0)*(( ((wZacc-zobs)*beta(wZacc,alpha,colat,Bobs)))+(2/gamma(wZacc,alpha,colat,Bobs))+discriminant)/( (beta(wZacc,alpha,colat,Bobs))**2  )*(2*m_e/((t-t0)**2))
        return y

    timeRange = np.linspace(timeMin, timeMax, N)
    totalEnergies = []
    for wZaccChoice in ZaccChoice:
        totalEnergies.append(np.array([ totalEnergyData(tme,t0,wZaccChoice,wPitch,colat,BobsMean) for tme in timeRange]))


    # --- Make the plot ---
    for i, data in enumerate(totalEnergies):
        plt.plot(timeRange, data, color=colors[i])

    plt.title('PLOT OF ' + r'$W_{tot}$' + '\n' +
              r'$\alpha_{obs}$' + f'= {wPitch}' + r'$^{\circ}$' + '\n'
              r'$\theta = $' + f'{colat}' + r'$^{\circ}$, ' +
              r'B0=' + f'{B0}, ' + '$B_{obs}=$' + f'{BobsMean},' + ' $z_{obs}$=' + f'{zobs}km')
    plt.ylabel(r'$W_{tot}$'+' [eV]')
    plt.xlabel('t [seconds]')
    plt.legend([r'$z_{acc}$= ' + f'{num}' for num in zaccels])
    plt.ylim(20,1000)
    plt.show()
