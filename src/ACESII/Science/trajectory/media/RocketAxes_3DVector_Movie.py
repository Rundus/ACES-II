# --- csv_to_cdf_attitude.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Turn the .cdf files of the TRICE attitude data into cdf files
import numpy as np

# --- --- --- --- ---
from ACESII_code.myImports import *
start_time = time.time()
# --- --- --- --- ---

# --- --- --- ---
# --- TOGGLES ---
# --- --- --- ---
wRocket = 4

unSpinRocket = True # toggle to "reverse spin" the rocket to see how effective the RingCore_L1_to_L2_despin.py code is

# --- --- --- ---
# --- import ---
# --- --- --- ---
import pyIGRF
from matplotlib import pyplot as plt, animation
from ACESII_code.class_var_func import color, RotationAboutAxes

print(color.BOLD + color.CYAN + 'RocketAxes_3DVector_Movie.py' + color.END + color.END)


def ACESIIplotting():

    rocketAttrs,b,c = ACES_mission_dicts()

    # --- Create Attitude Data Dict ---
    prgMsg('Collecting Attitude data')
    inputFile = glob(f'{ACES_data_folder}attitude\\{fliers[wRocket-4]}\\*Attitude_Solution*')
    data_dicts_attitude = loadDictFromFile(inputFile[0], {})
    Done(start_time)

    prgMsg('Preparing Data')

    from src.ACESII.Science.trajectory.plotting.plottingStyles import AttitudeMovie
    def sphere2cart(r, theta, phi):
        return [
            r * np.sin(np.radians(theta)) * np.cos(np.radians(phi)),
            r * np.sin(np.radians(theta)) * np.sin(np.radians(phi)),
            r * np.cos(np.radians(theta))
        ]



    # --- --- --- --- --- --- --- --- -
    # --- PREPARE DATA FOR PLOTTING ---
    # --- --- --- --- --- --- --- --- -

    # convert all Epoch data to datetimes to display them
    Epoch = np.array([pycdf.lib.datetime_to_tt2000(tme) for tme in data_dicts_attitude["Epoch"][0]])
    timeOffset = pycdf.lib.datetime_to_tt2000(dt.datetime(2022, 11, 20, 17, 20, 00, 107800)) - Epoch[0]  # number of nanoseconds from 17:20:00 each was launched
    datetimesEpoch = np.array([pycdf.lib.tt2000_to_datetime(Epoch[i] + timeOffset).strftime("%H:%M:%S:%f") for i in range(len(Epoch))])

    Alt = data_dicts_attitude['Alt'][0]
    Lat = data_dicts_attitude['Latgd'][0]
    Long = data_dicts_attitude['Long'][0]
    X_Az = data_dicts_attitude['X_Az'][0]
    X_El = data_dicts_attitude['X_El'][0]
    Y_Az = data_dicts_attitude['Y_Az'][0]
    Y_El = data_dicts_attitude['Y_El'][0]
    Z_Az = data_dicts_attitude['Z_Az'][0]
    Z_El = data_dicts_attitude['Z_El'][0]

    # Convert all the axes data into cartesian vectors
    xAxisData = np.array([sphere2cart(1, 90 - X_El[i],  X_Az[i]) for i in range(len(X_Az))])
    yAxisData = np.array([sphere2cart(1, 270 - Y_El[i], 180 + Y_Az[i]) for i in range(len(Y_Az))])
    zAxisData = np.array([sphere2cart(1, 90 - Z_El[i],  Z_Az[i]) for i in range(len(Z_Az))])

    Done(start_time)

    if unSpinRocket:

        prgMsg('Unspinnig Rocket Axes')

        initialPhase = Y_Az[12000]
        spinFreq = 0.6442441031179716
        Epoch_seconds = np.array([(tme - Epoch[0])/1E9 for tme in Epoch])

        Rolls = np.array([1 * np.degrees(2 * np.pi * spinFreq * tme) + initialPhase for tme in Epoch_seconds])
        Yaws = np.array([0 for tme in Epoch_seconds])
        Pitchs = np.array([0 for tme in Epoch_seconds])

        yAxisData = np.array([ np.matmul(RotationAboutAxes(Rolls[i], *xAxisData[i]), yAxisData[i]) for i in range(len(yAxisData))])
        zAxisData = np.array([ np.matmul(RotationAboutAxes(Rolls[i], *xAxisData[i]), zAxisData[i]) for i in range(len(zAxisData))])

        Done(start_time)


    prgMsg('Initializing Plot')


    # Get the magnetic field data
    # -- Output order forpyIGRF.igrf_value ---
    # [0] Declination (+ E | - W)
    # [1] Inclination (+ D | - U)
    # [2] Horizontal Intensity
    # [3] North Comp (+ N | - S)
    # [4] East Comp (+ E | - W)
    # [5] Vertical Comp (+ D | - U)
    # [6] Total Field
    date = 2022 + 323 / 365  # Corresponds to 11/20/2022
    IGRF = [pyIGRF.igrf_value(Lat[i], Long[i], Alt[i], date) for i in range(len(Epoch))]
    Bgeo = np.array([[IGRF[i][3], IGRF[i][4], -1*IGRF[i][5]] for i in range(len(Epoch))])
    bgeo = np.array([b/np.linalg.norm(b) for b in Bgeo]) # convert geomag field into unit vector

    # plot the attitude data at one time
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim([-1.5, 1.5])
    ax.set_zlim([-1.5, 1.5])

    def get_arrow(i, data):
        x, y, z = 0, 0, 0
        u = data[i][0]
        v = data[i][1]
        w = data[i][2]
        return x, y, z, u, v, w


    # --- INITIALIZE ANIMATION ---

    testval = 0

    # initial direction vectors
    xQ = ax.quiver(*get_arrow(testval, xAxisData), color='red', label='$X$')
    yQ = ax.quiver(*get_arrow(testval, yAxisData), color='green', label='$Y$')
    zQ = ax.quiver(*get_arrow(testval, zAxisData), color='blue', label='$Z$')
    BgeoQ = ax.quiver(*get_arrow(testval, bgeo), color='cyan', label='$IGRF$')
    # global xQ, yQ, zQ, BgeoQ

    title = ax.set_title(f'ACESII {rocketAttrs.rocketID[wRocket - 4]} \n {datetimesEpoch[0]}')
    ax.set_xlabel('North')
    ax.set_ylabel('East')
    ax.invert_yaxis() # makes it east instead of west
    ax.set_zlabel('Up')
    plt.legend(loc='best')
    Done(start_time)

    # --- --- --- --- --- --- --
    # --- ANIMATION FUNCTION ---
    # --- --- --- --- --- --- --
    def quiver_data_to_segments(X, Y, Z, u, v, w, length=1):
        segments = (X, Y, Z, X + v * length, Y + u * length, Z + w * length)
        segments = np.array(segments).reshape(6, -1)
        return [[[x, y, z], [u, v, w]] for x, y, z, u, v, w in zip(*list(segments))]

    def animatePlot(i):

        xQ.set_segments(quiver_data_to_segments(*get_arrow(i, xAxisData),length=1))
        yQ.set_segments(quiver_data_to_segments(*get_arrow(i, yAxisData), length=1))
        zQ.set_segments(quiver_data_to_segments(*get_arrow(i, zAxisData), length=1))
        BgeoQ.set_segments(quiver_data_to_segments(*get_arrow(i, bgeo), length=1))

        # update Epoch title
        title.set_text(f'ACESII {rocketAttrs.rocketID[wRocket - 4]} \n {datetimesEpoch[i]}')


    # --- --- --- --- --- --- -
    # --- PERFORM ANIMATION ---
    # --- --- --- --- --- --- -
    prgMsg('Creating Movie')
    ax.view_init(40, -40)
    # locations = [i for i in range(len(X_Az))]
    # locations = [i for i in range(12000, 17497)]
    locations = [i for i in range(12000, 13000,2)]
    anim = animation.FuncAnimation(fig=fig, func=animatePlot, interval=1000 / AttitudeMovie.fps, frames=locations)
    anim.save(f'C:\\Data\ACESII\\trajectories\\trajectory_plots\\ACESII_{rocketAttrs.rocketID[wRocket - 4]}_Trajectory.mp4')
    Done(start_time)



# --- --- --- ---
# --- EXECUTE ---
# --- --- --- ---
ACESIIplotting()