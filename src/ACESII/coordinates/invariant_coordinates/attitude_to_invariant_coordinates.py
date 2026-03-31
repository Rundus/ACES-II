# --- attitude_to_invariant_coorrdinates.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: loads the ACESII attitude data and calculates the invariant coordinates
# for a specific chosen invariant altitude



#################
# --- TOGGLES ---
#################
just_print_file_names_bool = False
rocket_str = 'high'
dict_file_path ={ # FORMAT: Data Name: [Str modifier to ACESII Data Folder Path, Which Datafile Indices in directory [[High flyer], [Low flyer]]]
    f'attitude':['/', [[0],[0]]],
}
outputData = True
targetProjectionAltitude = 150 # altitude you want to project B (in km). Should be ~100km
useAndoya = True

#################
# --- IMPORTS ---
#################
from src.ACESII.data_tools.my_imports import *
import time
start_time = time.time()
from spacepy import pycdf
from spacepy import coordinates as coord
coord.DEFAULTS.set_values(use_irbem=False, itol=5)  # maximum separation, in seconds, for which the coordinate transformations will not be recalculated. To force all transformations to use an exact transform for the time, set ``itol`` to zero.
from spacepy.time import Ticktock #used to determine the time I'm choosing the reference geomagentic field

lat_to_meter = 111.319488 # 1 deg latitude to kilometers on Earth
def long_to_meter(lat):
    return 111.319488 * np.cos(np.radians(lat))


def attitude_to_invariant_coordinates(data_dicts):

        # --- Load the Data ---
        data_dict_attitude = deepcopy(data_dicts[0])


        # --- prepare the output ---
        data_dict_output = {
            'Epoch':deepcopy(data_dict_attitude['Epoch']),
            'Alt': deepcopy(data_dict_attitude['Alt']),
            'Lat': deepcopy(data_dict_attitude['Lat']),
            'Long': deepcopy(data_dict_attitude['Long']),

            'ILong': [[], {'DEPEND_0':'Epoch','LABLAIXS': 'invariant ionospheric longitude','UNITS':'Degrees'}],
            'ILat': [[], {'DEPEND_0':'Epoch','LABLAIXS': 'invariant ionospheric latitude','UNITS':'Degrees'}],

            'mAlt': [[], {'DEPEND_0':'Epoch','LABLAIXS': 'magnetic altitude','UNITS':'Ree'}],
            'mLat': [[], {'DEPEND_0':'Epoch','LABLAIXS': 'magnetic latitude','UNITS':'Degrees'}],
            'mLong': [[], {'DEPEND_0':'Epoch','LABLAIXS': 'magnetic longitude','UNITS':'Degrees'}],

            'Lat_km': [[], {'DEPEND_0':'Epoch','LABLAIXS': 'latitude','UNITS':'km'}],
            'Long_km': [[], {'DEPEND_0':'Epoch','LABLAIXS': 'longitude','UNITS':'km'}],

            'mILat': [[], {'DEPEND_0':'Epoch','LABLAIXS': 'magnetic invariant ionospheric latitude','UNITS':'Degrees'}],
            'mILong': [[], {'DEPEND_0':'Epoch','LABLAIXS': 'magnetic invariant ionospheric longitude','UNITS':'Degrees'}],

            'ILat_km': [[], {'DEPEND_0':'Epoch','LABLAIXS': 'magnetic invariant ionospheric latitude','UNITS':'km'}],
            'ILong_km': [[], {'DEPEND_0':'Epoch','LABLAIXS': 'magnetic invariant ionospheric longitude','UNITS':'km'}],

            'distance_from_launch_site': [[], {'DEPEND_0':'Epoch','LABLAIXS': 'Distance from launc site','UNITS':'km'}]
        }
        Epoch = deepcopy(data_dict_attitude['Epoch'][0])

        ################################
        # --- Calculate IGRF B-Field ---
        ################################

        stl.prgMsg('Getting CHAOS Field')
        # -- Output order forpyIGRF.igrf_value ---
        # [0] Declination (+ E | - W)
        # [1] Inclination (+ D | - U), should be ~78deg for what we're doing
        # [2] Horizontal Intensity
        # [3] North Comp (+ N | - S)
        # [4] East Comp (+ E | - W)
        # [5] Vertical Comp (+ D | - U)
        # [6] Total Field
        pos = np.array([data_dict_attitude['Lat'][0], data_dict_attitude['Long'][0], data_dict_attitude['Alt'][0]/stl.m_to_km]).T
        Bgeo = stl.CHAOS(pos[:,0],pos[:,1], pos[:,2], Epoch)

        stl.Done(start_time)

        #################################
        # --- I-LAT I-LONG PROJECTION ---
        #################################

        stl.prgMsg('Projecting B-Field')
        intersection_pos = [] # [[long,lat,alt],...]
        bhat = []

        # Perform Triangulation Projection

        for i in range(len(Epoch)):
            #Coordiantes reported in (Lat (x) , Long (y), Alt (z))
            loc = pos[i] # IGRF vector Location, should be rocket coordinates
            bDir = Bgeo[i] # IGRF vector Direction. -1 added in third elemental due to down being positive in IGRF given
            bDirNorm = bDir / np.linalg.norm(bDir) # Normalize IGRF to get its direction only. This will make t larger, but that's fine
            bhat.append(bDirNorm)

            # Determine the Delta-Latitude and Longitutde
            # Theta_dip = np.arctan(bDirNorm[1]/np.abs(bDirNorm[2])) # latitude declination
            # Theta_inc = np.arctan(bDirNorm[0]/np.abs(bDirNorm[2])) # longitude declination
            h = loc[2] - targetProjectionAltitude

            deltaLat = (1/lat_to_meter) * h*np.tan(np.radians(90-Theta_dip))
            deltaLong = (1/long_to_meter(loc[0])) * h*np.tan(np.radians(Theta_inc))

            intersection_pos.append([loc[0] + deltaLong, loc[1] + deltaLat, 0])

        intersection_pos = np.array(intersection_pos)
        data_dict_output['ILat'][0] = np.array(intersection_pos[:,0])
        data_dict_output['ILong'][0] = np.array(intersection_pos[:,1])
        stl.Done(start_time)

        ############################################
        # --- CONVERT TO GEOMAGNETIC COORDINATES ---
        ############################################
        stl.prgMsg('Converting to Geomagnetic Coordinates')

        ISOtime = [tVal.isoformat() for tVal in Epoch]

        # Convert the rocket position to geomagnetic coordinates
        cvals_GDZ = coord.Coords(pos, 'GDZ', 'sph')
        cvals_GDZ.ticks = Ticktock(ISOtime, 'ISO')
        cvals_GDZ_MAG = cvals_GDZ.convert('MAG', 'sph')
        data_dict_output['mAlt'][0] = cvals_GDZ_MAG.radi
        data_dict_output['mLat'][0] = cvals_GDZ_MAG.lati
        data_dict_output['mLong'][0] = cvals_GDZ_MAG.long

        # Convert the invariant coordinates to geomagnetic coordinates
        pos_invar = np.array([intersection_pos[:, 1], intersection_pos[:, 0], intersection_pos[:, 2]]).T # in Lat, Long, Alt format
        cvals_GDZ = coord.Coords(pos_invar, 'GDZ', 'sph')
        cvals_GDZ.ticks = Ticktock(ISOtime, 'ISO')
        cvals_GDZ_MAG_intersects = cvals_GDZ.convert('MAG', 'sph')
        data_dict_output['mILat'][0] = cvals_GDZ_MAG_intersects.lati
        data_dict_output['mILong'][0] = cvals_GDZ_MAG_intersects.long

        stl.Done(start_time)

        # #######################
        # # --- CONVERT TO KM ---
        # #######################
        #
        # # Determine the distance in km from some lat/long reference point
        # if useAndoya:
        #     refLat = rocketAttrs.Andoya_Space_Lat_Long[0]
        #     refLong = rocketAttrs.Andoya_Space_Lat_Long[1]
        # else:
        #     refLat = 0
        #     refLong = 0
        #
        # geodeticLatIntersects_km = [[], []]
        # geodeticLongIntersects_km = [[], []]
        # LatDS_km = [[], []]
        # LongDS_km = [[], []]
        #
        #
        # for i in range(rangelen):
        #     # Convert lat/long data to KM
        #     for j in range(len(geodeticLatIntersects[i])):
        #         geodeticLatIntersects_km[i].append(lat_to_meter * (geodeticLatIntersects[i][j] - refLat))
        #         geodeticLongIntersects_km[i].append(long_to_meter(refLat) * (geodeticLongIntersects[i][j] - refLong))
        #         LatDS_km[i].append(lat_to_meter * (LatDS[i][j] - refLat))
        #         LongDS_km[i].append(long_to_meter(refLat) * LongDS[i][j] - long_to_meter(refLat) * refLong)
        #
        # # Determine the magntitude distance from the launch point
        # distanceFromLaunchPoint = [
        #     [np.sqrt((LatDS_km[0][i])**2 + (LongDS_km[0][i])**2) for i in range(len(LatDS_km[0]))],
        #     [np.sqrt((LatDS_km[1][i]) ** 2 + (LongDS_km[1][i]) ** 2) for i in range(len(LatDS_km[1]))]
        # ]

        if outputData:
            file_name = f'ACESII_{ACESII.fliers_dict[rocket_str]}_invariant_coordinates.cdf'
            outputPath = f'{DataPaths.ACES_data_folder}/coordinates/invariant_coordinates/{rocket_str}/' + file_name
            stl.outputDataDict(outputPath, data_dict_output)






#################
# --- EXECUTE ---
#################
DataClasses().ACEII_file_executor(attitude_to_invariant_coordinates,
                                  dict_file_path,
                                  rocket_str,
                                  just_print_file_names_bool)



