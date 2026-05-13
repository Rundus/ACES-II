#!/usr/bin/env python
'''Generate a KML string that matches the altitudemode example.

References:
http://code.google.com/apis/kml/documentation/kmlreference.html#gxaltitudemode
http://code.google.com/apis/kml/documentation/kmlfiles/altitudemode_reference.kml
'''

# --- Load the spatial data ---
import spaceToolsLib as stl
from src.ACESII.data_tools.data_paths import DataPaths
from src.ACESII.mission_attributes import ACESII

for i in range(2):

    wFlyer = ACESII.fliers[i]
    data_dict_spatial = stl.loadDictFromFile(rf'{DataPaths.ACES_data_folder}/attitude/{wFlyer}/ACESII_{ACESII.fliers_dict[wFlyer]}_Attitude_Solution.cdf')

    # Format the data into the KML scheme
    coordinates = [(data_dict_spatial["Long"][0][idx],
                    data_dict_spatial["Lat"][0][idx],
                    data_dict_spatial["Alt"][0][idx]) for idx in range(len(data_dict_spatial['Alt'][0]))]

    import simplekml
    kml = simplekml.Kml()

    ls = kml.newlinestring(name=f'ACESII_{wFlyer}_Flyer_Trajectory',
                      coords=coordinates)
    ls.altitudemode = simplekml.AltitudeMode.relativetoground
    ls.style.linestyle.width = 10
    ls.style.linestyle.color= simplekml.Color.red if i == 0 else simplekml.Color.blue

    kml.save(rf'/home/connor/PycharmProjects/ACES-II/src/ACESII/Science/GoogleEarth/KML_files/kml/ACESII_{wFlyer}_Flyer.kml')