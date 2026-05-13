#!/usr/bin/env python
'''Generate a KML string that matches the altitudemode example.

References:
http://code.google.com/apis/kml/documentation/kmlreference.html#gxaltitudemode
http://code.google.com/apis/kml/documentation/kmlfiles/altitudemode_reference.kml
'''

# --- Load the spatial data ---
import spaceToolsLib as stl
import numpy as np
from src.ACESII.data_tools.data_paths import DataPaths
data_dict = stl.loadDictFromFile(rf'/home/connor/Data/MODELS/ACESII_ionosphere/spatial_environment/spatial_environment.cdf')


# Format the data into the KML scheme
coordinates = [(data_dict["grid_long"][0][idx1][idx2],
                data_dict["grid_lat"][0][idx1][idx2],
                data_dict["grid_alt"][0][idx1][idx2]) for idx1 in range(0, len(data_dict['grid_alt'][0]),50) for idx2 in range(len(data_dict['grid_alt'][0][0]))]

import simplekml
kml = simplekml.Kml()

ls = kml.newlinestring(name='ACESII_SimulationRegion',
                  coords=coordinates,
                 visibility=1)
ls.altitudemode = simplekml.AltitudeMode.relativetoground
ls.style.linestyle.width = 5
ls.style.linestyle.color= simplekml.Color.cyan


kml.save(rf'/home/connor/PycharmProjects/ACES-II/src/ACESII/Science/GoogleEarth/KML_files/kml/ACESII_simulation_region.kml')