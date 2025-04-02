#!/usr/bin/env python
'''Generate a KML string that matches the altitudemode example.

References:
http://code.google.com/apis/kml/documentation/kmlreference.html#gxaltitudemode
http://code.google.com/apis/kml/documentation/kmlfiles/altitudemode_reference.kml
'''

# --- Load the spatial data ---
import spaceToolsLib as stl
data_dict_spatial = stl.loadDictFromFile(r'C:\Data\physicsModels\alfvenic_auroral_acceleration_AAA\spatial_environment\spatial_environment.cdf')


# Format the data into the KML scheme
coordinates = [(data_dict_spatial["simLong"][0][idx], data_dict_spatial["simLat"][0][idx], data_dict_spatial["simAlt"][0][idx]) for idx in range(len(data_dict_spatial['simAlt'][0]))]

import simplekml
kml = simplekml.Kml()

ls = kml.newlinestring(name='Bgeo_example',
                  coords=coordinates)
ls.altitudemode = simplekml.AltitudeMode.relativetoground
ls.style.linestyle.width = 5
ls.style.linestyle.color= simplekml.Color.blue

kml.save(r'C:\Users\cfelt\Desktop\projects\GoogleEarth\ACESII\pyKML_Bgeo_example.kml')