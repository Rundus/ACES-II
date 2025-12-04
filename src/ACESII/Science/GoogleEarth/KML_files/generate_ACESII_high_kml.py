#!/usr/bin/env python
'''Generate a KML string that matches the altitudemode example.

References:
http://code.google.com/apis/kml/documentation/kmlreference.html#gxaltitudemode
http://code.google.com/apis/kml/documentation/kmlfiles/altitudemode_reference.kml
'''

# --- Load the spatial data ---
import spaceToolsLib as stl
data_dict_spatial = stl.loadDictFromFile(r'C:\Data\ACESII\attitude\high\ACESII_36359_Attitude_Solution.cdf')


# Format the data into the KML scheme
coordinates = [(data_dict_spatial["Long"][0][idx],
                data_dict_spatial["Lat"][0][idx],
                data_dict_spatial["Alt"][0][idx]) for idx in range(len(data_dict_spatial['Alt'][0]))]

import simplekml
kml = simplekml.Kml()

ls = kml.newlinestring(name='ACESII_High_Flyer_Trajectory',
                  coords=coordinates)
ls.altitudemode = simplekml.AltitudeMode.relativetoground
ls.style.linestyle.width = 10
ls.style.linestyle.color= simplekml.Color.red

kml.save(r'C:\Users\cfelt\Desktop\projects\GoogleEarth\ACESII\ACESII_High_Flyer.kml')