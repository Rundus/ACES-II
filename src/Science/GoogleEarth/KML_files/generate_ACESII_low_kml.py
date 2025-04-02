#!/usr/bin/env python
'''Generate a KML string that matches the altitudemode example.

References:
http://code.google.com/apis/kml/documentation/kmlreference.html#gxaltitudemode
http://code.google.com/apis/kml/documentation/kmlfiles/altitudemode_reference.kml
'''

# --- Load the spatial data ---
import spaceToolsLib as stl
data_dict = stl.loadDictFromFile(r'C:\Data\ACESII\attitude\low\ACESII_36364_Attitude_Solution.cdf')


# Format the data into the KML scheme
coordinates = [(data_dict["Long"][0][idx],
                data_dict["Lat"][0][idx],
                data_dict["Alt"][0][idx]) for idx in range(len(data_dict['Alt'][0]))]

import simplekml
kml = simplekml.Kml()

ls = kml.newlinestring(name='ACESII_Low_Flyer_Trajectory',
                  coords=coordinates)
ls.altitudemode = simplekml.AltitudeMode.relativetoground
ls.style.linestyle.width = 5
ls.style.linestyle.color= simplekml.Color.blue

kml.save(r'C:\Users\cfelt\Desktop\projects\GoogleEarth\ACESII\ACESII_Low_Flyer.kml')