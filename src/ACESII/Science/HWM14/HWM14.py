from pyhwm2014 import HWM14

# Define parameters
year = 2023
day_of_year = 150  # Approximately May 30
universal_time = 12.0  # 12:00 UT (noon)
altitude_km = 300.0  # 300 km altitude
latitude = 40.0  # 40°N
longitude = -105.0  # 105°W
ap_index = 10  # Geomagnetic activity index

# Retrieve wind values
hwm14 = HWM14(
    alt=altitude_km,
    altlim=[altitude_km, altitude_km],
    altstp=1,
    year=year,
    day=day_of_year,
    ut=universal_time,
    glat=latitude,
    glon=longitude,
    ap=[-1, ap_index],
    option=1,
    verbose=False
)

# Access results
zonal_wind = hwm14.Uwind[0]  # m/s
meridional_wind = hwm14.Vwind[0]  # m/s

print(f"Zonal wind: {zonal_wind:.2f} m/s")
print(f"Meridional wind: {meridional_wind:.2f} m/s")