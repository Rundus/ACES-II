import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# File path
file_path = "/home/connor/PycharmProjects/ACES-II/src/ACESII/Science/Solar_wind_conditions/ACE_satellite_data/ace_mag_data_.txt"

# Read the data while skipping header/comment lines
cols = [
    "YR", "MO", "DA", "HHMM", "Julian", "SecDay", "Status",
    "Bx", "By", "Bz", "Bt", "Lat", "Long"
]

df = pd.read_csv(
    file_path,
    delim_whitespace=True,
    comment="#",
    names=cols,
    skiprows=20
)

# Remove bad/missing data
df = df[df["Status"] == 0]
df = df[(df["Bx"] > -999) & (df["By"] > -999) & (df["Bz"] > -999)]

# Create UTC datetime column
df["HHMM"] = df["HHMM"].astype(str).str.zfill(4)

df["UTC"] = pd.to_datetime(
    df["YR"].astype(str) + "-" +
    df["MO"].astype(str).str.zfill(2) + "-" +
    df["DA"].astype(str).str.zfill(2) + " " +
    df["HHMM"].str[:2] + ":" +
    df["HHMM"].str[2:],
    utc=True
)

# Plot
fig, ax = plt.subplots(figsize=(12, 5))

ax.plot(df["UTC"], df["Bx"], label="Bx")
ax.plot(df["UTC"], df["By"], label="By")
ax.plot(df["UTC"], df["Bz"], label="Bz")

ax.set_xlabel("UTC Time")
ax.set_ylabel("Magnetic Field (nT)")
ax.set_title("ACE Magnetometer Data")
ax.legend()

# Format time axis
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
fig.autofmt_xdate()

plt.tight_layout()
plt.show()