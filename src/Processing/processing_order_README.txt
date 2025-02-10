--- NOTES ON THE DATA PROCESSING FOR THE ACESII DATA ---

The order that files must be executed

[0] csv_to_cdf_trajectories.py
[1] Tad_to_tmCDF.py
[2] tmCDF_to_L0.py
[3] L0_to_L1.py


--- Langmuir Probes ---
[4] L1_to_Langmuir.py

--- ESAs ---
[4] L1_&_Traject_to_Traject_ILatILong.py
[5] L1_to_L2.py
[6] L2_to_DistFunc.py


....lol this document really needs to be updated