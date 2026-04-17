Launcher.sh executes the commands to run climatology of O2 from a BGC-Argo dataset.

Details at: 
https://github.com/CarolinaAmadio/bitsea_complement/wiki/Yr_Climfloat_netcdf.py

In the launcher I(Camadio) have used a subset of the SUPERFLOAT dataset 2012-2024-01-01 as ONLINE_REPO in order to avoid recents profiles in the profilelist.

This script computes the yearly vertical climatologies of a specific BGC-Argo variable across various Mediterranean sub-basins.

Input: BGC-Argo Float Profiles: Data for a specific variable (--variable) and time range selected for each sub-basin.

Method: Profile Interpolation: For each sub-basin, profiles are interpolated onto the model depth levels to standardize measurements across all sub-basins. Computation of Averages & Standard Deviations: The script calculates vertical averages and standard deviations for each sub-basin, providing a robust representation of the data.

Output:


NetCDF Files:

yr_Avg_.nc: Contains the mean values of the interpolated profiles for each sub-basin.
yr_Std_.nc: Contains the standard deviations of the profiles, offering insight into the variability across the sub- basins.

optional:
PNG Figures: The script generates vertical profile plots for each sub-basin, which are saved as PNG files, offering a visual comparison of the profiles.
CSV Files: The interpolated profile data for each sub-basin is saved as CSV files, enabling further analysis and validation.
