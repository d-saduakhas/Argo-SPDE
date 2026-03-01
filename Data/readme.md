# Data

Code for fetching and pre-processing the Argo profiles dataset. This code can be used independently to scrape the data and convert it to .mat files.

## Pre-processing pipeline

1. `Data_fetching_cleaning_org_git.ipynb` - Python notebook to automatically download data from the GDAC data repository (creates all sub-folders for each year in a hierarchical way)

2. `health_check_downloaded2.ipynb` - checks if all files were downloaded correctly

3. `Merging_lower_upper_text.ipynb` - identifies files where variables are in lower or upper case due to different encoding standards (pre-processing to merge them into one file later)

4. `preProcess.m` - modified code from Kuusela et al. (2018), cleans and converts all netCDF files into .mat format

5. `interpolateToPressureLevel.m` - modified code from Kuusela et al. (2018), interpolates the data from profiles to a specified pressure level of interest

6. `estimateMeanField.m` - fits the local polynomial regression

7. `subtractMeanField.m` - subtracts the mean field computed in the previous step

## Climatology

`RG_climatology/RG_ArgoClim_Temperature_2019.nc` - can be downloaded from the [Roemmich-Gilson Argo Climatology page](https://sio-argo.ucsd.edu/RG_Climatology.html)

## Total interpolated profiles

| Pressure level (dbar) | Profiles  |
|------------------------|-----------|
| 10                     | 1,296,996 |
| 1000                   | 1,182,667 |
