All the code used for fetching, pre-processing Argo profiles dataset.

This code can be used independently to scrape the data and convert it .mat files.

Data_fetching_cleaning_org_git.ipynb - Python notebook, fully automatic (creates all the sub-folders for each year in hierarchical way) to download data from the GADR data repository directly

health_check_downloaded2.ipynb - execute to check if all files were downloaded correctly

Merging_lower_upper_text.ipynb - gets the names of files where variables are either in lower or upper case, due to different encryption standards(pre-processing to merge them into one file later)

preProcess.m - modified code from Kuusela et.al (2018), to clean and convert all netCDF files into .mat format

interpolateToPressureLevel.m - modified code from Kuusela et.al (2018) to interpolate the data from profiles to specified pressure level of interest.

estimateMeanField.m - fit the local polynomial regression

subtractMeanField.m - subtract the mean field computed in the previous step

Data/RG_climatology/RG_ArgoClim_Temperature_2019.nc - can be downloaded from

total interpolated profiles:

pressure level 10 - 1296996

300 - 1349863

1000 - 1182667
