# feral-cat-otways-2017-SMR

This is the data and code to build MLE spatial mark-resight models (sighting-only) for feral cats in the Otway Ranges, VIC, AUS - from 1 survey in 2017. This is not the exact code used for https://doi.org/10.1016/j.biocon.2019.108287, but an updated script with the following changes:
 1) blurry tabby cats correctly placed as tn - mark status unknown (thanks to Joanne Potts for picking this error up)
 2) camera-trap grids seperated by being different sessions, rather than using a habitat mask with covariate
 3) "fastproximity = FALSE" arguement in models with time covariates as well as models which are being compared to time covariates (necessary due to change in default settings in the new secr package) 
NOTE: these changes have no effect on the results - all results presented in Rees et al. 2019 are correct. 

If you would like to use this data in another analysis, please contact me. 
