#!/usr/bin/env python3
# ---------------------------------
# Code created by Daniel Ford (d.ford@exeter.ac.uk) - 18/08/2022
# Code designed to run the FluxEngine configurations not/accounting for vertical temperature gradients
# and saving the data in a single netcdf, alongside the input data.
#-------- Version 1.2
# - Added fixed skin correction version back to script
# - Corrected issue with Yang COARE and Donlon fluxengine runs.
#-------- Version 1.1
# - Modifications to script to run fluxengine using Donlon 2002 and NOAA Coare skin corrections.
# - Modified to output all gas transfer parameterisations seperately, instead of mean of three.
# - Updated netcdf saving to include all data from orginial netcdf.
#--------  Version 1.0
# - Initial Version
# ---------------------------------
import pandas as pd
from fluxengine.core import fe_setup_tools as fluxengine
from netCDF4 import Dataset
import numpy as np

def run_fluxengine(config_file):
    returnCode, fe = fluxengine.run_fluxengine(config_file, 2010, 2010, singleRun=True,verbose=True)
    c = Dataset('DATA/FLUXENGINE_OUTPUTS/VT_CO/2010/01/OceanFluxGHG-month01-jan-2010-v0.nc','r')
    # Loading the flux data from the output location
    flux = np.squeeze(np.array(c.variables['OF'][:])/12.0107*1000) # Convert flux from g C m-2 d-1 to mmol m-2 d-1 (12.0107 = carbon wieght - 1000 conversion from mol to mmol)
    flux[flux < -2000] = np.nan # Remove all values below -2000 as these are fill values
    c.close()
    # Save the flux and k values in there variables
    flux_o = np.stack((flux[:,0],np.std(flux,axis=1)),axis=1)
    return flux_o


outfile = 'DATA/FLUXENGINE_OUT.nc'
#Loading the 20 minute EC flux and uncertianity data from the input file netcdf.
#Also loading the latitude data
input = 'DATA/AMT28_table_20min_version2_i.nc'

it = Dataset(input,'r')
#Squeezing the variables to remove extra dimensions
EC_flux = np.squeeze(np.array(it.variables['CO2Flux_EC'][:]))
headers = list(it.variables)

#Start running the FluxEngine configurations for with/without vertical temp gradients
#for the different gas transfer parameterisations
##### Need to write code to do an ensemble runs for uncertainity
#-------------------------------------------------------------------------------
# MIGHTINGALE 2000
#-------------------------------------------------------------------------------
# #Nightingale2000 with vertical temp COARE
flux_vt = run_fluxengine('FLUXENGINE_CONFIGS/NIGHT2000_VT_COARE.conf')

fl = run_fluxengine('FLUXENGINE_CONFIGS/NIGHT2000_VT_DON.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/NIGHT2000_VT_FIXED.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/NIGHT2000_VT_SKIN_ONLY_COARE.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/NIGHT2000_VT_SKIN_ONLY_DON.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/NIGHT2000_VT_SKIN_ONLY_FIXED.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

# #-------------------------------------------------------------------------------
# # Wanninkhof2013
# #-------------------------------------------------------------------------------

fl = run_fluxengine('FLUXENGINE_CONFIGS/WANNINK2013_VT_COARE.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/WANNINK2013_VT_DON.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/WANNINK2013_VT_FIXED.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/WANNINK2013_VT_SKIN_ONLY_COARE.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/WANNINK2013_VT_SKIN_ONLY_DON.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/WANNINK2013_VT_SKIN_ONLY_FIXED.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

# #-------------------------------------------------------------------------------
# # HO 2006
# #-------------------------------------------------------------------------------
fl = run_fluxengine('FLUXENGINE_CONFIGS/HO2006_VT_COARE.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/HO2006_VT_DON.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/HO2006_VT_FIXED.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/HO2006_VT_SKIN_ONLY_COARE.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/HO2006_VT_SKIN_ONLY_DON.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/HO2006_VT_SKIN_ONLY_FIXED.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

# #-------------------------------------------------------------------------------
# # YANG 2021
# #--------------------------------------------------------------------------------
fl = run_fluxengine('FLUXENGINE_CONFIGS/YANG2021_VT_COARE.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/YANG2021_VT_DON.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/YANG2021_VT_FIXED.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/YANG2021_VT_SKIN_ONLY_COARE.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/YANG2021_VT_SKIN_ONLY_DON.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/YANG2021_VT_SKIN_ONLY_FIXED.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

# #-------------------------------------------------------------------------------
# # SATELLITE RUNS
# #--------------------------------------------------------------------------------

#MERGED Dataset
fl = run_fluxengine('FLUXENGINE_CONFIGS/SAT/MERGED_WANNINK2013_VT_COARE.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/SAT/MERGED_WANNINK2013_VT_DON.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/SAT/MERGED_WANNINK2013_VT_FIXED.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

#METOP=A Dataset
fl = run_fluxengine('FLUXENGINE_CONFIGS/SAT/METOPA_WANNINK2013_VT_COARE.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/SAT/METOPA_WANNINK2013_VT_DON.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/SAT/METOPA_WANNINK2013_VT_FIXED.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

#S3A Dataset
fl = run_fluxengine('FLUXENGINE_CONFIGS/SAT/S3A_WANNINK2013_VT_COARE.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/SAT/S3A_WANNINK2013_VT_DON.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/SAT/S3A_WANNINK2013_VT_FIXED.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

#S3AB Dataset
fl = run_fluxengine('FLUXENGINE_CONFIGS/SAT/S3AB_WANNINK2013_VT_COARE.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/SAT/S3AB_WANNINK2013_VT_DON.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/SAT/S3AB_WANNINK2013_VT_FIXED.conf')
flux_vt = np.concatenate((flux_vt,fl),axis=1)

# ###-----------------------------------------------------------------------------
# #-------------------------------------------------------------------------------
# # NO VERTICAL TEMPERATURE GRADIENTS RUNS
# #-------------------------------------------------------------------------------
flux_nvt = run_fluxengine('FLUXENGINE_CONFIGS/NIGHT2000_NVT.conf')

fl = run_fluxengine('FLUXENGINE_CONFIGS/WANNINK2013_NVT.conf')
flux_nvt = np.concatenate((flux_nvt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/HO2006_NVT.conf')
flux_nvt = np.concatenate((flux_nvt,fl),axis=1)

fl = run_fluxengine('FLUXENGINE_CONFIGS/YANG2021_NVT.conf')
flux_nvt = np.concatenate((flux_nvt,fl),axis=1)
#-------------------------------------------------------------------------------
###EXTRACTING THE DATA AND SAVING TO NETCDF
l=flux_vt.shape[1]
vt = flux_vt[:,list(range(0,l,2))]
vts = flux_vt[:,list(range(1,l,2))] * 2
print(np.nanmean(vts/np.abs(vt)))
vts = np.abs(vt)*np.sqrt((vts/vt)**2 + (0.1)**2)
print(np.nanmean(vts/np.abs(vt)))
a = np.where(np.isnan(EC_flux))
vt[a,:] = np.nan # Where EC flux is nan, set vt flux to nan
#nvt = np.nanmean(flux_nvt,axis=1)
# Mean of nvt fluxes
nvt = flux_nvt[:,[0,2,4,6]]
nvts = flux_nvt[:,[1,3,5,7]]*2
nvts = np.abs(nvt)*np.sqrt((nvts/nvt)**2 + (0.1)**2)
nvt[a,:] = np.nan # vt includes nans of EC_flux and vt so use this to nan values

c = Dataset(outfile,'w','NETCDF4_CLASSIC')
c.createDimension('latitude',len(EC_flux))
c.createDimension('longitude',int(l/2))
c.createDimension('longitude2',4)
val = c.createVariable('longitude','f4',('longitude'))
val[:] = np.array(range(int(l/2)))
val = c.createVariable('longitude2','f4',('longitude2'))
val[:] = np.array(range(int(4)))
for i in range(len(headers)):
    if headers[i] != 'longitude':
        val = c.createVariable(headers[i],'f4',('latitude'))
        val[:] = np.squeeze(np.array(it.variables[headers[i]][:]))
val = c.createVariable('vt','f4',('latitude','longitude'))
val[:] = vt
val = c.createVariable('vts','f4',('latitude','longitude'))
val[:] = vts
val = c.createVariable('nvt','f4',('latitude','longitude2'))
val[:] = nvt
val = c.createVariable('nvts','f4',('latitude','longitude2'))
val[:] = nvts
c.close()
it.close()
