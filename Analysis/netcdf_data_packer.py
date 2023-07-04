#!/usr/bin/env python3

import pandas as pd
from netCDF4 import Dataset
import numpy as np
import datetime

# files = ['AMT28/DATA/AMT28_table_20min_version2_i.nc', 'AMT28/DATA/FLUXENGINE_OUT.nc']
# outfile = 'AMT28_20min_data_v3.nc'

files = ['AMT29/DATA/AMT28_table_20min_version2_i.nc', 'AMT29/DATA/FLUXENGINE_OUT.nc']
outfile = 'AMT29_20min_data_v3.nc'

var1 = ['Pres_mean',
    'fCO2_sw_mean',
    'fCO2_sw_subskin_co',
    'fCO2_sw_subskin_don',
    'fCO2_sw_subskin_fixed',
    'Salinity_mean',
    'SST_mean',
    'T_skine_co',
    'T_skine_don',
    'T_skine_fixed',
    'time',
    'Tskin_mean',
    'Tskin_uncertainty',
    'T_subskin_co',
    'T_subskin_don',
    'T_subskin_fixed',
    'U10n_2',
    'U10n_3',
    'U10n_mean',
    'xCO2atm_mean',
    'CO2Flux_EC',
    'Err_CO2flux_EC']
lon1 = ['Atmospheric pressure at 20 m above sea level',
    'Fugacity of underway seawater CO2 at SST depth',
    'Fugacity of underway seawater CO2 corrected to subskin temperature',
    'Fugacity of underway seawater CO2 corrected to subskin temperature',
    'Fugacity of underway seawater CO2 corrected to subskin temperature',
    'Sea surface salinity at depth',
    'Sea surface temperature at depth',
    'Skin temperature estimated from SST depth',
    'Skin temperature estimated from SST depth',
    'Skin temperature estimated from SST depth',
    'UTC time',
    'Skin temperature measured by the ISAR instrument',
    'Skin temperature uncertainity from ISAR uncertainty model',
    'Subskin temperature estimated from ISAR skin temperature',
    'Subskin temperature estimated from ISAR skin temperature',
    'Subskin temperature estimated from ISAR skin temperature',
    'Second moment of wind speed (U^2)',
    'Third moment of wind speed (U^3)',
    '10-m neutral wind speed',
    'Mixing ratio of dry atmospheric CO2',
    'Eddy Covariance CO2 flux',
    'Eddy Covaraince CO2 flux uncertainty']
com1 = ['',
    'fCO2 (sw,depth)',
    'fCO2 (sw,depth) corrected to subskin temperature estimated from ISAR skin temperature and NOAA COARE3.5',
    'fCO2 (sw,depth) corrected to subskin temperature estimated from ISAR skin temperature and Donlon et al. 2002',
    'fCO2 (sw,depth) corrected to subskin temperature estimated from ISAR skin temperature and Donlon et al. 1999',
    'Sea surface salinity measured underway at ~6 m',
    'Sea surface temperature measured underway at ~6m (Tdepth)',
    'Skin temperature estimated from Tdepth using NOAA COARE3.5 skin deviation (assuming no warm layer)',
    'Skin temperature estimated from Tdepth using Donlon et al. 2002 skin deviation (assuming no warm layer)',
    'Skin temperature estimated from Tdepth using Donlon et al. 1999 skin deviation (assuming no warm layer)',
    'Middle time point of 20 minute mean',
    'Skin temperature measured by the ISAR instrument',
    'Skin temperature uncertainty estimated by the ISAR uncertainty model described in Wimmer and Robinson (2016)',
    'Subskin temperature estimated from the ISAR skin temperature using the NOAA COARE 3.5 skin deviation',
    'Subskin temperature estimated from the ISAR skin temperature using the Donlon et al. 2002 skin deviation',
    'Subskin temperature estimated from the ISAR skin temperature using the Donlon et al. 1999 skin deviation',
    'Second moment of wind speed (U^2) for gas transfer parameterisations',
    'Third moment of wind speed (U^3) for gas transfer parameterisations',
    '10-m wind speed, bulk value derived from measured wind speed using COARE3.5 model',
    'Dry mixing ratio of atmospheric CO2',
    'Air-sea CO2 flux estimated by the eddy covariance system',
    'Uncertainity in the air-sea CO2 flux by the eddy convariance system, estimated following Dong et al. 2021']
uni1 = ['hPa',
    'uatm',
    'uatm',
    'uatm',
    'uatm',
    'psu',
    'degC',
    'K',
    'K',
    'K',
    'Seconds since 1970-01-01 00:00:00',
    'K',
    'K',
    'K',
    'K',
    'K',
    '',
    '',
    'ms-1',
    'ppm',
    'mmol m-2 d-1',
    'mmol m-2 d-1']

inps = Dataset(files[0],'r')
time = np.array(inps['time'])
print(time)
print(np.where(np.isnan(time) == 1))
outs = Dataset(outfile,'w')
outs.date_created = datetime.datetime.now().strftime(('%d/%m/%Y'))
outs.created_by = 'Daniel J. Ford (d.ford@exeter.ac.uk)'
outs.created_from = 'Data created from ' + files[0] + ' and ' + files[1]
outs.packed_with = 'Data created with netcdf_data_packer.py'
outs.version = '3'
outs.createDimension('time',time.shape[0])

for i in range(len(var1)):
    print(i)
    var = outs.createVariable(var1[i],'f8',('time'))
    v =np.array(inps[var1[i]][:])
    print(v)
    var[:] = v
    var.standard_name = lon1[i]
    var.comment = com1[i]
    var.units = uni1[i]

var = outs.createVariable('latitude','f8',('time'))
var[:] = np.array(inps['lat'])
var.standard_name = 'Latitude'
var.units = 'Degrees North'
var = outs.createVariable('longitude','f8',('time'))
var[:] = np.array(inps['lon'])
var.standard_name = 'Longitude'
var.units = 'Degrees East'
inps.close()
gas = ['Nightingale2000','Wanninkhof2014','Ho2006','Yang2022']
temp = ['vt_coare','vt_don','vt_fixed','vt_coare_sk','vt_don_sk','vt_fixed_sk']
skin = ['COARE3.5', 'Donlon et al. 2002','Donlon et al. 1999','COARE3.5', 'Donlon et al. 2002','Donlon et al. 1999']
warm = ['','','',' (No Warm Layer)',' (No Warm Layer)',' (No Warm Layer)']
inps = Dataset(files[1],'r')
vt = np.array(inps['vt'])[:,0:24]
print(vt.shape)
vts = np.array(inps['vts'])[:,0:24]
t = 0
l = 0
for i in range(24):
    print('bulk_'+temp[t]+'_'+gas[l])
    var = outs.createVariable('bulk_'+temp[t]+'_'+gas[l],'f8',('time'))
    var[:] = vt[:,i]
    var.standard_name = 'Bulk CO2 flux with temperature gradients using '+skin[t]+' cool skin and ' + gas[l] + ' gas transfer' + warm[t]
    var.units = 'mmol m-2 d-1'

    print('bulk_'+temp[t]+'_'+gas[l])
    var = outs.createVariable('bulk_'+temp[t]+'_'+gas[l]+'_unc','f8',('time'))
    var[:] = vts[:,i]
    var.standard_name = 'Bulk CO2 flux uncertainty with temperature gradients using '+skin[t]+' cool skin and ' + gas[l] + ' gas transfer' + warm[t]
    var.units = 'mmol m-2 d-1'
    t = t+1
    if t == 6:
        l = l+1
        t=0

nvt = np.array(inps['nvt'])
nvts = np.array(inps['nvts'])

for i in range(4):
    print('bulk_nvt_'+gas[i])
    var = outs.createVariable('bulk_nvt_'+gas[i],'f8',('time'))
    var[:] = nvt[:,i]
    var.standard_name = 'Bulk CO2 flux with no temperature gradients ' + gas[i] + ' gas transfer'
    var.units = 'mmol m-2 d-1'

    var = outs.createVariable('bulk_nvt_'+gas[i]+'_unc','f8',('time'))
    var[:] = nvts[:,i]
    var.standard_name = 'Bulk CO2 flux uncertainty with no temperature gradients ' + gas[i] + ' gas transfer'
    var.units = 'mmol m-2 d-1'
inps.close()
outs.close()
