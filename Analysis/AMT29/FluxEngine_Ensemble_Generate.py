#!/usr/bin/env python3
# ---------------------------------
# Code created by Daniel Ford (d.ford@exeter.ac.uk) - 23/08/2022
# Create ensemble uncertainity data to use as input to Fluxengine
# ------- Version 2.0
# Added Donlon, COARE and fixed skin corrections to ensemble generation
#--------  Version 1.0
# - Initial Version
# ---------------------------------

from netCDF4 import Dataset
import datetime
import numpy as np
from random import normalvariate

def noise_func(val,noise,perc,ensem):
    out = np.empty((ensem))
    out.shape
    if perc == 1:
        noise = val*noise
    for i in range(ensem):
        out[i] = val + normalvariate(0,noise)
    return out

def noise_array(val,noise,perc,ensem,file = '', var = ''):
    if perc == 2:
        c = Dataset(file,'r')
        print(var)
        noise = np.squeeze(np.array(c.variables[var]))
        print(noise)
        c.close()
    var = np.empty((val.shape[0],ensem))
    for i in range(val.shape[0]):
        if perc == 2:
            # print(perc)
            # print(noise[i])
            # print(val[i])
            v = noise_func(val[i],noise[i],perc,ensem)
        else:
            v = noise_func(val[i],noise,perc,ensem)
        var[i,:] = v
    var[:,0] = val
    return var

in_file = 'DATA/AMT28_table_20min_version2_i.nc'
out_file = 'DATA/AMT28_table_20min_version2.nc'
ensembles = 100
parameters = ['Tskin_mean', 'U10n_mean','Pres_mean','Salinity_mean','xCO2atm_mean','fCO2_sw_subskin_co','fCO2_sw_mean','T_subskin_co','SST_mean','fCO2_sw_subskin_don','T_subskin_don','fCO2_sw_subskin_fixed','T_subskin_fixed'
    ,'S3A_Tskin_mean','S3A_Tskin_mean_T_subskin_co','S3A_Tskin_mean_T_subskin_don','S3A_Tskin_mean_T_subskin_fixed','S3A_Tskin_mean_fCO2_sw_subskin_co','S3A_Tskin_mean_fCO2_sw_subskin_don','S3A_Tskin_mean_fCO2_sw_subskin_fixed'
    ,'CCI_avhrr_on_metop_a_mean','CCI_avhrr_on_metop_a_mean_T_subskin_co','CCI_avhrr_on_metop_a_mean_T_subskin_don','CCI_avhrr_on_metop_a_mean_T_subskin_fixed','CCI_avhrr_on_metop_a_mean_fCO2_sw_subskin_co',
    'CCI_avhrr_on_metop_a_mean_fCO2_sw_subskin_don','CCI_avhrr_on_metop_a_mean_fCO2_sw_subskin_fixed'
    ,'CCI_s3a+b_mean','CCI_s3a+b_mean_T_subskin_co','CCI_s3a+b_mean_T_subskin_don','CCI_s3a+b_mean_T_subskin_fixed','CCI_s3a+b_mean_fCO2_sw_subskin_co',
    'CCI_s3a+b_mean_fCO2_sw_subskin_don','CCI_s3a+b_mean_fCO2_sw_subskin_fixed'
    ,'CCI_merged_mean','CCI_merged_mean_T_subskin_co','CCI_merged_mean_T_subskin_don','CCI_merged_mean_T_subskin_fixed','CCI_merged_mean_fCO2_sw_subskin_co',
    'CCI_merged_mean_fCO2_sw_subskin_don','CCI_merged_mean_fCO2_sw_subskin_fixed'
    ,'T_skine_co','T_skine_don','T_skine_fixed']
unc = [0.2,0.03,0,0,1,4,4,0.1,0.1,4,0.1,4,0.1,
    0.2,0.1,0.1,0.1,4,4,4,
    0.2,0.1,0.1,0.1,4,4,4,
    0.2,0.1,0.1,0.1,4,4,4,
    0.2,0.1,0.1,0.1,4,4,4,
    0.1,0.1,0.1]
# if perc == 0 then the fixed value uncertainty is applied as noise.
# if perc == 1 then a percentage of the mean value is applied as the uncertainty
# if perc == 2 then a array of same length of the mean values is used for the uncertainty loaded from the netCDF.
perc = [2,1,0,0,0,0,0,0,0,0,0,0,0,
    2,0,0,0,0,0,0,
    2,0,0,0,0,0,0,
    2,0,0,0,0,0,0,
    2,0,0,0,0,0,0,
    0,0,0]
vari = ['Tskin_uncertainty','','','','','','','','','','','','',
    'S3A_Tskin_unc','','','','','','',
    'CCI_avhrr_on_metop_a_unc','','','','','','',
    'CCI_s3a+b_unc','','','','','','',
    'CCI_merged_unc','','','','','','',
    '','','']

c = Dataset(in_file ,'r')
l = np.array(c.variables['lat'])

d = Dataset(out_file ,'w',format='NETCDF4_CLASSIC')
d.createDimension('time',1)
d.createDimension('longitude',ensembles)
d.createDimension('latitude',l.shape[0])

val = d.createVariable('time','f4',('time'))
val[:] = 1
val = d.createVariable('latitude','f4',('latitude'))
val[:] = np.array(range(l.shape[0]))/l.shape[0]
val = d.createVariable('longitude','f4',('longitude'))
val[:] = np.array(range(ensembles))

for i in range(len(parameters)):
    print(vari[i])
    var = noise_array(np.squeeze(np.array(c.variables[parameters[i]])),unc[i],perc[i],ensembles,file=in_file,var=vari[i])
    out = d.createVariable(parameters[i],'f4',('latitude','longitude'))
    out[:] = var
    if parameters[i] == 'U10n_mean':
        out = d.createVariable('U10n_2','f4',('latitude','longitude'))
        out[:] = var**2
        out = d.createVariable('U10n_3','f4',('latitude','longitude'))
        out[:] = var**3
c.close()
d.close()
