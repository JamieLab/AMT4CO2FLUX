 #!/usr/bin/env python3
# ---------------------------------
# Code created by Daniel Ford (d.ford@exeter.ac.uk) - 18/08/2022
# Converts the AMT28 20 minute bin files into a netcdf for input into FluxEngine
#-------- Version 3.0
# Implemented Donlon (2002) cool skin correction
# Re-added fixed skin correciton of 0.17 K
# Added weighted mean to ERA5 data importing (mean of 4 nearest pixels, weighted by distance)
# Changed netcdf file writing so all variables in the Pandas data table are saved.
# Set time as seconds since 1970, consistent with 3 hour files.
#--------  Version 2.0
# Implementing COARE 3.5 cool skin bias
# Removed fixed skin correction approach
#--------  Version 1.0
# Initial Version
# ---------------------------------
import pandas as pd
from netCDF4 import Dataset
import datetime
import numpy as np
from NOAA_COARE import coare35vn as co
import cdsapi
from os.path import exists
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

#Function to convert the supplied dates (YYYY-MM-DD HH:MM:SS) into a number of seconds since 1990-01-01 00:00:00
#so time values can be saved in the netcdf file.
def datecalc(dates):
    d = datetime.datetime(1970,1,1,0,0,0)
    out = np.empty((len(dates),1))
    for i in range(len(dates)):
        out[i] = (datetime.datetime.strptime(dates[i],'%d/%m/%Y %H:%M')-d).total_seconds()
    return out

def dateback(date):
    d = []
    for i in range(date.shape[0]):
        d.append(datetime.datetime(1970,1,1,0,0,0)+datetime.timedelta(seconds=int(date[i])))
    return d

def hour_rounder(t):
    t = t[0]
    # Rounds to nearest hour by adding a timedelta hour if minute >= 30
    return (t.replace(second=0, microsecond=0, minute=0, hour=t.hour)
               +datetime.timedelta(hours=t.minute//30))

def era_retrieve(time,lat,lon):
    # time = time[0]
    #print(time.strftime('%Y%m%d_%H'))
    file = 'D:/DATA/ERA5/ERA5_HOURLY_'+time.strftime('%Y%m%d_%H')+'.nc'
    # print(time.strftime("%Y"))
    # print(time.strftime("%m"))
    # print(time.strftime("%d"))
    # print(time.strftime("%H:%M"))
    if exists(file) == 0:
        c = cdsapi.Client()
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'variable': [
                    'boundary_layer_height',
                    'mean_surface_downward_long_wave_radiation_flux',
                    'mean_surface_downward_short_wave_radiation_flux',
                ],
                'year': time.strftime("%Y"),
                'month': time.strftime("%m"),
                'day': time.strftime("%d"),
                'time': time.strftime("%H:00"),
                'format': 'netcdf',
            },
            file)
    else:
        print('Exists!')
    blh,Rs,Rl = ecmwf_load(file,lat,lon)

    return blh,Rs,Rl

def ecmwf_load(file,lat,lon):
    #print(lat)
    #print(lon)
    c = Dataset(file,'r')
    latg = np.array(c.variables['latitude'][:])
    long = np.array(c.variables['longitude'][:])
    # ECMWF longitude values are in degress E - therefore anything negative has 360 added to bring it
    # to correct longitude.
    if lon < 0:
        lon = lon+360
    la,dila = calc_nearest(lat,latg)
    lo,dilo = calc_nearest(lon,long)
    di = calc_dist_matrix(np.array(dilo),np.array(dila))
    blh = np.average(np.array(c.variables['blh'][0,int(la[0]):int(la[1])+1,int(lo[0]):int(lo[1])+1]),weights=di)
    Rl = np.average(np.array(c.variables['msdwlwrf'][0,int(la[0]):int(la[1])+1,int(lo[0]):int(lo[1])+1]),weights=di)
    Rs = np.average(np.array(c.variables['msdwswrf'][0,int(la[0]):int(la[1])+1,int(lo[0]):int(lo[1])+1]),weights=di)

    c.close()
    print(blh,Rs,Rl)
    return blh,Rs,Rl

def calc_dist_matrix(lo,la):
    la = np.transpose(np.tile(la,(2,1)))
    lo = np.tile(lo,(2,1))
    di = 1/np.sqrt(la**2 + lo**2)**2
    #print(di)
    return di

def calc_nearest(val,grid):
    dif = np.abs(grid - val)
    near_v = np.empty((2))
    difs = np.empty((2))
    l = np.where(np.min(dif) == dif)
    near_v[0] = l[0][0]
    difs[0] = np.min(dif)
    dif[int(near_v[0])] = np.nan
    l = np.where(np.nanmin(dif) == dif)
    near_v[1] = l[0][0]
    difs[1] = np.nanmin(dif)
    l = np.argsort(near_v)
    return np.array(near_v[l]),difs[l]

#-------------------------------------------------------------------------------
#Defining input and output file names
in_filename = 'DATA/AMT29_table_20min_v3_20221205_S3A_ADDED_CCISST_ADDED.txt'
out_filename = 'DATA/AMT28_table_20min_version2_i.nc'

#-------------------------------------------------------------------------------
#Loading the data table - delimiter is a tab (\t) and the header rows is 24 (needs to be adjusted depending on version)
data = pd.read_table(in_filename,sep='\t')#,skiprows=26)

#For creating the netcdf file
#Get header names from the pandas data table
header = list(data.columns)
#Time converison function to get times as seconds since 1970-01-01 00:00:00
time = datecalc(data[header[0]])
#Get the number of rows in the data table so this can be used to correctly define the size
#of variables in the netcdf
l = data.shape[0]
#-------------------------------------------------------------------------------
#Downloading and loading ERA5 data to input into COARE algorithm to get skin correction
blh = np.empty((l))
sw = np.empty((l))
lw = np.empty((l))

for i in range(l):
    blh[i],sw[i],lw[i] = era_retrieve(hour_rounder(dateback(time[i])),data['lat'][i],data['lon'][i])
data['Boundary_Layer_Height'] = blh
data['shortwave'] = sw
data['longwave'] = lw

out = co.coare35vn(data['U10n_mean'],data['Tair_mean'],data['RH_mean'],data['Tskin_mean']-273.15,P=data['Pres_mean'],Rs=sw,Rl=lw,zu=10,zt=18,zq=18,lat=data['lat'],zi=blh,jcool=0)
data['Coare_Cool_Skin'] = out[:,14]

sats = ['S3A_Tskin_mean','CCI_avhrr_on_metop_a_mean','CCI_s3a+b_mean','CCI_merged_mean'] #Calculating NOAA Coare for satellite Tskin values
for sat in sats:
    out2 = co.coare35vn(data['U10n_mean'],data['Tair_mean'],data['RH_mean'],data[sat]-273.15,P=data['Pres_mean'],Rs=sw,Rl=lw,zu=10,zt=18,zq=18,lat=data['lat'],zi=blh,jcool=0)
    data['Coare_Cool_Skin_'+sat] = out2[:,14]
data['Donlon_2002_cool_skin'] = -(-0.14 - 0.30 * np.exp(- (data['U10n_mean'] / 3.7)))

#-------------------------------------------------------------------------------
# Estimate of Tskin from Tdepth for just cool skin runs
data['T_skine_co'] = data['SST_mean'] - data['Coare_Cool_Skin']
data['T_skine_don'] = data['SST_mean'] - data['Donlon_2002_cool_skin']
data['T_skine_fixed'] = data['SST_mean'] - 0.17
#-------------------------------------------------------------------------------
#Reanalysing in situ fCO2sw at Tdepth to subskin SST
data['T_subskin_co'] = data['Tskin_mean'] + data['Coare_Cool_Skin'] - 273.15 #Converting the skin SST from the ISAR from kelvin to degrees, and then to subskin SST
data['T_subskin_don']= data['Tskin_mean'] + data['Donlon_2002_cool_skin'] - 273.15 # This becomes subskin temp
data['T_subskin_fixed'] = data['Tskin_mean'] + 0.17 - 273.15

for sat in sats:
    data[sat+'_T_subskin_co'] = data[sat] + data['Coare_Cool_Skin_'+sat] - 273.15
    data[sat+'_T_subskin_don']= data[sat] + data['Donlon_2002_cool_skin'] - 273.15
    data[sat+'_T_subskin_fixed'] = data[sat] + 0.17 - 273.15

#Calculate fco2_sw_subskin by converting the fco2_sw_depth using SST_depth and the subskin temp following Takahashi et al. (1993) using coefficients in Wanninkhof et al. 2022
data['fCO2_sw_subskin_co'] = data['fCO2_sw_mean'] * np.exp((0.0413*(data['T_subskin_co'] - data['SST_mean']))) #- (4.35e-5 * (data['T_subskin_co']**2 - data['SST_mean']**2)) )
data['fCO2_sw_subskin_don'] = data['fCO2_sw_mean'] * np.exp((0.0413*(data['T_subskin_don'] - data['SST_mean']))) #- (4.35e-5 * (data['T_subskin_don'] **2 - data['SST_mean']**2)) )
data['fCO2_sw_subskin_fixed'] = data['fCO2_sw_mean'] * np.exp((0.0413*(data['T_subskin_fixed'] - data['SST_mean'])))# - (4.35e-5 * (data['T_subskin_fixed'] **2 - data['SST_mean']**2)) )
for sat in sats:
    data[sat+'_fCO2_sw_subskin_co'] = data['fCO2_sw_mean'] * np.exp((0.0413*(data[sat+'_T_subskin_co'] - data['SST_mean']))) #- (4.35e-5 * (data['T_subskin_co']**2 - data['SST_mean']**2)) )
    data[sat+'_fCO2_sw_subskin_don'] = data['fCO2_sw_mean'] * np.exp((0.0413*(data[sat+'_T_subskin_don'] - data['SST_mean']))) #- (4.35e-5 * (data['T_subskin_don'] **2 - data['SST_mean']**2)) )
    data[sat+'_fCO2_sw_subskin_fixed'] = data['fCO2_sw_mean'] * np.exp((0.0413*(data[sat+'_T_subskin_fixed'] - data['SST_mean'])))# - (4.35e-5 * (data['T_subskin_fixed'] **2 - data['SST_mean']**2)) )
# fig = plt.figure()
# gs = GridSpec(2,2, figure=fig, wspace=0.15,hspace=0.2,bottom=0.08,top=0.95,left=0.06,right=0.95) # Setup a 2x2 subplot grid
# ax1 = fig.add_subplot(gs[0,0]) # Set top row into a single plot
# ax2 = fig.add_subplot(gs[0,1])
# ax3 = fig.add_subplot(gs[1,0]) # Set bottom row into two plots
# ax4 = fig.add_subplot(gs[1,1])
#
# ax1.plot(data['lat'],out[:,14])
# ax2.plot(data['lat'],data['fCO2_sw_mean'] - data['fCO2_sw_subskin'])
# ax3.plot(data['lat'],data['Tskin_mean']+data['Coare_Cool_Skin']-273.15 - data['SST_mean'])
# ax4.plot(data['lat'],data['Tskin_mean']-273.15,'r-')
# ax4.plot(data['lat'],data['Tair_mean'],'b-')
# plt.show()
#-------------------------------------------------------------------------------
#Create second and thrid order wind products
data['U10n_2'] = data['U10n_mean']**2
data['U10n_3'] = data['U10n_mean']**3
#-------------------------------------------------------------------------------
#Creating netcdf file

header = list(data.columns)
#Start creating netcdf file in write mode
outp = Dataset(out_filename,'w',format='NETCDF4_CLASSIC')
#Set the dimension of time - this in the only dimension of the data
outp.createDimension('time',1)
outp.createDimension('longitude',1)
outp.createDimension('latitude',l)
#Manually input the time variable as the values are in there own variable
val = outp.createVariable('time','f4',('latitude'))
#Set units - this is particularly important for this variable
val.units = 'Seconds since 1990-01-01 00:00:00'
val[:] = time
val = outp.createVariable('latitude','f4',('latitude'))
val[:] =np.array(range(l))/l
val = outp.createVariable('longitude','f4',('longitude'))
val[:] = 1
# Run through each column in the pandas data table, using the column header as the variable name in the netCDF
# and save the data to the netCDF. Units are not supplied as these can be retrieved from the intial data file
# if required.
for i in range(1,len(header)):
    val = outp.createVariable(header[i],'f4',('latitude','longitude'))
    val[:] = np.array(data[header[i]])
    print(val.shape)
# Close the netcdf once all variables are written
outp.close()
