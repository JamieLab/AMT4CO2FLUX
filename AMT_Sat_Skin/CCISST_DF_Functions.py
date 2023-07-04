
#!/usr/bin/env python3

# CCI-SST functions for handling retrieving, loading and matching the daily data to in situ observations
# Uses the Copernicus Climate Data Store - attempts to retrieve the Skin SST from all available satellites
# in the archive or specific satellites can be specified.
# Produced by Daniel Ford (d.ford@exeter.ac.uk)
# ===============================================
# Version history
# ===============================================
# v1
# -> Initial version

import cdsapi
import datetime
import numpy as np
from netCDF4 import Dataset
import os
import glob
import zipfile
import matplotlib.pyplot as plt
data_save_loc = 'D:/Data/SST-CCI/daily_l3'
sats = np.array(['avhrr_on_noaa_07','avhrr_on_noaa_09','avhrr_on_noaa_11','avhrr_on_noaa_12','avhrr_on_noaa_14','avhrr_on_noaa_15','avhrr_on_noaa_16','avhrr_on_noaa_17',
    'avhrr_on_noaa_18','avhrr_on_noaa_19','avhrr_on_metop_a','avhrr_on_metop_b','atsr1_on_ers_1','atrs2_on_ers_2','aatsr_on_envisat','slstr_on_sentinel_3a','slstr_on_sentinel_3b'])
sat_dates = np.array([
    [datetime.datetime(1981,8,24,0,0,0), datetime.datetime(1985,2,18,0,0,0)], #AVHRR 7
    [datetime.datetime(1985,1,4,0,0,0),datetime.datetime(1988,11,7,0,0,0)], #AVHRR 9
    [datetime.datetime(1988,10,12,0,0,0),datetime.datetime(1994,9,13,0,0,0)], #AVHRR 11
    [datetime.datetime(1991,9,16,0,0,0),datetime.datetime(1998,12,14,0,0,0)], #AVHRR 12
    [datetime.datetime(1995,1,19,0,0,0),datetime.datetime(1999,12,31,0,0,0)], #AVHRR 14
    [datetime.datetime(1998,9,24,0,0,0),datetime.datetime(2009,12,31,0,0,0)], #AVHRR 15
    [datetime.datetime(2003,6,1,0,0,0),datetime.datetime(2006,12,31,0,0,0)], #AVHRR 16
    [datetime.datetime(2002,7,10,0,0,0),datetime.datetime(2009,12,31,0,0,0)], #AVHRR 17
    [datetime.datetime(2005,6,5,0,0,0),datetime.datetime(2009,12,31,0,0,0)], #AVHRR 18
    [datetime.datetime(2009,2,22,0,0,0),datetime.datetime(2018,12,31,0,0,0)], #AVHRR 19
    [datetime.datetime(2006,11,21,0,0,0),datetime.datetime(2021,9,30,0,0,0)], #AVHRR on METOPA
    [datetime.datetime(2021,10,1,0,0,0),datetime.datetime(2022,12,3,0,0,0)], #AVHRR on METOPB
    [datetime.datetime(1991,11,1,0,0,0),datetime.datetime(1996,1,9,0,0,0)], #ATSR on ERS1
    [datetime.datetime(1995,8,1,0,0,0),datetime.datetime(2003,6,22,0,0,0)], #ATSR on ERS2
    [datetime.datetime(2002,7,24,0,0,0),datetime.datetime(2012,4,8,0,0,0)], #AASTR on envisat
    [datetime.datetime(2017,1,1,0,0,0),datetime.datetime(2022,12,3,0,0,0)], #SLSTR on S3A
    [datetime.datetime(2019,1,1,0,0,0),datetime.datetime(2022,12,3,0,0,0)] # SLSTR on S3B
    ])

bad = np.array([
    [], #AVHRR 7
    [], #AVHRR 9
    [], #AVHRR 11
    [], #AVHRR 12
    [], #AVHRR 14
    [], #AVHRR 15
    [], #AVHRR 16
    [], #AVHRR 17
    [], #AVHRR 18
    [], #AVHRR 19
    [datetime.datetime(2007,3,12,0,0,0),datetime.datetime(2007,4,21,0,0,0),datetime.datetime(2007,4,22,0,0,0),datetime.datetime(2007,4,23,0,0,0),datetime.datetime(2007,5,4,0,0,0)
    ,datetime.datetime(2007,5,5,0,0,0),datetime.datetime(2007,9,18,0,0,0),datetime.datetime(2008,3,20,0,0,0),datetime.datetime(2008,4,9,0,0,0),datetime.datetime(2019,8,4,0,0,0)
    ,datetime.datetime(2019,8,3,0,0,0),datetime.datetime(2019,8,2,0,0,0),datetime.datetime(2019,8,1,0,0,0),datetime.datetime(2019,7,31,0,0,0)], #AVHRR on METOPA
    [], #AVHRR on METOPB
    [], #ATSR on ERS1
    [], #ATSR on ERS2
    [], #AASTR on envisat
    [datetime.datetime(2018,2,16,0,0,0),datetime.datetime(2018,2,17,0,0,0),datetime.datetime(2018,2,18,0,0,0),datetime.datetime(2018,2,19,0,0,0),datetime.datetime(2017,2,16,0,0,0)
        ,datetime.datetime(2018,2,20,0,0,0),datetime.datetime(2017,2,18,0,0,0),datetime.datetime(2017,2,17,0,0,0),datetime.datetime(2017,2,15,0,0,0),datetime.datetime(2017,7,31,0,0,0)
        ,datetime.datetime(2017,8,1,0,0,0),datetime.datetime(2017,8,2,0,0,0),datetime.datetime(2017,8,3,0,0,0),datetime.datetime(2017,8,4,0,0,0),datetime.datetime(2017,8,5,0,0,0)
        ,datetime.datetime(2018,9,14,0,0,0),datetime.datetime(2018,9,15,0,0,0),datetime.datetime(2018,9,16,0,0,0),datetime.datetime(2018,9,17,0,0,0),datetime.datetime(2018,9,18,0,0,0)
        ,datetime.datetime(2021,6,30,0,0,0),datetime.datetime(2022,8,7,0,0,0),datetime.datetime(2022,7,22,0,0,0)], #SLSTR on S3A
    [datetime.datetime(2019,3,1,0,0,0),datetime.datetime(2019,3,2,0,0,0),datetime.datetime(2019,3,3,0,0,0),datetime.datetime(2019,3,4,0,0,0),datetime.datetime(2019,3,5,0,0,0)
    ,datetime.datetime(2019,3,6,0,0,0),datetime.datetime(2019,3,7,0,0,0),datetime.datetime(2019,3,8,0,0,0),datetime.datetime(2019,3,9,0,0,0),datetime.datetime(2019,3,10,0,0,0)
    ,datetime.datetime(2019,3,12,0,0,0),datetime.datetime(2019,3,13,0,0,0),datetime.datetime(2019,3,14,0,0,0),datetime.datetime(2019,3,15,0,0,0),datetime.datetime(2019,3,16,0,0,0)
    ,datetime.datetime(2019,3,17,0,0,0),datetime.datetime(2019,3,18,0,0,0),datetime.datetime(2019,2,24,0,0,0),datetime.datetime(2019,2,25,0,0,0),datetime.datetime(2019,2,26,0,0,0)
    ,datetime.datetime(2019,2,27,0,0,0),datetime.datetime(2019,2,27,0,0,0),datetime.datetime(2020,11,12,0,0,0),datetime.datetime(2019,9,24,0,0,0),datetime.datetime(2019,9,23,0,0,0)
    ,datetime.datetime(2019,9,22,0,0,0),datetime.datetime(2019,9,21,0,0,0),datetime.datetime(2019,9,20,0,0,0),datetime.datetime(2019,4,16,0,0,0),datetime.datetime(2019,4,15,0,0,0)
    ,datetime.datetime(2019,4,14,0,0,0),datetime.datetime(2019,4,13,0,0,0),datetime.datetime(2019,4,12,0,0,0),datetime.datetime(2019,2,24,0,0,0),datetime.datetime(2019,2,25,0,0,0)
    ,datetime.datetime(2019,2,26,0,0,0),datetime.datetime(2019,2,27,0,0,0),datetime.datetime(2019,2,28,0,0,0),datetime.datetime(2021,6,30,0,0,0)] # SLSTR on S3B
    ])

def makefolder(folder):
    if os.path.exists(folder) == 0:
        os.mkdir(folder)
        print('Made folder: ' + folder)

makefolder(data_save_loc)

def check_sats(date,sat_dates=sat_dates,sat_names=sats):
    out=[]
    for i in range(0,len(sat_dates)):
        if (date <= sat_dates[i,1]) & (date >= sat_dates[i,0]):
            out = np.append(out,i)
    out = out.astype(np.int64)
    out = np.array(sat_names[out])
    return out

def cci_sst_request(date,sat,save_loc = data_save_loc,sat_t = sats, sat_dates = sat_dates,bads = bad):
    # Function to send a API request to the Copernicus Climate Data Store to retrieve a zip file
    # contain the L3C data for each satellite, or the specified satellite.
    f = np.where(sat == sats)
    if (date <= sat_dates[f,1]) & (date >= sat_dates[f,0]):
        date_d = datetime.datetime(date.year,date.month,date.day)
        # print(date_d)
        # print(bad[f][0])
        if any(date_d == x for x in bad[f][0]):
            print('Bad Date')
            return 0
        else:
            makefolder(os.path.join(save_loc,sat))
            makefolder(os.path.join(save_loc,sat,str(date.year)))
            check = check_satdata_exist(date,sat)
            if check != 2:
                print(sat)
                print(date)
                c = cdsapi.Client()
                c.retrieve('satellite-sea-surface-temperature',
                 {
                   "processinglevel": "level_3c",
                   "version": '2_1',
                   "sensor_on_satellite": sat,
                   "year": date.strftime("%Y"),
                   "month": date.strftime("%m"),
                   "day": date.strftime("%d"),
                   "variable": "all",
                   "format": "zip",
                 },
                os.path.join(save_loc,'output.zip'))
                with zipfile.ZipFile(os.path.join(save_loc,'output.zip'), 'r') as zip_ref:
                    zip_ref.extractall(os.path.join(save_loc,sat,str(date.year)))
            else:
                print('Files exist!')
            files = glob.glob(os.path.join(save_loc,sat,str(date.year),date.strftime("%Y%m%d*.nc")))
            return files
    else:
        print('CCI_Request nope')
        return 0

def check_satdata_exist(date,sat,save_loc = data_save_loc):
    fold = os.path.join(save_loc,sat,str(date.year),date.strftime("%Y%m%d*.nc"))
    files = glob.glob(fold)
    return len(files)

def convert_ins_time(date,ref = datetime.datetime(1981,1,1,0,0,0)):
    return (date-ref).total_seconds()

def convert_to_datetime(dates,d = datetime.datetime(1970,1,1,0,0,0)):
    out = []
    for i in range(len(dates)):
        #print(datetime.timedelta(seconds=int(dates[i])))
        #print(dates[i])
        out.append(d + datetime.timedelta(seconds=float(dates[i])))
    return np.array(out)

def load_cci_sst(file,i_time,i_lat,i_lon,timewin = 1,qual_le = 3,plot = 0 ):
    # Function to load the Tskin from the retrieved file and makes sure the delta time between in situ observation
    # and the satellite observation is less than the time window in hours
    pad = 40
    c = Dataset(file)
    s_ref = np.array(c.variables['time'])
    s_lat = np.array(c.variables['lat'])
    s_lon = np.array(c.variables['lon'])
    f,g = calc_nearest(i_lat,i_lon,s_lat,s_lon)
    #print(f,g)
    s_time = np.array(c.variables['sst_dtime'][0,f,g])+s_ref
    print(np.abs(i_time-s_time))
    if np.abs(i_time-s_time) < timewin*3600:
        sst_c = np.array(c.variables['sea_surface_temperature'][0,f,g])
        unc_c = np.array(c.variables['sea_surface_temperature_total_uncertainty'][0,f,g])
        qual_c = np.array(c.variables['quality_level'][0,f,g])
        if qual_c <= qual_le:
            sst_c = np.nan
            unc_c = np.nan
            qual_c = np.nan
            s_time = np.nan
        if plot == 1:
            f2 = np.arange(f-pad,f+pad,dtype=np.int64)
            g2 = np.arange(g-pad,g+pad,dtype=np.int64)
            sst_r = np.array(c.variables['sea_surface_temperature'][0,f2,g2])
            qual_r = np.array(c.variables['quality_level'][0,f2,g2])
            sst_r[qual_r <=qual_le] = np.nan
            sst_r[sst_r <= -200] = np.nan
            plt.figure()
            plt.pcolor(s_lon[g2],s_lat[f2],sst_r)
            plt.colorbar()
            plt.scatter(i_lon,i_lat)
            plt.show()
        return [sst_c, unc_c, qual_c, s_time, s_lat[f], s_lon[g]]
    else:
        print('Outside Time Range!')
        return [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
    c.close()

def calc_nearest(i_lat,i_lon,s_lat,s_lon):
    lat_d = np.abs(i_lat - s_lat)
    lon_d = np.abs(i_lon - s_lon)

    f = np.squeeze(np.argwhere(lat_d == np.min(lat_d)))
    g = np.squeeze(np.argwhere(lon_d == np.min(lon_d)))
    #print(f.size)
    #print(g.size)
    if (f.size > 1) | (g.size > 1):
        print('Help')
        lat_g,lon_g = np.meshgrid(lat_d[f],lon_d[g])
        f_g,g_g = np.meshgrid(f,g)
        f_g,g_g = (f_g).flatten(),(g_g).flatten()
        t = (np.sqrt(lat_g**2 + lon_g**2)).flatten()
        # print(t)
        # print(f_g)
        # print(g_g)
        p = np.squeeze(np.argwhere(np.min(t) == t))
        f = f_g[p]
        g = g_g[p]
        if (f.size > 1) | (g.size > 1):
            f = f[0]
            g = g[0]
    return f,g

def find_all_sats(time):
    sats_all = []
    for i in range(int(len(time))):
        date = time[i]
        sats = CCI.check_sats(date)
        print(sats)
        for k in sats:
            if k not in sats_all:
                sats_all.append(k)
    return sats_all
