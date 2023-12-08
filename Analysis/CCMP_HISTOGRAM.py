#!/usr/bin/env python3
import glob
from netCDF4 import Dataset
import numpy as np
import os
import datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def ccmp_load(loc,startyr,endyr):
    tot = ((endyr - startyr)+1) * 12
    yr = startyr
    mon = 1
    co = 0
    i = 0
    while co == 0:
        print(str(yr) + '/' + str(mon))
        filepath = file_loc_construct(loc,yr,mon)
        files = glob.glob(filepath)
        c = Dataset(files[0])
        if i == 0:
            lat = np.array(c.variables['latitude'])
            lon = np.array(c.variables['longitude'])
            ws = np.empty((len(lat),len(lon),tot))
            land_mask = load_land(lat)
            #print(land_mask.shape)
        wspd = np.squeeze(np.array(c.variables['w']))
        wspd[land_mask != 0] = np.nan
        # plt.figure()
        # plt.pcolor(lon,lat,wspd)
        # plt.show()
        #print(wspd.shape)
        ws[:,:,i] = wspd
        i = i+1
        mon = mon+1
        print(mon)
        print(yr)
        if mon == 13:
            yr = yr+1
            mon=1
        if yr == endyr+1:
            print('c')
            co = 1
        print(co)
    return ws,lat,lon

def file_loc_construct(loc,yr,mon):
    d = datetime.datetime(yr,mon,1)
    filepath = os.path.join(loc,'Y'+str(d.year),'M'+d.strftime('%m'),'*L4.0.nc')
    return filepath

def load_land(lat,loc = 'D:/Data/CCMP/land_sea_mask_0.25.nc'):
    c = Dataset(loc)
    lsm = np.squeeze(np.array(c.variables['lsm']))
    latg = np.array(c.variables['latitude']) + 0.125

    f = np.where((latg >= lat[0]) & (latg <= lat[-1]))
    lsm = lsm[f,:]
    lsm = np.flip(lsm,axis=1)
    return np.squeeze(lsm)

def plot_histogram(ax,data):
    #fig,ax = plt.subplots()
    data = data.reshape(-1,1)
    f = np.where((np.isnan(data) == 0))
    ax.hist(data[f],bins = np.arange(0,17.5,.5),edgecolor='k',label='CCMP',weights = np.ones_like(data[f])/float(len(data[f])))
    #ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=len(data)))
    #plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_ylabel('Distribution (\% per 0.5 ms$^{-1}$)')
    ax.set_xlabel('U$_{10}$ (ms$^{-1}$)')
    #plt.show()

def plot_ccmp(ax):
    #ccmp_load('D:/Data/CCMP/DAILY',2002,2018)
    loc = 'D:/Data/CCMP/DAILY'
    plot_histogram(ax,ccmp_load(loc,2002,2018))
