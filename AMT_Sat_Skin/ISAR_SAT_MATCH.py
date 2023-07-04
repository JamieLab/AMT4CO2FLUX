#!/usr/bin/env python3
import datetime
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import glob
import os
import weight_stats as ws
import csv

def load_s3a(c):
    lat_sat = np.array(c.variables['s3a_sl_2_wst____o_nr__lat'])
    lon_sat = np.array(c.variables['s3a_sl_2_wst____o_nr__lon'])
    sst_sat = np.array(c.variables['s3a_sl_2_wst____o_nr__sea_surface_temperature'])
    qual_sat =np.array(c.variables['s3a_sl_2_wst____o_nr__quality_level'])
    ind_sat = np.array(c.variables['s3a_sl_2_wst____o_nr__dynamic_target_center_index'])
    time_sat = np.array(c.variables['s3a_sl_2_wst____o_nr__time'])
    unc_sat = np.array(c.variables['s3a_sl_2_wst____o_nr__sst_theoretical_uncertainty'])

    return lat_sat,lon_sat,sst_sat,qual_sat,ind_sat,time_sat,unc_sat

def load_s3b(c):
    lat_sat = np.array(c.variables['s3b_sl_2_wst____o_nr__lat'])
    lon_sat = np.array(c.variables['s3b_sl_2_wst____o_nr__lon'])
    sst_sat = np.array(c.variables['s3b_sl_2_wst____o_nr__sea_surface_temperature'])
    qual_sat =np.array(c.variables['s3b_sl_2_wst____o_nr__quality_level'])
    ind_sat = np.array(c.variables['s3b_sl_2_wst____o_nr__dynamic_target_center_index'])
    time_sat = np.array(c.variables['s3b_sl_2_wst____o_nr__time'])
    unc_sat = np.array(c.variables['s3b_sl_2_wst____o_nr__sst_theoretical_uncertainty'])

    return lat_sat,lon_sat,sst_sat,qual_sat,ind_sat,time_sat,unc_sat

def time_correct(dates,d = datetime.datetime(1970,1,1,0,0,0)):
    out = []
    for i in range(len(dates)):
        #print(datetime.timedelta(seconds=int(dates[i])))
        #print(dates[i])
        out.append(d + datetime.timedelta(seconds=float(dates[i])))
    return np.array(out)

loc = ['AMT28/DATA/ISAR_UPDATED/2018/sat','AMT29/DATA/ISAR_UPDATED/2019/sat']
output = []
for p in loc:
    files = glob.glob(os.path.join(p,'s3a*.nc'))
    #i=4
    for i in range(0,len(files)):
        print(files[i])
        st = files[i].split('\\')
        c = Dataset(files[i])
        site = np.array(c.variables['felyx_site_name'])
        print(site)
        site_t = len(np.argwhere((site == 'JCR') | (site == 'DCY')))
        if site_t > 0:
            lat_isar = np.array(c.variables['ship4sstr1i1__lat'])
            lon_isar = np.array(c.variables['ship4sstr1i1__lon'])
            ind_isar = np.array(c.variables['ship4sstr1i1__dynamic_target_center_index'])
            time_isar = np.array(c.variables['ship4sstr1i1__time']).astype(np.float32)
            sst_isar = np.array(c.variables['ship4sstr1i1__sea_surface_temperature'])
            unc_isar = np.array(c.variables['ship4sstr1i1__sst_total_uncertainty'])
            #print(sst_isar.shape)
            time_isar[time_isar<0] = np.nan
            lat_isar[lat_isar>90] = np.nan
            lon_isar[lon_isar>90] = np.nan
            if st[-1][2] == 'a':
                lat_sat,lon_sat,sst_sat,qual_sat,ind_sat,time_sat,unc_sat = load_s3a(c)
            else:
                lat_sat,lon_sat,sst_sat,qual_sat,ind_sat,time_sat,unc_sat = load_s3b(c)
            sh = sst_sat.shape
            # print(ind_sat)
            #print(sh)
            # print(time_sat)
            # print(time_isar)
            sst_sat[qual_sat <= 3] = np.nan
            unc_sat[qual_sat <= 3] = np.nan
            c.close()
            #fig,ax = plt.subplots(int(np.ceil(site_t/2)),2)
            #ax = ax.flatten()
            #print(ax)
            #ax = [ax1,ax2,ax3,ax4]
            dt = 7200
            pl = 0
            for j in range(0,sh[0]):
                if site[j] == 'JCR' or site[j] == 'DCY':
                    #print(time_correct(time_isar[j,:],datetime.datetime(1981,1,1,0,0,0)))
                    #print(lat_isar[j,:])
                    # m = ax[pl].pcolor(lon_sat[j,:,:],lat_sat[j,:,:],sst_sat[j,:,:])
                    # plt.colorbar(m,ax=ax[pl])
                    # ax[pl].plot(lon_isar[j,:],lat_isar[j,:])
                    t = np.nanmean(time_sat[j]) - time_isar[j,:]
                    f = np.argwhere((t <=dt) & (t >=-dt))
                    # ax[pl].plot(lon_isar[j,f],lat_isar[j,f])
                    # pl = pl+1
                    for k in f:
                        dis = np.sqrt((lon_isar[j,k] - lon_sat[j,:,:])**2 + (lat_isar[j,k] - lat_sat[j,:,:])**2)
                        #print(np.min(dis))
                        g = np.argwhere(dis == np.min(dis))[0]
                        print(g)
                        #print(np.array([sst_isar[j,k][0],sst_sat[j,g[0],g[1]],t[k]]).shape)
                        #print(output.shape)
                        output = np.append(output,np.array([sst_isar[j,k][0],sst_sat[j,g[0],g[1]],t[k][0],unc_isar[j,k][0],unc_sat[j,g[0],g[1]],time_isar[j,k][0],lat_isar[j,k][0],lon_isar[j,k][0]]))
                #ax[j].scatter(lon_isar[j,ind_isar[j]],lat_isar[j,ind_isar[j]])
            # print(np.reshape(output,[-1,2]))
            # print(output)
            #plt.show()
output = np.reshape(output,[-1,8])
print(output)
f,(ax1,ax2) = plt.subplots(1,2)
m = ax1.scatter(output[:,0],output[:,1],c = output[:,2]/3600,zorder=3)
ax1.errorbar(output[:,0],output[:,1],yerr=output[:,4],xerr=output[:,3],linestyle='',zorder=1)
cbar = plt.colorbar(m,ax=ax1)
cbar.set_label('Time$_{satellite}$ - Time$_{in situ}$ (hours)')
g = np.argwhere((np.isnan(output[:,0]) == 0) & (np.isnan(output[:,1]) == 0))
out = []
stats=ws.unweighted_stats(np.squeeze(output[:,0]),np.squeeze(output[:,1]),'S3A_2hr')
out.append(stats)
stats=ws.weighted_stats(np.squeeze(output[:,0]),np.squeeze(output[:,1]), np.sqrt(np.squeeze(output[:,3])**2 + np.squeeze(output[:,4])**2),'S3A_2hr_weighted')
out.append(stats)
lim = ax1.get_xlim()
ax1.plot(lim,lim,'k--')
ax1.plot(lim,(np.array(lim)*stats['slope'])+stats['intercept'],'r--')
ax1.set_xlabel('in situ T$_{skin}$ (K)')
ax1.set_ylabel('Satellite T$_{skin}$ (K)')
ax1.grid()
p = np.squeeze(np.argwhere((output[:,2] <= 3600) & (output[:,2] >= -3600)))
#print(p)
output_1h = output[p,:]
m = ax2.scatter(output_1h[:,0],output_1h[:,1],c = output_1h[:,2]/3600,zorder=3)
ax2.errorbar(output_1h[:,0],output_1h[:,1],yerr=output_1h[:,4],xerr=output_1h[:,3],linestyle = '',zorder=1)
cbar = plt.colorbar(m,ax=ax2)
cbar.set_label('Time$_{satellite}$ - Time$_{in situ}$ (hours)')
g = np.argwhere((np.isnan(output_1h[:,0]) == 0) & (np.isnan(output_1h[:,1]) == 0))
stats2=ws.unweighted_stats(np.squeeze(output_1h[:,0]),np.squeeze(output_1h[:,1]),'S3A_1hr')
out.append(stats2)
stats2=ws.weighted_stats(np.squeeze(output_1h[:,0]),np.squeeze(output_1h[:,1]), np.sqrt(np.squeeze(output_1h[:,3])**2 + np.squeeze(output_1h[:,4])**2),'S3A_1hr_weighted')
out.append(stats2)
lim = ax2.get_xlim()
ax2.plot(lim,lim,'k--')
ax2.plot(lim,(np.array(lim)*stats2['slope'])+stats2['intercept'],'r--')
ax2.set_xlabel('in situ T$_{skin}$ (K)')
ax2.set_ylabel('Satellite T$_{skin}$ (K)')
ax2.grid()
#print(output_1h)
field_names = out[0].keys()
with open('SATELLITE/3A_STATS.csv','w',newline='') as l:
    writer = csv.DictWriter(l, fieldnames = field_names)
    writer.writeheader()
    writer.writerows(out)
np.savetxt("SATELLITE/2h_ISAR_S3A_MATCH.csv", output, delimiter=",")
np.savetxt("SATELLITE/1h_ISAR_S3A_MATCH.csv", output_1h, delimiter=",")
plt.show()
