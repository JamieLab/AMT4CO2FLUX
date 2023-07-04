#!/usr/bin/env python3
import pandas as pd
import datetime
import numpy as np
import sys
import scipy.stats as sk
def time_correct(dates,d = datetime.datetime(1970,1,1,0,0,0)):
    out = []
    for i in range(len(dates)):
        #print(datetime.timedelta(seconds=int(dates[i])))
        #print(dates[i])
        out.append(d + datetime.timedelta(seconds=float(dates[i])))
    return np.array(out)

def unc_prop(sst_u):
    g = np.sum(np.isnan(sst_u) == 0)
    out = np.nansum(sst_u**2)/g
    return out

ISAR = ['SATELLITE/2h_ISAR_S3A_MATCH.csv']
#print(np.array(cruise_time))
ISAR = pd.read_table(ISAR[0],sep=',',header=None)

ISAR_T = time_correct(ISAR[5],d = datetime.datetime(1981,1,1,0,0,0))
ISAR[5] = ISAR_T
#j=0
file20 = ['AMT28\DATA\AMT28_table_20min_v3_20221205.txt','AMT29\DATA\AMT29_table_20min_v3_20221205.txt']
for j in range(0,len(file20)):
    data = pd.read_table(file20[j],sep='\t',skiprows=26)
    cruise_time = []
    for i in range(0,len(data)):
        cruise_time.append(datetime.datetime.strptime(np.array(data['YearMonthDay_HHMMSS_Mid'])[i],'%d/%m/%Y %H:%M'))
    cruise_time = np.array(cruise_time)
    #print(data)
    #print(ISAR)
    out = np.empty((len(data),5))
    out[:] = np.nan
    for i in range(0,len(data)):
        f = np.squeeze(np.argwhere( (ISAR_T < cruise_time[i] + datetime.timedelta(minutes=10) ) & (ISAR_T > cruise_time[i] - datetime.timedelta(minutes=10) )))
        #print(f.size)
        if f.size!=0:
            sst = np.array(ISAR[1][f])
            sst_u = np.array(ISAR[3][f])
            print(sst_u)
            #print(np.array(ISAR[0][f]))
            print(data['lat'][i])
            print(data['CO2Flux_EC'][i])
            tout = [np.nanmean(sst),np.nanmedian(sst),np.nanstd(sst),sk.skew(sst,nan_policy='omit'),unc_prop(sst_u)]
            print(tout)
            out[i,:]=tout

    out[np.isnan(out[:,0]),:] = np.nan
    data['S3A_Tskin_mean'] = out[:,0]
    data['S3A_Tskin_med'] = out[:,1]
    data['S3A_Tskin_std'] = out[:,2]
    data['S3A_Tskin_skew'] = out[:,3]
    data['S3A_Tskin_unc'] = out[:,4]
    print(data.columns.values)
    outfile = file20[j].replace('.txt', '_S3A_ADDED.txt')
    data.to_csv(outfile,'\t',index=None,na_rep='NaN')
