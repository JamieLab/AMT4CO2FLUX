#!/usr/bin/env python3

import CCISST_DF_Functions as CCI
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

in_filenames = ['AMT28_table_20min_v3_20221205_S3A_ADDED.txt','AMT29_table_20min_v3_20221205_S3A_ADDED.txt']
for in_filename in in_filenames:
    data = pd.read_table(in_filename,sep='\t')#,skiprows=26)
    print(data)

    time = []
    #time[:] = np.nan
    for i in range(len(data)):
        time.append(datetime.datetime.strptime(data['YearMonthDay_HHMMSS_Mid'][i],'%d/%m/%Y %H:%M'))
    time = np.array(time)


    sats_all = ['avhrr_on_noaa_19', 'avhrr_on_metop_a','slstr_on_sentinel_3a','slstr_on_sentinel_3b']

    for k in range(len(sats_all)):
        out = []
        for i in range(0,len(data)):
            date = time[i]
            files = CCI.cci_sst_request(date,sats_all[k])
            if files == 0:
                print('Outside Time')
                out.append(np.nan)
                out.append(np.nan)
                out.append(np.nan)
                out.append(np.nan)
            else:
                for p in range(len(files)):
                    print(files[p])
                    #print(CCI.load_cci_sst(files[p],CCI.convert_ins_time(date),data['lat'][i],data['lon'][i],timewin=1))
                    put = CCI.load_cci_sst(files[p],CCI.convert_ins_time(date),data['lat'][i],data['lon'][i],timewin=2)
                    out.append(put[0])
                    out.append(put[1])
        out=np.array(out)
        out=out.reshape((-1,4))
        sst = np.nanmean(out[:,[0,2]],axis=1)
        unc = np.sqrt(np.nansum(out[:,[1,3]]**2,axis=1)) / np.sum(np.isnan(out[:,[1,3]]) == 0,axis=1)
        data['CCI_'+sats_all[k]+'_mean'] = sst
        data['CCI_'+sats_all[k]+'_unc'] = unc

    sst_m = np.nanmean(np.array(data[('CCI_' + s + '_mean' for s in sats_all)]),axis=1)
    unc_m = np.array(data[('CCI_' + s + '_unc' for s in sats_all)])
    unc_m = np.sqrt(np.nansum(unc_m**2,axis=1)) / np.sum(np.isnan(unc_m) == 0,axis=1)
    data['CCI_merged_mean'] = sst_m
    data['CCI_merged_unc'] = unc_m
    sst_m = np.nanmean(np.array(data[('CCI_' + s + '_mean' for s in sats_all[-2:])]),axis=1)
    unc_m = np.array(data[('CCI_' + s + '_unc' for s in sats_all[-2:])])
    unc_m = np.sqrt(np.nansum(unc_m**2,axis=1)) / np.sum(np.isnan(unc_m) == 0,axis=1)
    data['CCI_s3a+b_mean'] = sst_m
    data['CCI_s3a+b_unc'] = unc_m
    outfile = in_filename.replace('.txt', '_CCISST_ADDED.txt')
    data.to_csv(outfile,'\t',index=None,na_rep='NaN')
