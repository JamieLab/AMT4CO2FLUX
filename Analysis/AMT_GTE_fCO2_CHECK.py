#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import weight_stats as ws
import datetime
import matplotlib.transforms
font = {'weight' : 'normal',
        'size'   : 16}
matplotlib.rc('font', **font)
matplotlib.rcParams['text.usetex'] = True

def find_lat(time_s,time_u,lat):
    g = np.abs(time_s - time_u)
    g = np.argwhere(g == np.nanmin(g))
    return lat[g[0]].squeeze()

infile = 'AMT29/fCO2fromTADIC.txt'
min20_file = 'AMT29/DATA/AMT29_table_20min_v2_20221022.txt'

out = pd.read_table(infile)
out_time = []
for i in range(0,len(out)):
    out_time.append(datetime.datetime.strptime(out['Timestamp'][i],'%Y-%m-%d %H:%M:%S'))
out_time = np.array(out_time)
print(out)
min20 = pd.read_table(min20_file,skiprows=28)
min20_time = []
for i in range(0,len(min20)):
    min20_time.append(datetime.datetime.strptime(min20['YearMonthDay_HHMMSS_Mid'][i],'%Y-%m-%d %H:%M:%S'))
min20_time = np.array(min20_time)

#Extracting the rough latitude where the DIC/TA measurements were made for comparision.
#Are we seeing larger bias' in the South Atlantic.
lat = np.empty((len(out)))
for i in range(0,len(out)):
    lat[i]=find_lat(out_time[i],min20_time,np.array(min20['lat']))
print(lat)

#Plotting

cols = ['#332288','#44AA99','#DDCC77', '#117733', '#88CCEE','#999933','#CC6677']
lim = np.array([330,450])
fig1 = plt.figure(figsize=(15,7))
gs = GridSpec(1,2, figure=fig1, wspace=0.33,hspace=0.2,bottom=0.1,top=0.95,left=0.15,right=0.9)
ax1 = fig1.add_subplot(gs[0,0])
f = out['Depth'] == 2
ax1.scatter(out['fCO2w_calcTADIC'][f],out['fCO2w_meas'][f],color=cols[0],label='2m')
ax1.errorbar(out['fCO2w_calcTADIC'][f],out['fCO2w_meas'][f],yerr = 2.6,xerr=10,linestyle='',color=cols[0])
f = out['Depth'] == 5
w = ws.unweighted_stats(np.array(out['fCO2w_calcTADIC'])[f],np.array(out['fCO2w_meas'])[f],'TA_DIC')
w2 = ws.unweighted_stats(np.array(out['fCO2w_calcTADIC']),np.array(out['fCO2w_meas']),'TA_DIC')
ax1.scatter(out['fCO2w_calcTADIC'][f],out['fCO2w_meas'][f],color=cols[1],label='5m')
ax1.errorbar(out['fCO2w_calcTADIC'][f],out['fCO2w_meas'][f],yerr = 2.6,xerr=10,linestyle='',color=cols[1])
ax1.plot(lim,lim,'k-')
ax1.plot(lim,lim*w['slope'] + w['intercept'],linestyle='--',color=cols[1])
ax1.plot(lim,lim*w2['slope'] + w2['intercept'],linestyle='--',color='k')
ax1.set_xlabel('DIC/TA derived fCO$_{2 (sw,depth)}$ ($\mu$atm)')
ax1.set_ylabel('SFCE estimated fCO$_{2 (sw,depth)}$ ($\mu$atm)')
ax1.grid()
ax1.set_xlim(lim)
ax1.set_ylim(lim)
ax1.text(0.32,0.18,'All Data (n = ' + str(w2['n']) + ') \n RMSD = '+str(round(w2['rmsd'],1)) + ' $\mu$atm; Mean Bias = ' + str(round(w2['rel_bias'],2)) + ' $\mu$atm \n \n 5m data (n = ' + str(w['n']) + ') \n RMSD = '+str(round(w['rmsd'],2)) + ' $\mu$atm; Mean Bias = ' + str(round(w['rel_bias'],2)) + ' $\mu$atm' ,transform=ax1.transAxes,va='top',fontsize=12)
print(w)
print(w2)

f = out['Depth'] == 10
ax1.scatter(out['fCO2w_calcTADIC'][f],out['fCO2w_meas'][f],color=cols[2],label = '10m')
ax1.errorbar(out['fCO2w_calcTADIC'][f],out['fCO2w_meas'][f],yerr = 2.6,xerr=10,linestyle='',color=cols[2])
ax1.legend(loc = 3)
ax2 = fig1.add_subplot(gs[0,1])
dep = [2,5,10]
for i in range(0,len(dep)):
    f = out['Depth'] == dep[i]
    ax2.scatter(lat[f],out['fCO2w_calcTADIC'][f]-out['fCO2w_meas'][f],label=str(dep[i])+'m',color=cols[i])
    ax2.errorbar(lat[f],out['fCO2w_calcTADIC'][f]-out['fCO2w_meas'][f],yerr=np.sqrt(2.6**2 + 10**2),linestyle='',color=cols[i])
ax2.set_ylabel('DIC/TA derived fCO$_{2 (sw,depth)}$ - SFCE estimated fCO$_{2 (sw,depth)}$ ($\mu$atm)')
ax2.set_xlabel('Latitude ($^{o}$N)')
ax2.grid()
ax2.legend()
ax2.invert_xaxis()
ax2.text(0.03,0.95,r"{\fontsize{20}{22}\textbf{(b)}}",transform=ax2.transAxes,va='top')
ax1.text(0.03,0.95,r"{\fontsize{20}{22}\textbf{(a)}}",transform=ax1.transAxes,va='top')
fig1.savefig('FIGS/SUP_FIG_11_GTE_COMPARISION.png',format='png',dpi=300)
