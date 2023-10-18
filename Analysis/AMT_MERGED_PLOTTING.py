#!/usr/bin/env python3
# ---------------------------------
# Code created by Daniel Ford (d.ford@exeter.ac.uk) - 18/08/2022
# Code designed to plot the fluxengine output and Eddy covariance CO2 fluxes for comparision
#-------- Version 1.1 - 05/10/2022
# - Changed 3 hour binning to a static time window following Tom Bell
# - Changed data cutoff for 3 hour binning from 1/2 to 2/3 following Tom Bell
# -
#--------  Version 1.0
# - Initial Version
# ---------------------------------

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from netCDF4 import Dataset
import numpy as np
import numpy.matlib
from pylr2 import regress2
import matplotlib.transforms
from collections import defaultdict
import csv
import weight_stats as ws
import datetime
from random import normalvariate
import CCMP_HISTOGRAM as ccmp
from scipy.stats import mannwhitneyu
from scipy.stats import f_oneway
import AMT_mapo as mapo

font = {'weight' : 'normal',
        'size'   : 19}
matplotlib.rc('font', **font)
matplotlib.rcParams['text.usetex'] = True

def data_load_feout(file,time_i,window,cutoff = 0.5):
    # Loadinng the variables from the Fluxengin output netcdf. These are saved raw
    # in a structure, and then two 3 hour averaging approaches are conducted (one static times, one flexible 9 point averages)
    c = Dataset(file,'r')
    EC_flux = np.array(c.variables['CO2Flux_EC'])
    EC_flux_u = np.array(c.variables['Err_CO2flux_EC'])
    lat = np.array(c.variables['lat'])
    lon = np.array(c.variables['lon'])
    vt = np.array(c.variables['vt'])
    vts = np.array(c.variables['vts'])
    nvt = np.array(c.variables['nvt'])
    nvts = np.array(c.variables['nvts'])
    Tskin = np.array(c.variables['Tskin_mean'])
    Tsst = np.array(c.variables['SST_mean'])
    time = time_correct(np.array(c.variables['time']))
    tsubskin_co = np.array(c.variables['T_subskin_co'])
    tsubskin_don = np.array(c.variables['T_subskin_don'])
    fco2_atm = np.array(c.variables['fCO2_atm_mean'])
    fco2_sw = np.array(c.variables['fCO2_sw_mean'])
    Tskin_unc = np.array(c.variables['Tskin_uncertainty'])
    ws = np.array(c.variables['U10n_mean'])
    c.close()
    # Saving the raw data into a structure
    struct_m = {
        'EC_flux': EC_flux,
        'EC_flux_u': EC_flux_u,
        'latitude': lat,
        'longitude': lon,
        'vt': vt,
        'vts':vts,
        'nvt':nvt,
        'nvts':nvts,
        'Tskin': Tskin,
        'Tsst':Tsst,
        'Tsubskin_co':tsubskin_co,
        'Tsubskin_don':tsubskin_don,
        'fco2_sw':fco2_sw,
        'fco2_atm':fco2_atm,
        'time':time,
        'Tskin_unc':Tskin_unc,
        'U10n':ws
    }
    # Saving structure with data averaged to a fixed 3 hour time window
    # Currently data is averaged onto a time grid supplied (i.e from Tom Bell's data)
    E_m,E_se = mean_bin_fix(EC_flux,time,time_i,window,cutoff)
    vt_m,vt_se = mean_bin_fix(vt,time,time_i,window,cutoff)
    nvt_m,nvt_se =mean_bin_fix(nvt,time,time_i,window,cutoff)
    struct = {
        'EC_flux': E_m,
        'EC_SE': E_se,
        'EC_flux_u': unc_prop_fix(EC_flux_u,time,time_i,window,cutoff),
        'latitude': mean_bin_fix(lat,time,time_i,window,cutoff)[0],
        'longitude': mean_bin_fix(lon,time,time_i,window,cutoff)[0],
        'vt': vt_m,
        'vtse': vt_se,
        'vts':unc_prop_fix(vts,time,time_i,window,cutoff),
        'nvt': nvt_m,
        'nvtse':nvt_se,
        'nvts':unc_prop_fix(nvts,time,time_i,window,cutoff),
        'fco2_atm':mean_bin_fix(fco2_atm,time,time_i,window,cutoff)[0],
        'fco2_sw':mean_bin_fix(fco2_sw,time,time_i,window,cutoff)[0],
        'Tskin': mean_bin_fix(Tskin,time,time_i,window,cutoff)[0],
        'Tsst': mean_bin_fix(Tsst,time,time_i,window,cutoff)[0],
        'ws': mean_bin_fix(ws,time,time_i,window,cutoff)[0]
        }
    # Saving structure with data averaged to a 9 point average through the data.
    struct_t= {
        'EC_flux': mean_bin(EC_flux,9),
        'EC_flux_u': unc_bin(EC_flux_u,9),
        'latitude': mean_bin(lat,9),
        'vt': mean_bin(vt,9),
        'vts':unc_bin(vts,9),
        'nvt':mean_bin(nvt,9),
        'nvts':unc_bin(nvts,9),
        'fco2_atm':mean_bin(fco2_atm,9),
        'fco2_sw':mean_bin(fco2_sw,9),
        'Tskin': mean_bin(Tskin,9),
        'Tsst': mean_bin(Tsst,9)
        }
    # For each row, we removed data where any of the varaibles are nan so all statistical comparisions
    # are conducted on the same data.  (This is for the fixed window data)
    a = np.argwhere(np.isnan(struct['EC_flux'])  | np.isnan(struct['vt'][:,3]) | np.isnan(struct['nvt'][:,1]) | np.isnan(struct['EC_flux_u']))
    # # Cycle through all variables and remove the data where any of the above are nan
    for key in struct:
        #print(struct[key])
        t = struct[key].shape
        if len(t) > 1: # If the array is not just a vector, cycle through each column
            for j in range(t[1]):
                struct[key][a,j] = np.nan
        else: # Data is a vector and must be handled differently
            struct[key][a] = np.nan

    # For each row, we removed data where any of the varaibles are nan so all statistical comparisions
    # are conducted on the same data.  (This is for the 9 point average data)
    a = np.argwhere(np.isnan(struct_t['EC_flux'])  | np.isnan(struct_t['vt'][:,3]) | np.isnan(struct_t['nvt'][:,1]) | np.isnan(struct_t['EC_flux_u']))
    for key in struct_t:
        t = struct_t[key].shape
        if len(t) > 1: # If the array is not just a vector, cycle through each column
            for j in range(t[1]):
                struct_t[key][a,j] = np.nan
        else: # Data is a vector and must be handled differently
            struct_t[key][a] = np.nan
    #print(struct)
    return struct,struct_m,struct_t

def load_3hr(file):
    data = pd.read_table(file,sep=' ',skiprows=20,index_col=False,encoding='utf-8')
    #print(data)
    #print(data['latitude'])
    #print(data['DateTime'])
    time = time_correct(data['DateTime']*1e-9) # *by 1e-9 as these are in nanoseconds so converting to seconds to be consistent with 20 min files.
    return data,time

def time_correct(dates,d = datetime.datetime(1970,1,1,0,0,0)):
    out = []
    for i in range(len(dates)):
        #print(datetime.timedelta(seconds=int(dates[i])))
        #print(dates[i])
        out.append(d + datetime.timedelta(seconds=int(dates[i])))
    return np.array(out)

def mean_bin_fix(data,time,time_i,window,cutoff=0.5):
    #print((window/datetime.timedelta(minutes=20)))
    t = data.shape
    # print(time_i)
    # print(time)
    if len(t) > 1:
        t = t[1]
    else:
        t = 1
    out = np.empty((len(time_i),t))
    se = np.empty((len(time_i),t))
    out[:] = np.nan
    se[:] = np.nan
    if t > 1:
        for i in range(len(time_i)):
            f = np.argwhere((time >= time_i[i]) & (time <= time_i[i]+window))
            for j in range(t):
                #print(np.sum(np.isnan(data[f,j])))
                #print((window/datetime.timedelta(minutes=20))*cutoff)
                if np.sum(np.isnan(data[f,j])) <= (window/datetime.timedelta(minutes=20))*cutoff:
                    #print('Used')
                    #print(data[f,j])
                    out[i,j] = np.nanmean(data[f,j])
                    se[i,j] = np.nanstd(data[f,j])#/np.sqrt((len(f) - np.sum(np.isnan(data[f,j]))))
    else:
        for i in range(len(time_i)):
            f = np.argwhere((time >= time_i[i]) & (time <= time_i[i]+window))
            # print(f)
            if np.sum(np.isnan(data[f])) <= (window/datetime.timedelta(minutes=20))*cutoff:
                # print(np.sum(np.isnan(data[f])))
                #print(data[f])
                out[i] = np.nanmean(data[f])
                se[i] = np.nanstd(data[f])#/np.sqrt((len(f) - np.sum(np.isnan(data[f]))))
    # print(out)
    return np.squeeze(out),np.squeeze(se)

def unc_prop_fix(data,time,time_i,window,cutoff=0.5):
    t = data.shape
    if len(t) > 1:
        t = t[1]
    else:
        t = 1
    out = np.empty((len(time_i),t))
    out[:] = np.nan
    if t > 1:
        for i in range(len(time_i)):
            f = np.argwhere((time >= time_i[i]) & (time <= time_i[i]+window))
            for j in range(t):
                if np.sum(np.isnan(data[f,j])) <= (window/datetime.timedelta(minutes=20))*cutoff:
                    out[i,j] = np.sqrt(np.nansum(data[f,j]**2)) / (len(f) - np.sum(np.isnan(data[f,j])))
    else:
        for i in range(len(time_i)):
            f = np.argwhere((time >= time_i[i]) & (time <= time_i[i]+window))
            if np.sum(np.isnan(data[f])) <= (window/datetime.timedelta(minutes=20))*cutoff:
                out[i] = np.sqrt(np.nansum(data[f]**2)) / (len(f) - np.sum(np.isnan(data[f])))
    return np.squeeze(out)

def mean_bin(data,window):
    t = data.shape
    if len(t) > 1:
        t = t[1]
    else:
        t = 1
    l = range(0,len(data),window)
    # print(l)
    out = np.empty(((len(l)-1),t))
    if t > 1:
        for i in range(len(l)-1):
            for j in range(t):
                if np.sum(np.isnan(data[l[i]:l[i+1],j])) <= 3:
                    out[i,j] = np.nanmean(data[l[i]:l[i+1],j])
                else:
                    out[i,j] = np.nan
    else:
        for i in range(len(l)-1):
            if np.sum(np.isnan(data[l[i]:l[i+1]])) <= 3:
                out[i] = np.nanmean(data[l[i]:l[i+1]])
            else:
                out[i] = np.nan
    return np.squeeze(out)

def unc_bin(data,window):
    t = data.shape
    if len(t) > 1:
        t = t[1]
    else:
        t = 1
    l = range(0,len(data),window)
    # print(l)
    out = np.empty(((len(l)-1),t))

    if t > 1:
        for i in range(len(l)-1):
            for j in range(t):
                if np.sum(np.isnan(data[l[i]:l[i+1],j])) <= 3:
                    out[i,j] = np.sqrt(np.nansum(data[l[i]:l[i+1],j]**2))/np.sqrt(window-np.sum(np.isnan(data[l[i]:l[i+1],j])))
                else:
                    out[i,j] = np.nan
    else:
        for i in range(len(l)-1):
            if np.sum(np.isnan(data[l[i]:l[i+1]])) <= 3:
                out[i] = np.sqrt(np.nansum(data[l[i]:l[i+1]]**2))/np.sqrt(window-np.sum(np.isnan(data[l[i]:l[i+1]])))
            else:
                out[i] = np.nan
    return np.squeeze(out)

def plot_latitude(ax,AMT,let,vals = [7,1],ylim = [-25,15]):
    han = ax.errorbar(AMT['latitude'],AMT['EC_flux'],np.sqrt(AMT['EC_flux_u']**2),color='k',linestyle='None',marker= 'o')
    han2 = ax.errorbar(AMT['latitude'],AMT['vt'][:,vals[0]],np.sqrt(AMT['vts'][:,vals[0]]**2 + AMT['vtse'][:,vals[0]]**2),color='r',linestyle='None',marker= 'o')
    han3 = ax.errorbar(AMT['latitude'],AMT['nvt'][:,vals[1]],np.sqrt(AMT['nvts'][:,vals[1]]**2 + AMT['nvtse'][:,vals[1]]**2),color='b',linestyle='None',marker= 'o')
    opacity_bar(han)
    opacity_bar(han2)
    opacity_bar(han3)
    ax.invert_xaxis()
    ax.set_ylabel('CO$_2$ Flux (mmol m$^{-2}$ d$^{-1}$)')
    ax.set_xlabel('Latitude ($^o$N)')
    ax.grid()
    ax.text(0.03,0.95,r"{\fontsize{20}{22}\textbf{(" + let + ")}}",transform=ax.transAxes,va='top')
    #ax.text(58,18,,fontweight='bold',fontsize=18)
    ax.set_xlim([60,-60])
    ax.plot([60,-60],[0,0],'k--')
    ax.set_ylim(ylim)
    ax.legend([han,han2,han3],['Direct EC flux','Indirect flux with Vertical Temperature Gradients','Indirect flux with no Vertical Temperature Gradients'],loc=3,fontsize=13)

def plot_latitude_diff(ax,AMT,let,vals = [7,1]):
    han2 = ax.errorbar(AMT['latitude'],AMT['EC_flux']-AMT['vt'][:,vals[0]],np.sqrt(AMT['vts'][:,vals[0]]**2 + AMT['EC_flux_u']**2 ),color='r',linestyle='None',marker= 'o')
    han3 = ax.errorbar(AMT['latitude'],AMT['EC_flux']-AMT['nvt'][:,vals[1]],np.sqrt(AMT['nvts'][:,vals[1]]**2  + AMT['EC_flux_u']**2 ),color='b',linestyle='None',marker= 'o')
    ax.invert_xaxis()
    ax.set_ylabel('Direct EC - Indirect bulk CO$_2$ flux (mmol m$^{-2}$ d$^{-1}$)')
    ax.set_xlabel('Latitude ($^o$N)')
    ax.grid()
    ax.text(0.03,0.95,r"{\fontsize{20}{22}\textbf{(" + let + ")}}",transform=ax.transAxes,va='top')
    #ax.text(58,18,,fontweight='bold',fontsize=18)
    ax.set_xlim([60,-60])
    ax.plot([60,-60],[0,0],'k--')
    ax.set_ylim([-15,10])
    ax.legend([han2,han3],['Indirect flux with Vertical Temperature Gradients','Indirect flux with no Vertical Temperature Gradients'],loc=3,fontsize=13)

def plot_scatter(ax,AMT,let,file,lims=[-30,20],vals = [7,1]):
    lims = np.array(lims)
    if len(AMT['vt'].shape) >1:
        f = (np.isnan(AMT['EC_flux']) == 0) & (np.isnan(AMT['vt'][:,vals[0]]) ==0) &  (np.isnan(AMT['nvt'][:,vals[1]]) ==0)
    else:
        f = (np.isnan(AMT['EC_flux']) == 0) & (np.isnan(AMT['vt'][vals[0]]) ==0) &  (np.isnan(AMT['nvt'][vals[1]]) ==0)
    print(len(f))
    g = np.argwhere((AMT['EC_flux'][f] < 5) & (AMT['EC_flux'][f] > -5))
    #print(len(g)/len(AMT['EC_flux'][f]))

    if len(f) > 1:
        han =ax.errorbar(AMT['EC_flux'][f],AMT['vt'][f,vals[0]],xerr=np.sqrt(AMT['EC_flux_u'][f]**2),yerr= np.sqrt(AMT['vts'][f,vals[0]]**2),color='r',linestyle='none',linewidth=0.6,marker = 'o',alpha = 0.5,markersize = 4)
        han2 = ax.errorbar(AMT['EC_flux'][f],AMT['nvt'][f,vals[1]],xerr=np.sqrt(AMT['EC_flux_u'][f]**2),yerr= np.sqrt(AMT['nvts'][f,vals[1]]**2),color='b',linestyle='none',linewidth=0.6,marker = 'o',alpha = 0.5,markersize = 4)
        vt = ws.unweighted_stats(AMT['EC_flux'][f],AMT['vt'][f,vals[0]],'vt')
        nvt = ws.unweighted_stats(AMT['EC_flux'][f],AMT['nvt'][f,vals[1]],'nvt')
    else:
        han =ax.errorbar(AMT['EC_flux'][f],AMT['vt'][vals[0]],xerr=np.sqrt(AMT['EC_flux_u'][f]**2),yerr= np.sqrt(AMT['vts'][vals[0]]**2),color='r',linestyle='none',linewidth=0.6,marker = 'o',alpha = 0.5,markersize = 4)
        han2 = ax.errorbar(AMT['EC_flux'][f],AMT['nvt'][vals[1]],xerr=np.sqrt(AMT['EC_flux_u'][f]**2),yerr= np.sqrt(AMT['nvts'][vals[1]]**2),color='b',linestyle='none',linewidth=0.6,marker = 'o',alpha = 0.5,markersize = 4)
        vt = ws.unweighted_stats(AMT['EC_flux'][f],np.array(AMT['vt'][vals[0]]),'vt')
        nvt = ws.unweighted_stats(AMT['EC_flux'][f],np.array(AMT['nvt'][vals[1]]),'nvt')
    opacity_bar(han)
    opacity_bar(han2)
    ax.plot(lims,lims*vt['slope'] + vt['intercept'],'r--')
    ax.plot(lims,lims*nvt['slope'] + nvt['intercept'],'b--')
    ax.plot(lims,lims,'k-.')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.grid()
    ax.set_ylabel('Indirect CO$_2$ Flux (mmol m$^{-2}$ d$^{-1}$)')
    ax.set_xlabel('Direct EC CO$_2$ Flux (mmol m$^{-2}$ d$^{-1}$)')
    ax.legend([han,han2],['Vertical Temperature Gradients','No Vertical Temperature Gradients'],loc=4,fontsize=13)
    ax.text(0.03,0.95,r"{\fontsize{20}{22}\textbf{(" + let + ")}}" + '\n',transform=ax.transAxes,va = 'top')
    ax.text(0.11,0.95, r"{\fontsize{14}{15}\textbf{ Temperature Gradients (n = " + str(vt['n']) + ")}}" + '\n' + r"{\fontsize{14}{15}\textmd{RMSD =  "+str(round(vt['rmsd'],2)) + r" mmol m\textsuperscript{-2} d\textsuperscript{-1}} }" + '\n' \
    r'{\fontsize{14}{15}\textmd{Mean Bias = ' + str(round(vt['rel_bias'],2))+ r' mmol m\textsuperscript{-2} d\textsuperscript{-1} }}' + '\n' + r"{\fontsize{14}{15}\textbf{ No Temperature Gradients (n = " + str(nvt['n']) + ")}}" + '\n' + r'{\fontsize{14}{15}\textmd{RMSD =  '+str(round(nvt['rmsd'],2)) + r' mmol m\textsuperscript{-2} d\textsuperscript{-1} }}' + '\n' \
    r'{\fontsize{14}{15}\textmd{Mean Bias = ' + str(round(nvt['rel_bias'],2))+ r' mmol m\textsuperscript{-2} d\textsuperscript{-1} }}' + '\n', transform=ax.transAxes,va = 'top',linespacing=0.9)

    field_names = vt.keys()
    names = ['NIGHTINGALE','WANNINKHOF','HO','YANG','S3A_Fred','CCI_avhrr_on_metop_a','CCI_s3a+b','CCI_merged']
    var = ['_COARE','_DON','_FIXED','_SKIN_COARE','_SKIN_DON','_SKIN_FIXED']
    out = []

    if len(AMT['vt'].shape) > 1:
        tp = 0
        tp2 = 0
        # for i in range(AMT['vt'].shape[1]):
        #     #print(tp)
        #     out.append(ws.weighted_stats(AMT['EC_flux'][f],AMT['vt'][f,i],1/np.sqrt((AMT['vts'][f,i]/AMT['vt'][f,i])**2 +(AMT['vtse'][f,i]/AMT['vt'][f,i])**2 + (AMT['EC_flux_u'][f]/AMT['EC_flux'][f])**2 + (AMT['EC_SE'][f]/AMT['EC_flux'][f])**2),names[tp2]+var[tp]+'_vt_weighted_relative'))
        #     out.append(ws.weighted_stats(AMT['EC_flux'][f],AMT['vt'][f,i],1/np.sqrt((AMT['vts'][f,i])**2 + (AMT['vtse'][f,i])**2 + (AMT['EC_flux_u'][f])**2 + (AMT['EC_SE'][f])),names[tp2]+var[tp]+'_vt_weighted_absolute'))
        #     tp = tp+1
        #     if tp == 4:
        #         tp = 0
        #         tp2 = tp2+1
        #
        # for i in range(AMT['nvt'].shape[1]):
        #     out.append(ws.weighted_stats(AMT['EC_flux'][f],AMT['nvt'][f,i],1/np.sqrt((AMT['nvts'][f,i]/AMT['nvt'][f,i])**2 + (AMT['nvtse'][f,i]/AMT['nvt'][f,i])**2 + (AMT['EC_flux_u'][f]/AMT['EC_flux'][f])**2 +(AMT['EC_SE'][f]/AMT['EC_flux'][f])**2),names[i]+'_nvt_weighted_relative'))
        #     out.append(ws.weighted_stats(AMT['EC_flux'][f],AMT['nvt'][f,i],1/np.sqrt((AMT['nvts'][f,i])**2 + (AMT['nvtse'][f,i])**2 + (AMT['EC_flux_u'][f])**2 + (AMT['EC_SE'][f])**2),names[i]+'_nvt_weighted_absolute'))
        tp = 0
        tp2 = 0
        for i in range(AMT['vt'].shape[1]):
            #print(tp)
            out.append(ws.unweighted_stats(AMT['EC_flux'][f],AMT['vt'][f,i],names[tp2]+var[tp]+'_vt'))
            tp = tp+1
            if tp == len(var):
                tp = 0
                tp2 = tp2+1
        for i in range(AMT['nvt'].shape[1]):
            out.append(ws.unweighted_stats(AMT['EC_flux'][f],AMT['nvt'][f,i],names[i]+'_nvt'))
        out.reverse()
    else:
        out = []

    with open(file,'w',newline='') as f:
        writer = csv.DictWriter(f, fieldnames = field_names)
        writer.writeheader()
        writer.writerows(out)

def opacity_bar(han):
    [bar.set_alpha(0.5) for bar in han[2]]
    [cap.set_alpha(0.5) for cap in han[1]]

def merge_EC(AMT):
    AMT_o = {}
    for key in AMT[0]:
        AMT_o[key] = np.concatenate((AMT[0][key],AMT[1][key]))
    a = np.argwhere(np.isnan(AMT_o['EC_flux'])  | np.isnan(AMT_o['vt'][:,2]) | np.isnan(AMT_o['nvt'][:,1]))
    for key in AMT_o:
        if len(AMT_o[key].shape) > 1:
            AMT_o[key][a,:] = np.nan
        else:
            AMT_o[key][a] = np.nan
    return AMT_o

def plot_SST(ax,AMT,let = 'c'):
    lims = [-60,60]
    lab = ['AMT28','AMT29']
    col = ['r','b']
    for i in range(len(AMT)):
        ax.errorbar(AMT[i]['latitude'],(AMT[i]['Tsubskin_don']) - AMT[i]['Tsst'],yerr = np.sqrt((AMT[i]['Tskin_unc']**2) + (0.1**2)),color=col[i],linestyle='none',marker='.',label=lab[i],linewidth=0.3,markersize=2)
        print('Tsubskin - Tdepth = ',str(np.nanmean((AMT[i]['Tsubskin_don']) - AMT[i]['Tsst'])))
    ax.plot(lims,[0,0],'k--')
    #ax.plot(lims,[-0.17,-0.17],'k-.',label='-0.17 K')
    ax.set_xlim(lims)
    ax.grid()
    ax.set_ylim([-1.5,1.5])
    ax.invert_xaxis()
    ax.set_ylabel('T$_{subskin}$ - T$_{depth}$ (K)',fontsize=20)
    ax.set_xlabel('Latitude ($^o$N)')
    ax.legend(loc=3,fontsize=13)
    ax.text(0.03,0.95,r"{\fontsize{20}{22}\textbf{(" + let + ")}}",transform=ax.transAxes,va='top')

def lat_divide(AMT2,maxlat,minlat,ax,let,file):
    #print(AMT2['latitude'])
    f = np.where((AMT2['latitude'] < maxlat) & (AMT2['latitude'] > minlat))
    #print(f)
    AMT3 = split_data(AMT2,f)
    #print(AMT3['vt'].shape)
    plot_scatter(ax,AMT3,let,file)
    #print(AMT3)
    ax.set_title('Latitude = ' + lat_title(maxlat) + ' - ' + lat_title(minlat),fontsize=15)

def split_data(AMT,f):
    AMT2 = {}
    for key in AMT:
        if len(AMT[key].shape) > 1:
            AMT2[key] = np.squeeze(AMT[key][f,:])
        else:
            AMT2[key] = AMT[key][f]
    return AMT2

def lat_title(lat):
    l = np.sign(lat)
    if l == 1:
        return str(lat*l)+ '$^o$N'
    else:
        return str(lat*l)+ '$^o$S'

def produce_time_window(time_s,time_e,window):
    out = [time_s]
    c = 0
    while c == 0:
        if out[-1] < time_e:
            out.append(out[-1] + window)
        else:
            c = 1
    return out

def diurnal_sst_plot(AMT,ax,subskin,let=''):
    c = plt.cm.get_cmap('RdYlBu')
    col = ['b','r']
    LAB = ['AMT28','AMT29']
    for j in range(2):
        av = np.empty((23))
        hour = np.empty((len(AMT[j]['time'])))
        for i in range(len(AMT[j]['time'])):
            hour[i] = AMT[j]['time'][i].hour + AMT[j]['time'][i].minute/60 + AMT[j]['time'][i].second/3600 + (AMT[j]['longitude'][i]/15)
            if np.sign(hour[i]) == -1.:
                hour[i] = 24+hour[i]
        #print(hour)
        a = ax.scatter(hour,AMT[j][subskin] - AMT[j]['Tsst'],c = AMT[j]['U10n'],vmin=0,vmax=15,cmap =c)
        for i in range(23):
            f = np.where((hour <= i + 1) & (hour > i))
            av[i] = np.nanmedian(AMT[j][subskin][f] - AMT[j]['Tsst'][f])
        ax.plot(np.arange(.5,23.5,1),av,color=col[j],label=LAB[j])
    ax.plot([0,24],[0.14,0.14],'k--')
    ax.plot([0,24],[-0.14,-0.14],'k--')
    ax.set_ylabel('T$_{subskin}$ - T$_{depth}$ (K)')
    ax.set_xlabel('Local Hour')
    ax.set_ylim([-0.5,1.75])
    ax.set_xlim([0,24])
    ax.text(0.90,0.95,r"{\fontsize{20}{22}\textbf{(" + let + ")}}",transform=ax.transAxes,va='top')
    ax.grid()
    c2 = plt.colorbar(a)
    c2.set_label('Wind Speed (ms$^{-1}$)')
    ax.legend(loc = 2)

def load_data_merge(window,cutoff):
    AMT28_3,time28 = load_3hr('AMT28/DATA/AMT28_3hrAvg_v2_2022-09-21.txt')
    #time28 = time28 - datetime.timedelta(hours=1.5)
    time28 = produce_time_window(time28[0],time28[-1],datetime.timedelta(hours = window))
    AMT28,AMT28_r,_ = data_load_feout('AMT28/DATA/FLUXENGINE_OUT.nc',time28,datetime.timedelta(hours=window),cutoff)
    # print(AMT28_3['latitude'])
    # print(AMT28['latitude'])

    AMT29_3,time29 = load_3hr('AMT29/DATA/AMT29_3hrAvg_v2_2022-09-21.txt')
    #time29 = time29 - datetime.timedelta(hours=1.5)
    time29 = produce_time_window(time29[0],time29[-1],datetime.timedelta(hours = window))
    AMT29,AMT29_r,_ = data_load_feout('AMT29/DATA/FLUXENGINE_OUT.nc',time29,datetime.timedelta(hours=window),cutoff)
    AMT = [AMT28,AMT29]
    AMT_r = [AMT28_r,AMT29_r]
    AMT_m = merge_EC(AMT)
    return AMT,AMT_r,AMT_m

cutoff = 0.33 # Smaller value is more restrictive (i.e 0.33 = 67% of data avaiable)
window = 3
AMT,AMT_r,AMT_m = load_data_merge(window,cutoff)

na = np.where((np.isnan(AMT_m['EC_flux']) == 0) & (np.isnan(AMT_m['vt'][:,7]) ==0) )
nvt = np.transpose(AMT_m['EC_flux'][na] - AMT_m['nvt'][na,1])
vt = np.transpose(AMT_m['EC_flux'][na] - AMT_m['vt'][na,7])
print(vt.shape)
print('Mann_Whitney')
u1,p = mannwhitneyu(vt,nvt)
print(p)
f = f_oneway(vt,nvt)
print(f)
vt = np.matlib.repmat(vt,25,1)
nvt = np.matlib.repmat(nvt,25,1)
print(vt.shape)
f = f_oneway(vt,nvt)
print(f)
#--------------------------------------------------------------------------------------------------------------------------------
# Diurnal Tsubskin - Tdepth plots
fig11 = plt.figure(figsize=(15,7))
gs = GridSpec(1,2, figure=fig11, wspace=0.33,hspace=0.2,bottom=0.1,top=0.95,left=0.1,right=0.95)
ax1 = fig11.add_subplot(gs[0,0])
diurnal_sst_plot(AMT_r,ax1,'Tsubskin_don',let='a')
ax2 = fig11.add_subplot(gs[0,1])
diurnal_sst_plot(AMT_r,ax2,'Tsubskin_co',let='b')
fig11.savefig('FIGS/Diurnal_Tsubskin_check.png',format='png',dpi=300)
#ax1.set_title('Donlon et al. (2002) Cool Skin')
#ax2.set_title('NOAA COARE 3.5 Cool Skin')

# AMT_r2 = []
# for i in range(2):
#     f = np.where((AMT_r[i]['latitude'] < 60) & (AMT_r[i]['latitude'] > 30))
#     AMT_r2.append(split_data(AMT_r[i],f))
#
# ax1 = fig11.add_subplot(gs[1,0])
# diurnal_sst_plot(AMT_r2,ax1,'Tsubskin_don')
# ax2 = fig11.add_subplot(gs[1,1])
# diurnal_sst_plot(AMT_r2,ax2,'Tsubskin_co')
# ax1.set_title('60N - 30N')
# ax2.set_title('60N - 30N')
#
# AMT_r2 = []
# for i in range(2):
#     f = np.where((AMT_r[i]['latitude'] < 30) & (AMT_r[i]['latitude'] > -30))
#     AMT_r2.append(split_data(AMT_r[i],f))
#
# ax1 = fig11.add_subplot(gs[2,0])
# diurnal_sst_plot(AMT_r2,ax1,'Tsubskin_don')
# ax2 = fig11.add_subplot(gs[2,1])
# diurnal_sst_plot(AMT_r2,ax2,'Tsubskin_co')
# ax1.set_title('30N - 30S')
# ax2.set_title('30N - 30S')
#
# AMT_r2 = []
# for i in range(2):
#     f = np.where((AMT_r[i]['latitude'] < -30) & (AMT_r[i]['latitude'] > -60))
#     AMT_r2.append(split_data(AMT_r[i],f))
#
# ax1 = fig11.add_subplot(gs[3,0])
# diurnal_sst_plot(AMT_r2,ax1,'Tsubskin_don')
# ax2 = fig11.add_subplot(gs[3,1])
# diurnal_sst_plot(AMT_r2,ax2,'Tsubskin_co')
# ax1.set_title('30S - 60S')
# ax2.set_title('30S - 60S')
#--------------------------------------------------------------------------------------------------------------------------------
# Main figure 1

fig = plt.figure(figsize=(21,15))
gs = GridSpec(2 ,3, figure=fig, wspace=0.2,hspace=0.2,bottom=0.05,top=0.97,left=0.07,right=0.95) # Setup a 2x2 subplot grid
ax1 = fig.add_subplot(gs[0,1]) # Set top row into a single plot
ax2 = fig.add_subplot(gs[1,1])
ax3 = fig.add_subplot(gs[1,0]) # Set bottom row into two plots
ax4 = fig.add_subplot(gs[0,0])
ax5 = fig.add_subplot(gs[0,2])
ax6 = fig.add_subplot(gs[1,2])
ax1.set_title('AMT28')
ax2.set_title('AMT29')
ax5.set_title('AMT28')
ax6.set_title('AMT29')
plot_latitude(ax1,AMT[0],'b',ylim=[-15,10])
plot_latitude(ax2,AMT[1],'e')
plot_latitude_diff(ax5,AMT[0],'c')
plot_latitude_diff(ax6,AMT[1],'f')
plot_scatter(ax3,AMT_m,'d','STATS/ALL.csv')
mapo.plotmap(ax4)
fig.savefig('FIGS/FIG_2_EC_FE_FLUX_COMPARISION.png',format='png',dpi=300)

#--------------------------------------------------------------------------------------------------------------------------------
# Latitudinal splits

# fig2 = plt.figure(figsize=(19.20,10.80))
# gs = GridSpec(2,3, figure=fig2, wspace=0.3,hspace=0.3,bottom=0.08,top=0.95,left=0.06,right=0.95) # Setup a 2x2 subplot grid
# ax1 = fig2.add_subplot(gs[0,0]) # Set top row into a single plot
# ax2 = fig2.add_subplot(gs[0,1])
# ax3 = fig2.add_subplot(gs[0,2]) # Set bottom row into two plots
# ax4 = fig2.add_subplot(gs[1,0])
# ax5 = fig2.add_subplot(gs[1,1])
# #ax6 = fig2.add_subplot(gs[1,2])
# AMT_2 = AMT_m
# lat_divide(AMT_2,60,40,ax1,'a','STATS/60_45.csv')
# lat_divide(AMT_2,40,12,ax2,'b','STATS/45_12.csv')
# lat_divide(AMT_2,12,-10,ax3,'c','STATS/12_-5.csv')
# lat_divide(AMT_2,-10,-35,ax4,'d','STATS/-5_-35.csv')
# lat_divide(AMT_2,-35,-60,ax5,'e','STATS/-35_-60.csv')
# fig2.savefig('FIGS/FIG_3_LAT_SPLIT_COMPARISION.png',format='png',dpi=300)
#--------------------------------------------------------------------------------------------------------------------------------
# Cool Skin effects plot

fig3 = plt.figure(figsize=(10,15))
gs = GridSpec(2,1, figure=fig3, wspace=0.15,hspace=0.2,bottom=0.08,top=0.95,left=0.1,right=0.95) # Setup a 2x2 subplot grid
ax1 = fig3.add_subplot(gs[0,0])
ax2 = fig3.add_subplot(gs[1,0])
ax1.scatter(AMT_r[0]['latitude'],AMT_r[0]['Tsubskin_co'] - AMT_r[0]['Tskin']+273.15,color='r')
ax1.scatter(AMT_r[0]['latitude'],AMT_r[0]['Tsubskin_don'] - AMT_r[0]['Tskin']+273.15,color='b')
ax1.plot([-50,50],[0.17,0.17],'k--')
ax2.scatter(AMT_r[1]['latitude'],AMT_r[1]['Tsubskin_co'] - AMT_r[1]['Tskin']+273.15,color='r')
ax2.scatter(AMT_r[1]['latitude'],AMT_r[1]['Tsubskin_don'] - AMT_r[1]['Tskin']+273.15,color='b')
ax2.plot([-50,50],[0.17,0.17],'k--')
ax1.invert_xaxis()
ax1.set_ylabel('AMT28 Cool Skin (K)')
ax1.set_xlabel('Latitude ($^o$N)')
ax1.grid()
ax1.set_xlim([60,-60])
ax1.set_ylim([0,0.6])
ax1.legend(['NOAA COARE 3.5','Donlon et al. (2002)','Fixed Cool Skin'],loc=1,fontsize=13)
ax1.text(0.03,0.95,r"{\fontsize{20}{22}\textbf{(a)}}",transform=ax1.transAxes,va='top')
ax2.text(0.03,0.95,r"{\fontsize{20}{22}\textbf{(b)}}",transform=ax2.transAxes,va='top')

ax2.invert_xaxis()
ax2.set_ylabel('AMT29 Cool Skin (K)')
ax2.set_xlabel('Latitude ($^o$N)')
ax2.grid()
ax2.set_ylim([0,0.6])
ax2.set_xlim([60,-60])
ax2.legend(['NOAA COARE 3.5','Donlon et al. (2002)','Fixed Cool Skin'],loc=1,fontsize=13)
fig3.savefig('FIGS/SUP_FIG_1_COARE_COOL_SKIN.png',format='png',dpi=300)
#--------------------------------------------------------------------------------------------------------------------------------
# Per cruise plots
##
fig4 = plt.figure(figsize=(15,5))
gs = GridSpec(1,2, figure=fig4, wspace=0.33,hspace=0.2,bottom=0.14,top=0.95,left=0.1,right=0.95) # Setup a 2x2 subplot grid
ax1 = fig4.add_subplot(gs[0,0])
ax2 = fig4.add_subplot(gs[0,1])
plot_scatter(ax1,AMT[0],'a','STATS/AMT28.csv')
plot_scatter(ax2,AMT[1],'b','STATS/AMT29.csv')
fig4.savefig('FIGS/SUP_FIG_4_PER_CRUISE_PLOT.png',format='png',dpi=300)
#--------------------------------------------------------------------------------------------------------------------------------
# # Low dpCO2 regions
# ##
# fig5 = plt.figure(figsize=(7,5))
# gs = GridSpec(1,1, figure=fig5, wspace=0.1,hspace=0.1,bottom=0.14,top=0.95,left=0.18,right=0.95)
# ax1 = fig5.add_subplot(gs[0,0])
# f = np.squeeze(np.argwhere(np.abs(AMT_m['fco2_atm'] - AMT_m['fco2_sw']) < 5))
# AMT_dpco2 = split_data(AMT_m,f)
# plot_scatter(ax1,AMT_dpco2,'a','STATS/low_dpco2.csv',lims=[-8,8])
# fig5.savefig('FIGS/SUP_FIG_5_DPCO2_0.png',format='png',dpi=300)
#--------------------------------------------------------------------------------------------------------------------------------
# # Weighted stats, weightings plot.
#
# fig6 = plt.figure(figsize=(15,15))
# gs = GridSpec(2,2, figure=fig6, wspace=0.33,hspace=0.2,bottom=0.1,top=0.95,left=0.1,right=0.9)
# ax1 = fig6.add_subplot(gs[0,0])
# ax2 = fig6.add_subplot(gs[0,1])
# ax3 = fig6.add_subplot(gs[1,0])
# ax4 = fig6.add_subplot(gs[1,1])
# ax1.scatter(AMT_m['EC_flux'],np.sqrt(AMT_m['EC_flux_u']**2 + AMT_m['EC_SE']**2))
# ax1.scatter(AMT_m['vt'][:,4],np.sqrt(AMT_m['vts'][:,4]**2 + AMT_m['vtse'][:,4]**2))
# ax1.set_ylabel('Flux Uncertainty (mmol m$^{-2}$ d$^{-1}$)')
# ax1.set_xlabel('CO$_2$ Flux (mmol m$^{-2}$ d$^{-1}$)')
# ax1.grid()
# ax1.legend(['Flux$_{EC}$', 'Flux$_{Bulk}$'],loc=2)
# ax2.scatter(AMT_m['EC_flux'],np.sqrt((AMT_m['EC_flux_u']/AMT_m['EC_flux'])**2 + (AMT_m['EC_SE']/AMT_m['EC_flux'])**2)*100)
# ax2.scatter(AMT_m['vt'][:,4],np.sqrt((AMT_m['vts'][:,4]/AMT_m['vt'][:,4])**2 + (AMT_m['vtse'][:,4]/AMT_m['vt'][:,4])**2)*100)
# ax2.set_ylabel('Relative Flux Uncertainty (\%)')
# ax2.set_xlabel('CO$_2$ Flux (mmol m$^{-2}$ d$^{-1}$)')
# ax2.set_ylim([0,300])
# ax2.grid()
# ax2.legend(['Flux$_{EC}$', 'Flux$_{Bulk}$'],loc=2)
#
# weight = 1/np.sqrt(AMT_m['EC_flux_u']**2 + AMT_m['vts'][:,4]**2 + AMT_m['EC_SE']**2 + AMT_m['vtse'][:,4]**2)
# weight = weight/np.nansum(weight)
# ax3.scatter(AMT_m['EC_flux'],weight)
#
# weight = 1/np.sqrt((AMT_m['EC_flux_u']/AMT_m['EC_flux'])**2 + (AMT_m['vts'][:,4]/AMT_m['vt'][:,4])**2 + (AMT_m['vtse'][:,4]/AMT_m['vt'][:,4])**2 +(AMT_m['EC_SE']/AMT_m['EC_flux'])**2)
# weight = weight/np.nansum(weight)
# ax3.scatter(AMT_m['EC_flux'],weight)
# ax3.legend(['Absolute Uncertainty','Relative Uncertainty'],loc=2)
# ax3.grid()
# ax3.set_ylabel('Weight')
# ax3.set_xlabel('EC CO$_2$ Flux (mmol m$^{-2}$ d$^{-1}$)')
#
# ncp = pd.read_excel('NCP_GPP_R_AMT_22_06_2020_DANIEL_FORD.xlsx',sheet_name=1)
# h2 = ax4.scatter(ncp['NCP'][7:],np.abs(ncp['error.1'][7:]),label='Uncertainty')
# ax5 = ax4.twinx()
#
# weights = 1/np.abs(ncp['error.1'][7:])
# weights = weights/np.nansum(weights)
# h = ax5.scatter(ncp['NCP'][7:],weights,label='Weights',color='tab:orange')
# lns = [h,h2]
# labs = [l.get_label() for l in lns]
# ax4.legend(lns, labs, loc=2)
#
# ax5.set_ylabel('Weights')
# ax4.set_xlim([-250,250])
# ax4.set_xlabel('NCP (mmol O$_2$ m$^{-2}$ d$^{-1}$)')
# ax4.set_ylabel('NCP Uncertainty (mmol O$_2$ m$^{-2}$ d$^{-1}$)')
# ax4.grid()
#
# ax1.text(0.90,0.90,r"{\fontsize{20}{22}\textbf{(a)}}",fontweight='bold',fontsize=22,transform=ax1.transAxes,horizontalalignment='center')
# ax2.text(0.90,0.90,r"{\fontsize{20}{22}\textbf{(b)}}",fontweight='bold',fontsize=22,transform=ax2.transAxes,horizontalalignment='center')
# ax3.text(0.90,0.90,r"{\fontsize{20}{22}\textbf{(c)}}",fontweight='bold',fontsize=22,transform=ax3.transAxes,horizontalalignment='center')
# ax4.text(0.90,0.90,r"{\fontsize{20}{22}\textbf{(d)}}",fontweight='bold',fontsize=22,transform=ax4.transAxes,horizontalalignment='center')
# fig6.savefig('FIGS/SUP_FIG_7_EC_FLUX_UNC.png',format='png',dpi=300)
#--------------------------------------------------------------------------------------------------------------------------------
# # EC uncertainty discussion plots
# ##
# fig7 = plt.figure(figsize=(15,7))
# gs = GridSpec(1,2, figure=fig7, wspace=0.33,hspace=0.2,bottom=0.1,top=0.95,left=0.1,right=0.9)
# ax1 = fig7.add_subplot(gs[0,0])
# ax2 = fig7.add_subplot(gs[0,1])
#
# ax1.scatter(AMT_m['EC_flux'],AMT_m['EC_flux_u'],label='Uncertainty')
# ax1.scatter(AMT_m['EC_flux'],AMT_m['EC_SE'],label='Standard Deviation')
# ax1.set_xlabel('EC CO$_2$ Flux (mmol m$^{-2}$ d$^{-1}$)')
# ax1.set_ylabel('EC Uncertainty (mmol m$^{-2}$ d$^{-1}$)')
# ax1.legend()
# ax1.grid()
# ax2.scatter(AMT_m['EC_flux_u'],AMT_m['EC_SE'])
# ax2.plot([0,8],[0,8],'k--')
# ax2.set_xlabel('EC Uncertainty (mmol m$^{-2}$ d$^{-1}$)')
# ax2.set_ylabel('EC Standard Deviation (mmol m$^{-2}$ d$^{-1}$)')
# ax2.grid()
# ax2.set_ylim([0,8])
# ax2.set_xlim([0,8])
# fig7.savefig('FIGS/SUP_FIG_8_EC_FLUX_VARIABILITY.png',format='png',dpi=300)
#--------------------------------------------------------------------------------------------------------------------------------
# Main wind speed split plot
##
def wind_speed_split(ax,window,split_val,let=['','',''],val=[7,1]):
    AMT,AMT_r,AMT_m = load_data_merge(window,cutoff)
    f = np.where((AMT_m['ws'] <= split_val[0]))
    AMT2 = split_data(AMT_m,f)
    plot_scatter(ax[0],AMT2,let[0],'STATS/ws' + str(split_val[0]) + '_' + str(window) +'h.csv',vals=val)
    f = np.where((AMT_m['ws'] > split_val[0]) & (AMT_m['ws'] <= split_val[1]))
    AMT2 = split_data(AMT_m,f)
    plot_scatter(ax[1],AMT2,let[1],'STATS/ws' + str(split_val[0]) + '_' + str(split_val[1]) + '_' + str(window) +'h.csv',vals=val)
    f = np.where((AMT_m['ws'] > split_val[1]))
    AMT2 = split_data(AMT_m,f)
    plot_scatter(ax[2],AMT2,let[2],'STATS/ws'+str(split_val[1]) + '_' + str(window) +'h.csv',vals=val)

    ax[0].set_title('U$_{10}$ $\leq$' +str(split_val[0]) + 'ms$^{-1}$')
    ax[1].set_title(str(split_val[0]) + 'ms$^{-1}$ $<$ U$_{10}$ $\leq$ ' +str(split_val[1]) + 'ms$^{-1}$')
    ax[2].set_title( str(split_val[1]) + 'ms$^{-1}$ $<$ U$_{10}$')
    ax[0].set_ylabel('Bulk CO$_2$ Flux (mmol m$^{-2}$ d$^{-1}$)')

fig10 = plt.figure(figsize=(21,7))
gs = GridSpec(1,3, figure=fig10, wspace=0.33,hspace=0.2,bottom=0.1,top=0.95,left=0.1,right=0.95)
ax1 = fig10.add_subplot(gs[0,0])
ax2 = fig10.add_subplot(gs[0,1])
ax3 = fig10.add_subplot(gs[0,2])
# ax4 = fig10.add_subplot(gs[1,0])
# ax5 = fig10.add_subplot(gs[1,1])
# ax6 = fig10.add_subplot(gs[1,2])
# ax7 = fig10.add_subplot(gs[2,0])
# ax8 = fig10.add_subplot(gs[2,1])
# ax9 = fig10.add_subplot(gs[2,2])
wind_speed_split([ax1,ax2,ax3],window,[5,11],let = ['a','b','c'],val = [7,1])
# wind_speed_split([ax4,ax5,ax6],3,[5,11])
# wind_speed_split([ax7,ax8,ax9],3,[5,12])
fig10.savefig('FIGS/SUP_FIG_12_WIND_SPEED_REGIMES_IAN.png',format='png',dpi=300)
#--------------------------------------------------------------------------------------------------------------------------------
# # Wind speed splits Ian
# fig12 = plt.figure(figsize=(21,28))
# gs = GridSpec(4,3, figure=fig12, wspace=0.33,hspace=0.23,bottom=0.05,top=0.97,left=0.1,right=0.95)
# ax1 = fig12.add_subplot(gs[0,0])
# ax2 = fig12.add_subplot(gs[0,1])
# ax3 = fig12.add_subplot(gs[0,2])
# #ax4 = fig12.add_subplot(gs[1,:])
# ax4 = fig12.add_subplot(gs[1,0])
# ax5 = fig12.add_subplot(gs[1,1])
# ax6 = fig12.add_subplot(gs[1,2])
# ax7 = fig12.add_subplot(gs[2,0])
# ax8 = fig12.add_subplot(gs[2,1])
# ax9 = fig12.add_subplot(gs[2,2])
# ax10 = fig12.add_subplot(gs[3,0])
# ax11 = fig12.add_subplot(gs[3,1])
# ax12 = fig12.add_subplot(gs[3,2])
# # ax13 = fig12.add_subplot(gs[4,0])
# # ax14 = fig12.add_subplot(gs[4,1])
# # ax15 = fig12.add_subplot(gs[4,2])
# wind_speed_split([ax1,ax2,ax3],3,[5,11],let=['a','b','c'])
# wind_speed_split([ax4,ax5,ax6],6,[5,11],let=['d','e','f'])
# wind_speed_split([ax7,ax8,ax9],12,[5,11],let=['g','h','i'])
# wind_speed_split([ax10,ax11,ax12],24,[5,11],let=['j','k','l'])
# fig12.savefig('FIGS/SUP_FIG_12_WIND_SPEED_REGIMES_JAMIE.png',format='png',dpi=300)

#--------------------------------------------------------------------------------------------------------------------------------
# #Plotting uncertainty and standard deviation at different averaging windows.
# # ##
# wins = np.arange(1,24,1.0)
# wins = np.concatenate((np.array(0.6),wins),axis=None)
# print(wins)
# lap = [-60,40,15,-15,-40,-60]
# lapt = [60,60,40,15,-15,-40]
# out = np.empty((len(wins),6,len(lap)))
# for i in range(0,len(wins)):
#     AMT,AMT_r,AMT_m = load_data_merge(wins[i],cutoff)
#     for j in range(0,len(lap)):
#         f = np.argwhere((AMT_m['latitude'] < lapt[j]) & (AMT_m['latitude'] > lap[j]))
#         out[i,:,j] = np.array((np.nanmean(AMT_m['EC_flux_u'][f]),np.nanstd(AMT_m['EC_flux_u'][f]),np.nanmean(AMT_m['EC_SE'][f]),np.nanstd(AMT_m['EC_SE'][f]),np.nanmean(np.sqrt(AMT_m['EC_SE'][f]**2 + AMT_m['EC_flux_u'][f]**2)),np.nanstd(np.sqrt(AMT_m['EC_SE'][f]**2 + AMT_m['EC_flux_u'][f]**2))))
# #print(out)
#
# fig8 = plt.figure(figsize=(21,15))
# gs = GridSpec(2,3, figure=fig8, wspace=0.3,hspace=0.2,bottom=0.05,top=0.95,left=0.05,right=0.95)
# ax1 = fig8.add_subplot(gs[0,0])
# ax2 = fig8.add_subplot(gs[0,1])
# ax3 = fig8.add_subplot(gs[0,2])
# ax4 = fig8.add_subplot(gs[1,0])
# ax5 = fig8.add_subplot(gs[1,1])
# ax6 = fig8.add_subplot(gs[1,2])
# ax = [ax1,ax2,ax3,ax4,ax5,ax6]
# #
# tit = ['ALL','60N - 40N', '40N - 15N', '15N - 15S', '15S - 40S', '40S - 60S']
# for i in range(0,len(lap)):
#     h = ax[i].errorbar(wins,out[:,0,i],yerr=out[:,1,i],label='Measurement Unc')
#     h2 = ax[i].errorbar(wins,out[:,2,i],yerr=out[:,3,i],label='Standard Deviation')
#     h3 = ax[i].errorbar(wins,out[:,4,i],yerr=out[:,5,i],label='Total Unc')
#     ax[i].grid()
#     ax[i].legend()
#     ax[i].set_ylim([0,8])
#     ax[i].set_title(tit[i])
#     ax[i].set_xlabel('Averaging Windows (hours)')
#     ax[i].set_ylabel('Uncertainty (mmol m$^{-2}$ d$^{-1}$)')
#     ax[i].set_xlim([0,24])
#     ap = ax[i].twinx()
#     # j = np.empty((len(wins)))
#     # for k in range(0,len(wins)):
#     #     j[k] = unc_snr(out[k,0,i],out[k,1,i],out[k,2,i],out[k,3,i],100)
#     h4 = ap.errorbar(wins,out[:,2,i] / out[:,0,i],c='r',yerr = (out[:,2,i] / out[:,0,i]) * np.sqrt((out[:,1,i] / out[:,0,i])**2 + (out[:,3,i]/out[:,2,i])**2),linestyle='',marker='o',label='S:N')
#     ap.set_ylim([0,20])
#     ap.plot(wins[[0,-1]],[3,3],'k--')
#     ap.set_ylabel('Signal:Noise Ratio')
#     lns = [h,h2,h3,h4]
#     labs = [l.get_label() for l in lns]
#     ax[i].legend(lns, labs, loc=0)
# fig8.savefig('FIGS/SUP_FIG_9_EC_FLUX_VARIABILITY_TIME.png',format='png',dpi=300)
#--------------------------------------------------------------------------------------------------------------------------------
# # # Plot different temporal averages as one plot
# wins = [1,3,6,9,12,15,18,24]
# #wins = [1,3,6,9]
# fig9 = plt.figure(figsize=(28,15))
# gs = GridSpec(2,4, figure=fig9, wspace=0.3,hspace=0.2,bottom=0.05,top=0.95,left=0.05,right=0.95)
# ax1 = fig9.add_subplot(gs[0,0])
# ax2 = fig9.add_subplot(gs[0,1])
# ax3 = fig9.add_subplot(gs[0,2])
# ax4 = fig9.add_subplot(gs[0,3])
# ax5 = fig9.add_subplot(gs[1,0])
# ax6 = fig9.add_subplot(gs[1,1])
# ax7 = fig9.add_subplot(gs[1,2])
# ax8 = fig9.add_subplot(gs[1,3])
# ax = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
# let = ['a','b','c','d','e','f','g','h']
# for i in range(0,len(wins)):
#     AMT,AMT_r,AMT_m = load_data_merge(wins[i],cutoff)
#     plot_scatter(ax[i],AMT_m,let[i],'STATS/test.csv')
#     ax[i].set_title(str(wins[i])+' Hours')
# fig9.savefig('FIGS/SUP_FIG_10_RMSD_BIAS_CHECK.png',format='png',dpi=300)

fig13 = plt.figure(figsize=(15,15))
gs = GridSpec(2,2, figure=fig13, wspace=0.3,hspace=0.2,bottom=0.1,top=0.95,left=0.15,right=0.95)
ax1 = fig13.add_subplot(gs[0,0])
ax2 = fig13.add_subplot(gs[0,1])
ax3 = fig13.add_subplot(gs[1,0])
ax4 = fig13.add_subplot(gs[1,1])
AMT,AMT_r,AMT_m = load_data_merge(window,cutoff)
plot_scatter(ax1,AMT_m,'a','STATS/ALL_SAT.csv',vals=[31,1])
ax1.set_title('S3A_Tskin')
plot_scatter(ax2,AMT_m,'b','STATS/ALL_SAT.csv',vals=[28,1])
ax2.set_title('CCI_avhrr_on_metop_a')
plot_scatter(ax3,AMT_m,'c','STATS/ALL_SAT.csv',vals=[34,1])
ax3.set_title('CCI_s3a+b')
plot_scatter(ax4,AMT_m,'d','STATS/ALL_SAT.csv',vals=[25,1])
ax4.set_title('CCI_merged')
fig13.savefig('FIGS/S3A_SAT_CHECK.png',format='png',dpi=300)
# plt.show()

#---------------------------------------------------------------------------------------
def wind_speed_split_mapping(ax,window,split_val,let=['','',''],val=[5,1]):
    AMT,AMT_r,AMT_m = load_data_merge(window,cutoff)
    f = np.where((AMT_m['ws'] <= split_val[0]))
    AMT2 = split_data(AMT_m,f)
    g = np.where((np.isnan(AMT2['EC_flux']) == 0) & (np.isnan(AMT2['vt'][:,5]) == 0))
    mapo.plotbase_scatter(ax[0],AMT2['longitude'][g],AMT2['latitude'][g],let=let[0])
    f = np.where((AMT_m['ws'] > split_val[0]) & (AMT_m['ws'] <= split_val[1]))
    AMT2 = split_data(AMT_m,f)
    g = np.where((np.isnan(AMT2['EC_flux']) == 0) & (np.isnan(AMT2['vt'][:,5]) == 0))
    mapo.plotbase_scatter(ax[1],AMT2['longitude'][g],AMT2['latitude'][g],let=let[1])
    f = np.where((AMT_m['ws'] > split_val[1]))
    AMT2 = split_data(AMT_m,f)
    g = np.where((np.isnan(AMT2['EC_flux']) == 0) & (np.isnan(AMT2['vt'][:,5]) == 0))
    mapo.plotbase_scatter(ax[2],AMT2['longitude'][g],AMT2['latitude'][g],let=let[2])

    ax[0].set_title('U$_{10}$ $\leq$' +str(split_val[0]) + 'ms$^{-1}$')
    ax[1].set_title(str(split_val[0]) + 'ms$^{-1}$ $<$ U$_{10}$ $\leq$ ' +str(split_val[1]) + 'ms$^{-1}$')
    ax[2].set_title( str(split_val[1]) + 'ms$^{-1}$ $<$ U$_{10}$')
    #ax[0].set_ylabel('Bulk CO$_2$ Flux (mmol m$^{-2}$ d$^{-1}$)')

fig14 = plt.figure(figsize=(21,10))
gs = GridSpec(1,3, figure=fig14, wspace=0.33,hspace=0.2,bottom=0.1,top=0.95,left=0.1,right=0.95)
ax1 = fig14.add_subplot(gs[0,0])
ax2 = fig14.add_subplot(gs[0,1])
ax3 = fig14.add_subplot(gs[0,2])
# ax4 = fig10.add_subplot(gs[1,0])
# ax5 = fig10.add_subplot(gs[1,1])
# ax6 = fig10.add_subplot(gs[1,2])
# ax7 = fig10.add_subplot(gs[2,0])
# ax8 = fig10.add_subplot(gs[2,1])
# ax9 = fig10.add_subplot(gs[2,2])
wind_speed_split_mapping([ax1,ax2,ax3],window,[5,11],let = ['a','b','c'])#,val = [17,3])
# wind_speed_split([ax4,ax5,ax6],3,[5,11])
# wind_speed_split([ax7,ax8,ax9],3,[5,12])
fig14.savefig('FIGS/SUP_FIG_wind_spped_maps.png',format='png',dpi=300)

#------------------------------------------------------------------------------
# for i in range(1,24,2):
#     AMT,AMT_r,AMT_m = load_data_merge(i,cutoff)
#     print('Mann_Whitney')
#     print('Hours = '+str(i))
#     na = np.where((np.isnan(AMT_m['EC_flux']) == 0) & (np.isnan(AMT_m['vt'][:,4]) ==0) )
#     u1,p = mannwhitneyu(np.transpose(AMT_m['EC_flux'][na] - AMT_m['vt'][na,4]), np.transpose(AMT_m['EC_flux'][na] - AMT_m['nvt'][na,1]))
#     print(u1)
#     print(p)
