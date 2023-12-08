#!/usr/bin/env python3
import numpy as np
from math import sqrt, cos, sin, radians
from netCDF4 import Dataset
import CCMP_HISTOGRAM as ccmp
import matplotlib.pyplot as plt

def oneDegreeArea(latDegrees):
    re, rp = 6378.137, 6356.7523
    dtor = radians(1.)
    #Taken from FLuxengine v4.0 - flux budgets code.
    # area of a 1x1 degree box at a given latitude in radians
    latRadians = latDegrees * dtor
    cosLat, sinLat = cos(latRadians), sin(latRadians)
    rc, rs = re * cosLat, rp * sinLat
    r2c, r2s = re * re * cosLat, rp * rp * sinLat
    earth_radius = sqrt((r2c * r2c + r2s * r2s) / (rc * rc + rs * rs))
    erd = earth_radius * dtor
    return erd * erd * cosLat

def load_land_area(file):
    c = Dataset(file)
    land = np.array(c.variables['land_proportion'])
    c.close()
    land = np.abs(land-1)
    return np.squeeze(land)

def latitude_band(lat1,lat2,lat_g,area):
    f = np.argwhere( (lat_g >= lat1) & (lat_g < lat2))
    out = np.nansum(area[f,:])
    return out

def box_band(lat1,lat2,lat_g,lon1,lon2,lon_g,area):
    f = np.squeeze(np.argwhere( (lat_g >= lat1) & (lat_g < lat2)))
    g = np.squeeze(np.argwhere( (lon_g >= lon1) & (lon_g < lon2)))
    area = area[f,:]
    area = area[:,g]
    print(area.shape)
    out = np.sum(area)
    return out

def grid_average(latg,long,data,nlat,nlon,res):
    d = data.shape
    out = np.empty((len(nlat),len(nlon),d[2]))
    out[:] = np.nan
    res=res/2
    if len(d) != 3:
        for i in range(0,len(nlat)):
            for j in range(0,len(nlon)):
                f = np.squeeze(np.argwhere( (latg >= nlat[i]-res) & (latg < nlat[i]+res) ))
                g = np.squeeze(np.argwhere( (long >= nlon[j]-res) & (long < nlon[j]+res) ))
                if (len(f) != 0) & (len(g) != 0):
                    out[i,j] = np.nanmean(data[f[:,None],g])
    else:
        for k in range(0,d[2]):
            print(k)
            for i in range(0,len(nlat)):
                for j in range(0,len(nlon)):
                    f = np.squeeze(np.argwhere( (latg >= nlat[i]-res) & (latg < nlat[i]+res) ))
                    g = np.squeeze(np.argwhere( (long >= nlon[j]-res) & (long < nlon[j]+res) ))
                    if (len(f) != 0) & (len(g) != 0):
                        out[i,j,k] = np.nanmean(data[f[:,None],g,k])
    return out

def grid_flip(data):
    d = data.shape
    out = np.empty((d))
    if len(d) != 3:
        out[:,0:int(np.floor(d[1]/2))] = data[:,int(np.floor(d[1]/2)):]
        out[:,int(np.floor(d[1]/2)):] = data[:,0:int(np.floor(d[1]/2))]
    else:
        for k in range(0,d[2]):
            out[:,0:int(np.floor(d[1]/2)),k] = data[:,int(np.floor(d[1]/2)):,k]
            out[:,int(np.floor(d[1]/2)):,k] = data[:,0:int(np.floor(d[1]/2)),k]
    return out

def wind_bias(ws,ws_split,bias):
    out = np.empty((ws.shape))
    out[:] = np.nan
    f = np.squeeze(np.argwhere( (ws <= ws_split[0])))
    #print(f.shape)
    out[f[:,0],f[:,1]] = bias[0]
    f = np.squeeze(np.argwhere( (ws > ws_split[0]) & (ws <= ws_split[1]) ))
    #print(f.shape)
    out[f[:,0],f[:,1]] = bias[1]
    f = np.squeeze(np.argwhere( (ws > ws_split[1])))
    #print(f)
    out[f[:,0],f[:,1]] = bias[2]
    return out

def load_sea_mask(file):
    c = Dataset(file)
    out = np.squeeze(np.array(c.variables['sea-mask']))
    c.close()
    return out

def load_gcb(file,var):
    c = Dataset(file)
    out = np.squeeze(np.array(c.variables[var]))
    c.close()
    #Convert flux from mol C m-2 s-1 to g C m-2 yr-1
    out = out * 12.01 * (3600 * 24 * 365)
    return out


loc = 'D:/Data/CCMP/v3.1/'
ws,latw,lonw = ccmp.ccmp_load(loc,2018,2018)
lonw = lonw-180
lat = np.arange(89.5,-90.5,-1).tolist()
lon = np.arange(-179.5,180.5,1).tolist()
ws = grid_average(latw,lonw,ws,np.array(lat),np.array(lon),1)
ws = grid_flip(ws)
ws[np.isnan(ws) == 1] = 7



land = load_land_area('C:/Users/df391/Anaconda3/envs/FluxEngine/Lib/site-packages/fluxengine/data/onedeg_land.nc')
mask = load_sea_mask('C:/Users/df391/Anaconda3/envs/FluxEngine/Lib/site-packages/fluxengine/data/World_Seas-IHO-mask.nc')
gcb_mod = load_gcb('GLOBAL_CARBON_BUDGET/Ocean_carbon_uptake_GOBMs_gridded_GCB2022_2012-2021_mean.nc','fgco2_A_avg')
gcb_mod = grid_flip(np.flipud(gcb_mod))


gcb_prod = load_gcb('GLOBAL_CARBON_BUDGET\Ocean_carbon_uptake_dataproducts_gridded_GCB2022_2012-2021_mean.nc','fgco2_ensemble_avg')
gcb_prod = grid_flip(np.flipud(gcb_prod))

area = np.empty((land.shape[0]))
for i in range(len(lat)):
    area[i] = oneDegreeArea(lat[i])
area = np.transpose(np.tile(area,(land.shape[1],1)))
lat_g = np.transpose(np.tile(lat,(land.shape[1],1)))
area_l = area*land * 1000000
area_lg = np.copy(area_l)
global_area = np.nansum(area_l)
area_l[mask != 30.0] = np.nan
area_l[lat_g>50.0] = np.nan
area_l[lat_g < -50] = np.nan
area_tot = np.nansum(area_l)
print(area_tot)

gcb_mod[mask!=30.0] = np.nan
gcb_prod[mask!=30.0] = np.nan
at_val = [np.nansum(gcb_mod*area_l)/1e15,np.nansum(gcb_prod*area_l)/1e15]
print('GCB_MOD = ' + str(at_val[0]))
print('GCB_PROD = ' + str(at_val[1]))
print('GCB_mean = ' + str(np.mean(at_val)))
# plt.figure()
# ws2 = ws.reshape(-1,1)
# plt.hist(ws2)
# #plt.show()
ws_b = np.empty((ws.shape))
ws_b2 = np.empty((ws.shape))
ws_b3 = np.empty((ws.shape))
for i in range(0,ws.shape[2]):
    ws_b[:,:,i] = wind_bias(ws[:,:,i],[5,11],[0.03,0.11,0.2])
    ws_b2[:,:,i] = wind_bias(ws[:,:,i],[5,11],[0.03,0.25,0.5])
    ws_b3[:,:,i] = wind_bias(ws[:,:,i],[5,11],[0.11,0.11,0.11])
# plt.figure()
# plt.pcolor(lon,lat,ws)
# plt.colorbar()
plt.figure()
for i in range(0,12):
    plt.subplot(3,4,i+1)
    plt.pcolor(lon,lat,ws_b[:,:,i])
plt.colorbar()
plt.show()
area_l2 = np.repeat(area_lg[:, :, np.newaxis], ws_b.shape[2], axis=2)
out = np.nansum(ws_b  * 30.5 / 1000 * 12.01 * area_l2) /1e15
print('Wind Speed Correct Global: ' + str(out))
out = np.nansum(ws_b2  * 30.5 / 1000 * 12.01 * area_l2) /1e15
print('Wind Speed Correct Global Skin: ' + str(out))
out = np.nansum(ws_b3  * 30.5 / 1000 * 12.01 * area_l2) /1e15
print('Wind Speed Correct Global Constant: ' + str(out))
area_band = []
area_band.append(latitude_band(-90,90,np.array(lat),area_l))
# area_band.append(latitude_band(-90,90,np.array(lat),area_l))
# area_band.append(latitude_band(-90,90,np.array(lat),area_l))
# area_band.append(latitude_band(-90,90,np.array(lat),area_l))
# area_band.append(latitude_band(-90,90,np.array(lat),area_l))
# area_band.append(latitude_band(-90,90,np.array(lat),area_l))
# #area_band.append(box_band(-44,90,np.array(lat),-70,25,np.array(lon),area_l))
#area_band.append(box_band(-44,90,np.array(lat),-70,25,np.array(lon),area_l))
# area_band.append(latitude_band(40,60,np.array(lat),area_l))
# area_band.append(latitude_band(15,40,np.array(lat),area_l))
# area_band.append(latitude_band(-15,15,np.array(lat),area_l))
# area_band.append(latitude_band(-40,-15,np.array(lat),area_l))
# area_band.append(latitude_band(-60,-40,np.array(lat),area_l))
print(np.array(area_band)/1000000)
print(global_area)
bias = np.array([0.19-0.09,0.19--0.07,-0.07-0.09,
    0.19-0.08, 0.19--0.07,-0.07-0.08,
    0.19-0.12,0.19--0.13,-0.13-0.12])
# bias_w = np.array([-0.29--0.4,-1.5--1.7,1.1-1.1,-0.52--0.52,-0.94--0.71,-0.32--0.054])
#print(bias)
bias_g = bias * 365.25 / 1000 * 12.01
# bias_w = bias_w *365.25 / 1000 * 12.01

co = bias_g * area_band / 1e15
# co2 = bias_w * area_band / 1e15
print(co)
print(co / np.mean(at_val)*100)

bias_g = bias * 365.25 / 1000 * 12.01
# bias_w = bias_w *365.25 / 1000 * 12.01

co = bias_g * global_area / 1e15
# co2 = bias_w * area_band / 1e15
print(co)
# print(co2)
# print(np.sum(co[1:]))
# print(np.sum(co2[1:]))
# print(np.sum(area_band[1:]))
