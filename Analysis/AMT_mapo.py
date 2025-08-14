#!/usr/bin/env python3
# ---------------------------------
# Code created by Daniel Ford (d.ford@exeter.ac.uk) - 18/08/2022
#
#--------  Version 1.0
# - Initial Version
# ---------------------------------

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from netCDF4 import Dataset
import numpy as np
from pylr2 import regress2
import matplotlib.transforms
from mpl_toolkits.basemap import Basemap
import matplotlib.transforms
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.patches import PathPatch
import matplotlib.patheffects as pe

font = {'weight' : 'normal',
        'size'   : 19}
matplotlib.rc('font', **font)
matplotlib.rcParams['text.usetex'] = True

def load_data(file):

    c = Dataset(file,'r')
    lat = np.array(c.variables['lat'])
    lon = np.array(c.variables['lon'])
    ws = np.array(c.variables['U10n_mean'])
    fco2_sw = np.array(c.variables['fCO2_sw_mean'])
    fco2_atm = np.array(c.variables['fCO2_atm_mean'])
    Tskin = np.array(c.variables['Tskin_mean'])
    c.close()
    out = {
        'lat':lat,
        'lon':lon,
        'ws':ws,
        'fco2_sw':fco2_sw,
        'fco2_atm':fco2_atm,
        'Tskin': Tskin
    }
    return out

def plotmap(ax1):
    AMT28 = load_data('AMT28/DATA/AMT28_table_20min_version2_i.nc')
    AMT29 = load_data('AMT29/DATA/AMT28_table_20min_version2_i.nc')

    m = Basemap(projection='merc',llcrnrlat=-60.1,urcrnrlat=60.1,\
                llcrnrlon=-80,urcrnrlon=20,lat_ts=20,resolution='l',ax=ax1)
    # m.readshapefile('E:/Data/Longhurst/Longhurst_world_v4_2010','Long',drawbounds=False)
    #
    #
    #
    # provs = ['NADR','NASE','NATR','WTRA','SATL','SSTC','FKLD',]
    # dz = np.arange(0,len(provs)+2,1).astype(np.int64)
    # norm = plt.Normalize()
    # colors = plt.cm.plasma(norm(dz))
    # vals = np.array([-25,50,-20,35,-50,18,-30,3,-15,-25,-15,-42,-60,-55])
    # vals = np.reshape(vals,[-1,2])
    # x,y = m(vals[:,0],vals[:,1])
    # print(colors)
    # t = 0
    # for prov in provs:
    #     patches   = []
    #     vals = []
    #     for info, shape in zip(m.Long_info, m.Long):
    #         #print(info)
    #         if info['ProvCode'] == prov:
    #             patches.append( Polygon(np.array(shape), True) )
    #             #vals = np.append(vals,np.array(shape))
    #             #print(shape)
    #     print(prov)
    #     #print(vals)
    #     vals = np.reshape(vals,[-1,2])
    #     ax1.add_collection(PatchCollection(patches, facecolor= colors[t+1,:], edgecolor='k', linewidths=2., zorder=2))
    #     ax1.text(x[t],y[t],prov,ha='center',fontsize=16,fontweight='bold',color='w')
    #     t = t+1


    x, y = m(AMT28['lon'],AMT28['lat'])
    m.plot(x,y,label='AMT28',linestyle='--',color='r',linewidth=3,path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
    x, y = m(AMT29['lon'],AMT29['lat'])
    m.plot(x,y,color='b',label='AMT29',linestyle='-.',linewidth=3,path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
    x,y = m(-60,55)
    #ax1.annotate('(a)', xy=(0.24, .93), xycoords='axes fraction',fontsize=24,fontweight='bold')]
    ax1.text(-0.1,0.98,r"{\fontsize{20}{22}\textbf{(a)}}",transform=ax1.transAxes,va='top')
    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='#C2c3c3')
    m.drawcoastlines()
    parallels = np.arange(-50.,51,25.)
    meridion = [-50,-25,0]
    m.drawparallels(parallels,labels=[True,False,False,False],zorder=-2)
    m.drawmeridians(meridion,labels=[False,False,False,True],zorder=-2)
    ax1.legend(loc=4)
    ax1.set_xlabel('Longitude',labelpad=30)
    ax1.set_ylabel('Latitude',labelpad=30)
    _,parallels = m(parallels,parallels)
    ax1.set_yticks(parallels,direction='in')
    ax1.set_yticklabels(['','','','',''])

    meridion,_ = m(meridion,meridion)
    ax1.set_xticks(meridion,direction='in')
    ax1.set_xticklabels(['','',''])

def plotbase_scatter(ax1,lon,lat,let=[]):
    m = Basemap(projection='merc',llcrnrlat=-60.1,urcrnrlat=60.1,\
                llcrnrlon=-80,urcrnrlon=20,lat_ts=20,resolution='l',ax=ax1)
    # m.readshapefile('E:/Data/Longhurst/Longhurst_world_v4_2010','Long',drawbounds=False)
    #
    # provs = ['NADR','NASE','NATR','WTRA','SATL','SSTC','FKLD',]
    # dz = np.arange(0,len(provs)+2,1).astype(np.int64)
    # norm = plt.Normalize()
    # colors = plt.cm.plasma(norm(dz))
    # vals = np.array([-25,50,-20,35,-50,18,-30,3,-15,-25,-15,-42,-60,-55])
    # vals = np.reshape(vals,[-1,2])
    # x,y = m(vals[:,0],vals[:,1])
    # print(colors)
    # t = 0
    # for prov in provs:
    #     patches   = []
    #     vals = []
    #     for info, shape in zip(m.Long_info, m.Long):
    #         #print(info)
    #         if info['ProvCode'] == prov:
    #             patches.append(Polygon(np.array(shape), True))
    #             #vals = np.append(vals,np.array(shape))
    #             #print(shape)
    #     print(prov)
    #     #print(vals)
    #     vals = np.reshape(vals,[-1,2])
    #     ax1.add_collection(PatchCollection(patches, facecolor= colors[t+1,:], edgecolor='k', linewidths=2., zorder=2))
    #     ax1.text(x[t],y[t],prov,ha='center',fontsize=16,fontweight='bold',color='w')
    #     t = t+1
    x,y = m(lon,lat)
    ax1.scatter(x,y, zorder=3,s=48,color='r',marker='d',edgecolor='k')
    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='#C2c3c3')
    m.drawcoastlines()

    parallels = np.arange(-50.,51,25.)
    meridion = [-50,-25,0]
    m.drawparallels(parallels,labels=[True,False,False,False])
    m.drawmeridians(meridion,labels=[False,False,False,True])
    ax1.text(0.25,0.95,r"{\fontsize{20}{22}\textbf{("+let+")}}",transform=ax1.transAxes,va='top')
