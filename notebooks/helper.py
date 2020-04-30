import xarray as xr
import datetime as dt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import geopandas as gpd
import pandas as pd
import numpy as np

def get_gdf(file, msurf):
    gdf =gpd.read_file(file)
    pd.set_option('mode.chained_assignment', None)
    gdf_ps = gdf[gdf['hru_id'].isin(msurf.var['nhm_id'].data)]
#     print(type(msurf.var['nhm_id'].data))
    dindex = np.zeros(np.shape(gdf_ps.hru_id.values), dtype=np.int8)
    for index, val in np.ndenumerate(msurf.var['nhm_id'].data):
        tind = np.int(np.where(gdf_ps['hru_id'].values == msurf.var['nhm_id'].data[index])[0])
    #     print(type(tind), tind)
        dindex[tind] = np.array([index])
    # print(dindex)
    gdf_ps.loc[:,'tindex'] = dindex
    gdf_ps.sort_values(by=['tindex'], inplace=True)
    return gdf_ps

def plot_climate(clim_file, gdf_ps, msurf):
    
    clim = xr.open_dataset(clim_file)
    ptime = msurf.var['nowtime'].data
    timesel = dt.datetime(ptime[0], ptime[1], ptime[2])
    start_date = timesel
    gdf_ps['tmax'] = clim.tmax.sel(time=timesel)
    gdf_ps['tmin'] = clim.tmin.sel(time=timesel)
    gdf_ps['prcp'] = clim.prcp.sel(time=timesel)
    fig, ax = plt.subplots(ncols=3)
    divider0 = make_axes_locatable(ax[0])
    divider1 = make_axes_locatable(ax[1])
    divider2 = make_axes_locatable(ax[2])
    cax0 = divider0.append_axes("right", size="5%", pad=0.1)
    cax1 = divider1.append_axes("right", size="5%", pad=0.1)
    cax2 = divider2.append_axes("right", size="5%", pad=0.1)
    h_tmax = gdf_ps.tmax.max()
    l_tmax = gdf_ps.tmax.min()
    h_tmin= gdf_ps.tmin.max()
    l_tmin= gdf_ps.tmin.min()
    h_tmax = gdf_ps.tmax.max()
    l_tmax = gdf_ps.tmax.min()
    h_ppt= gdf_ps.prcp.max()
    l_ppt= gdf_ps.prcp.min()

    gdf_ps.plot(column='tmax', ax=ax[0], vmin=l_tmax, vmax=h_tmax, legend=True,
                label='tmax',cax=cax0)
    gdf_ps.plot(column='tmin', ax=ax[1], vmin=l_tmin, vmax=l_tmin, legend=True, 
                label='tmin',cax=cax1)
    gdf_ps.plot(column='prcp', ax=ax[2], vmin=l_ppt, vmax=l_ppt, legend=True,
                label='prcp',cax=cax2)
    for i in range(3):
        ax[i].set_xticklabels([])
        ax[i].set_yticklabels([])
        if i == 0:
            ax[i].set_title('tmax')
        elif i == 1:
            ax[i].set_title('tmin')
        elif i == 2:
            ax[i].set_title('prcp')
    plt.tight_layout()
    
    return clim


def example_plot(clim, gdf_ps, msurf, msoil, j, timesel):
    gdf_ps['tmax'] = clim.tmax.sel(time=timesel)
    gdf_ps['tmin'] = clim.tmin.sel(time=timesel)
    gdf_ps['prcp'] = clim.prcp.sel(time=timesel)

    gdf_ps['infil'] = msurf.var['infil'].data
    gdf_ps['sroff'] = msurf.var['sroff'].data
    gdf_ps['soil_moist_tot'] = msoil.var['soil_moist_tot'].data

    fig, ax = plt.subplots(ncols=6, figsize = (12,2))
    divider0 = make_axes_locatable(ax[0])
    divider1 = make_axes_locatable(ax[1])
    divider2 = make_axes_locatable(ax[2])
    divider3 = make_axes_locatable(ax[3])
    divider4 = make_axes_locatable(ax[4])
    divider5 = make_axes_locatable(ax[5])
    cax0 = divider0.append_axes("right", size="5%", pad=0.1)
    cax1 = divider1.append_axes("right", size="5%", pad=0.1)
    cax2 = divider2.append_axes("right", size="5%", pad=0.1)
    cax3 = divider3.append_axes("right", size="5%", pad=0.1)
    cax4 = divider4.append_axes("right", size="5%", pad=0.1)
    cax5 = divider5.append_axes("right", size="5%", pad=0.1)
    
    gdf_ps.plot(column='tmax', vmin=20.0, vmax=65.0, ax=ax[0], legend=True, cax=cax0)
    gdf_ps.plot(column='tmin', vmin=20.0, vmax=65.0, ax=ax[1], legend=True, cax=cax1)
    gdf_ps.plot(column='prcp', vmin=0.0, vmax=0.7, ax=ax[2], legend=True, cax=cax2)
    gdf_ps.plot(column='infil', vmin=0.0, vmax=0.7, ax=ax[3], legend=True, cax=cax3)
    gdf_ps.plot(column='sroff', vmin=0.0, vmax=0.25, ax=ax[4], legend=True, cax=cax4)
    gdf_ps.plot(column='soil_moist_tot', vmin=0.25, vmax=1.75, ax=ax[5], legend=True, cax=cax5)
    for i in range(6):
        ax[i].set_xticklabels([])
        ax[i].set_yticklabels([])
        if j == 0:
            if i == 0:
                ax[i].set_title('tmax')
            elif i == 1:
                ax[i].set_title('tmin')
            elif i == 2:
                ax[i].set_title('prcp')
            elif i == 3:
                ax[i].set_title('infil')
            elif i == 4:
                ax[i].set_title('sroff')
            elif i == 5:
                ax[i].set_title('soil_moist_tot')
    plt.tight_layout()
    
def gm_example_plot(gdf_ps, gmdata, msurf, msoil, j, timesel):
    gdf_ps['tmax'] = (gmdata.tmax.data[j,:]*(9/5))+32.0
    gdf_ps['tmin'] = (gmdata.tmin.data[j,:]*(9/5))+32.0
    gdf_ps['prcp'] = gmdata.precip.data[j,:]*.0393701
#     print(gmdata.precip[j,:]*.0393701)
#     print(msurf.var['tmax'].data)
#     print(msurf.var['tmin'].data)
#     print(msurf.get_value('hru_ppt'))

    gdf_ps['infil'] = msurf.var['infil'].data
    gdf_ps['sroff'] = msurf.var['sroff'].data
    gdf_ps['soil_moist_tot'] = msoil.var['soil_moist_tot'].data

    fig, ax = plt.subplots(ncols=6, figsize = (12,2))
    divider0 = make_axes_locatable(ax[0])
    divider1 = make_axes_locatable(ax[1])
    divider2 = make_axes_locatable(ax[2])
    divider3 = make_axes_locatable(ax[3])
    divider4 = make_axes_locatable(ax[4])
    divider5 = make_axes_locatable(ax[5])
#     divider6 = make_axes_locatable(ax[6])
    cax0 = divider0.append_axes("right", size="5%", pad=0.1)
    cax1 = divider1.append_axes("right", size="5%", pad=0.1)
    cax2 = divider2.append_axes("right", size="5%", pad=0.1)
    cax3 = divider3.append_axes("right", size="5%", pad=0.1)
    cax4 = divider4.append_axes("right", size="5%", pad=0.1)
    cax5 = divider5.append_axes("right", size="5%", pad=0.1)
#     cax6 = divider6.append_axes("right", size="5%", pad=0.1)
    
    gdf_ps.plot(column='tmax', vmin=50.0, vmax=70.0, ax=ax[0], legend=True, cax=cax0)
    gdf_ps.plot(column='tmin', vmin=20.0, vmax=45.0, ax=ax[1], legend=True, cax=cax1)
    gdf_ps.plot(column='prcp', vmin=0.0, vmax=.75, ax=ax[2], legend=True, cax=cax2)
    gdf_ps.plot(column='infil', vmin=0.0, vmax=0.7, ax=ax[3], legend=True, cax=cax3)
    gdf_ps.plot(column='sroff', vmin=0.0, vmax=0.25, ax=ax[4], legend=True, cax=cax4)
    gdf_ps.plot(column='soil_moist_tot', vmin=0.25, vmax=1.75, ax=ax[5], legend=True, cax=cax5)
#     gdf_ps.plot(column='soil_moist_tot', vmin=0.0, vmax=1.5, ax=ax[6], legend=True, cax=cax6)
    for i in range(6):
        ax[i].set_xticklabels([])
        ax[i].set_yticklabels([])
        if j == 0:
            if i == 0:
                ax[i].set_title('tmax')
            elif i == 1:
                ax[i].set_title('tmin')
            elif i == 2:
                ax[i].set_title('prcp')
            elif i == 3:
                ax[i].set_title('soil_to_gw')
            elif i == 4:
                ax[i].set_title('ssr_to_gw')
            elif i == 5:
                ax[i].set_title('soil_moist_tot')
#             elif i == 6:
#                 ax[i].set_title('soil_moist_tot')
    plt.tight_layout()
    
def example_plot2(gdf_ps, msurf, msoil, j, timesel):
    gdf_ps['tmax'] = (msurf.get_value('tmax')*(9.0/5.0)) + 32.0
    gdf_ps['tmin'] = (msurf.get_value('tmin')*(9.0/5.0)) + 32.0
    gdf_ps['prcp'] = msurf.get_value('hru_ppt')

    gdf_ps['soil_to_gw'] = msoil.var['soil_to_gw'].data
    gdf_ps['ssr_to_gw'] = msoil.var['ssr_to_gw'].data
    gdf_ps['ssres_flow'] = msoil.var['ssres_flow'].data
    gdf_ps['soil_moist_tot'] = msoil.var['soil_moist_tot'].data

    fig, ax = plt.subplots(ncols=7, figsize = (14,2))
    divider0 = make_axes_locatable(ax[0])
    divider1 = make_axes_locatable(ax[1])
    divider2 = make_axes_locatable(ax[2])
    divider3 = make_axes_locatable(ax[3])
    divider4 = make_axes_locatable(ax[4])
    divider5 = make_axes_locatable(ax[5])
    divider6 = make_axes_locatable(ax[6])
    cax0 = divider0.append_axes("right", size="5%", pad=0.1)
    cax1 = divider1.append_axes("right", size="5%", pad=0.1)
    cax2 = divider2.append_axes("right", size="5%", pad=0.1)
    cax3 = divider3.append_axes("right", size="5%", pad=0.1)
    cax4 = divider4.append_axes("right", size="5%", pad=0.1)
    cax5 = divider5.append_axes("right", size="5%", pad=0.1)
    cax6 = divider6.append_axes("right", size="5%", pad=0.1)
    
    gdf_ps.plot(column='tmax', vmin=20.0, vmax=75.0, ax=ax[0], legend=True, cax=cax0)
    gdf_ps.plot(column='tmin', vmin=20.0, vmax=75.0, ax=ax[1], legend=True, cax=cax1)
    gdf_ps.plot(column='prcp', vmin=0.0, vmax=0.7, ax=ax[2], legend=True, cax=cax2)
    gdf_ps.plot(column='soil_to_gw', vmin=0.0, vmax=0.1, ax=ax[3], legend=True, cax=cax3)
    gdf_ps.plot(column='ssr_to_gw', vmin=0.0, vmax=0.15, ax=ax[4], legend=True, cax=cax4)
    gdf_ps.plot(column='ssres_flow', vmin=0.0, vmax=0.1, ax=ax[5], legend=True, cax=cax5)
    gdf_ps.plot(column='soil_moist_tot', vmin=0.25, vmax=3.0, ax=ax[6], legend=True, cax=cax6)
    for i in range(7):
        ax[i].set_xticklabels([])
        ax[i].set_yticklabels([])
        if j == 0:
            if i == 0:
                ax[i].set_title('tmax')
            elif i == 1:
                ax[i].set_title('tmin')
            elif i == 2:
                ax[i].set_title('prcp')
            elif i == 3:
                ax[i].set_title('soil_to_gw')
            elif i == 4:
                ax[i].set_title('ssr_to_gw')
            elif i == 5:
                ax[i].set_title('ssres_flow')
            elif i == 6:
                ax[i].set_title('soil_moist_tot')
    plt.tight_layout()