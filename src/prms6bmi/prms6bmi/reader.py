"""
Created on Thu Dec 12 08:00:48 2019

@author:rmcd build on pangeo package by Steve Markstrom - USGS
"""

import xarray as xr
import glob
import os
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def get_DataSet_prms6(summary, myparam):
    # merge spatial locations of hru and segments into summary file
    ds = xr.open_dataset(summary)
    param = xr.open_dataset(myparam)
    hru_lat = param.get("hru_lat")
    ds['hru_lat'] = hru_lat
    hru_lon = param.get("hru_lon")
    ds['hru_lon'] = hru_lon
    seg_lat = param.get("seg_lat")
    ds['seg_lat'] = seg_lat
    seg_lon = param.get("seg_lon")
    ds['seg_lon'] = seg_lon
    return ds


def bmi_prms6_value_splot(gdf, mbmi, value, tvmin, tvmax, index, timesel, pax = None):
    tax = pax or plt.gca()
    gdf[value] = mbmi.get_value(value)
    divider = make_axes_locatable(tax)
    tcax = divider.append_axes(position='right', size='5%', pad=0.1)
    gdf.plot(column=value, vmin=tvmin, vmax=tvmax, ax=tax, legend=True, cax=tcax)
    tax.set_title(value)


def plot_climate(c_xarray, hru_index, val, start, end, tax=None):
    tax = tax or plt.gca()
    hru_ids = c_xarray.hru.values
    simclimate = c_xarray.sel(time=slice(start, end))

    line, = simclimate.sel(hru=hru_ids[hru_index])[val].plot(ax=tax)
    tax.set_title(val)

def bmi_prms6_value_plot(data, n_index, val, label, start, end, tax = None):
    tax = tax or plt.gca()
    #test if val exists in both and get nhru or nsegment
    dim_type = None
    try:
        dim_type = data[val].dims[1]

        if dim_type == 'nhru':
            data_val = data[val].sel(nhru=n_index, time=slice(start, end)).to_pandas()
            # dprms_val = dprms[val].sel(nhru=n_index, time=slice(start, end))
            data_val.plot.line(ax=tax, label=label)
            tax.legend()
            # line1, = dprms_val.plot.line(x='time', ax=tax, add_legend=True)

        elif dim_type == 'nsegment':
            data_val = data[val].sel(nsegment=n_index, time=slice(start, end)).to_pandas()
            # dprms_val = dprms[val].sel(nsegment=n_index, time=slice(start, end)).to_pandas()

            data_val.plot(ax=tax, label=label)
            tax.legend()
            # line1, = dprms_val.plot(label='PRMS6')

        tax.set_title(f'{val} {n_index}')

    except Exception as err:
        print('Error', {err})

def bmi_prms6_residual_plot(dbmi, dprms, n_index, val, label, start, end, tax = None):
    tax = tax or plt.gca()
    dim_type = dbmi[val].dims[1]
    try:
        if dim_type == 'nhru':
            data_val = dbmi[val] - dprms[val]
            data = data_val.sel(nhru=n_index, time=slice(start, end)).to_pandas()
            # bmi = dbmi[val]
            # prms = dprms.sel(nhru=n_index, time=slice(start, end))[val]
        elif dim_type == 'nsegment':
            data_val = dbmi[val] - dprms[val]
            data = data_val.sel(nsegment=n_index, time=slice(start, end)).to_pandas()
            # bmi = dbmi.sel[val]
            # prms = dprms.sel(nsegment=n_index, time=slice(start, end))[val]
        # res = prms-bmi
        data.plot(ax=tax, label=label)
        plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
        tax.legend()
        tax.set_title('Residual (prms-bmi)')
    except Exception as err:
        print('Error', {err})

def get_feat_coord(feat, data_set, feat_id):
    lat_da = data_set[feat + '_lat']
    lat = lat_da[feat_id-1].values
    lon_da = data_set[feat + '_lon']
    lon = lon_da[feat_id-1].values
    return lat,lon


def get_hrus_for_box(ds, lat_min, lat_max, lon_min, lon_max):
    sel = ds.hru_lat.sel(hruid=((ds.hru_lat.values >= lat_min)
                            & (ds.hru_lat.values <= lat_max)))
    ids_1 = sel.hruid.values
    sel_1 = ds.hru_lon.sel(hruid=ids_1)
    sel_2 = sel_1.sel(hruid=((sel_1.values >= lon_min) & (sel_1.values <= lon_max)))
    ids_2 = sel_2.hruid.values
    return ids_2


def get_segs_for_box(ds, lat_min, lat_max, lon_min, lon_max):
    sel = ds.seg_lat.sel(segid=((ds.seg_lat.values >= lat_min)
                            & (ds.seg_lat.values <= lat_max)))
    ids_1 = sel.segid.values
    sel_1 = ds.seg_lon.sel(segid=ids_1)
    sel_2 = sel_1.sel(segid=((sel_1.values >= lon_min) & (sel_1.values <= lon_max)))
    ids_2 = sel_2.segid.values
    return ids_2


def get_values_for_DOY(ds, timestamp, hru_ids, var_name):
    if (timestamp < pd.Timestamp('1979-10-01') or timestamp > pd.Timestamp('1980-09-30')):
        print("The date you provided is outside of range 1979-10-01 to 1980-09-30")
        return None
        
    time_range = pd.date_range(timestamp, freq='1Y', periods=40)
    dif = timestamp - time_range[0]
    time_range = time_range + dif
    # print(time_range)

    date_list = []
    val_list = []
    for ts in time_range:
        try:
            date_str = str(ts.year).zfill(4) + '-' + str(ts.month).zfill(2) + '-' + str(ts.day).zfill(2)
            ds_sel = ds[var_name].sel(hruid=hru_ids, time=date_str)
            val = ds_sel.values[0][0]
            date_list.append(date_str + 'T05:00:00')
            val_list.append(val)
        except:
            pass
        
    val_np = np.asarray(val_list, dtype=np.float64)
    val_np = val_np.reshape((1, val_np.shape[0]))
    hru_ids_np = np.asarray(hru_ids, dtype=np.int32)
    date_np = np.asarray(date_list, dtype='datetime64[ns]')
    
    attrs = ds[var_name].attrs
    da_new = xr.DataArray(data=val_np, dims=['hruid','time'],
                          coords={'hruid':hru_ids_np,'time':date_np},
                          attrs=attrs)

    return da_new
    
