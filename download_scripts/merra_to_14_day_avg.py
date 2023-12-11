import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import os
from scipy.interpolate import interp2d
import glob

res = 9;
res = 36;
in_dir = '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/data_global/merra_forcing/ease2_'+str(res)+'km/'


file_list = glob.glob(in_dir+'*.nc')
file_list.sort()
def open_1day(filename):
    ds = xr.open_dataset(filename);
    return ds

avg_window_size = 15;
range = int(avg_window_size/2.);

print len(file_list)
print ''
for i in np.arange(0,len(file_list)):
    print file_list[i];
    if i<range:
        start_idx = 0;
    else:
        start_idx = i-range;
    if i>len(file_list)-range:
        end_idx = len(file_list);
    else:
        end_idx = i+range;
    ds = open_1day(file_list[i]);

    time_dim = end_idx - start_idx + 1;
    temp_data = np.ones((ds.temperature_mean.shape[0],ds.temperature_mean.shape[1],time_dim)); temp_data[:,:,:]=np.nan;
    ea_data = np.ones((ds.temperature_mean.shape[0],ds.temperature_mean.shape[1],time_dim)); ea_data[:,:,:]=np.nan;
    k = 0;
    for j in np.arange(start_idx,end_idx):
        ds_i = open_1day(file_list[j]);
        temp_data[:,:,k]=ds_i.temperature_max.data;
        ea_data[:,:,k]=ds_i.ea2m_daily_mean.data;
        k+=1;

    ds['Ea_'+str(avg_window_size)+'_day_mean']   = (('lat', 'lon'), np.nanmean(ea_data,2))
    ds['Tmax_'+str(avg_window_size)+'_day_mean'] = (('lat', 'lon'), np.nanmean(temp_data,2))
    os.system('rm ' +file_list[i])
    ds.to_netcdf(file_list[i])
    print 'file overwritten '+file_list[i]
