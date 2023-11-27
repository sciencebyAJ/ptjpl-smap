# _run_ptjpl_smap_global_.py
# Global PTJPL SMAP ET 
__author__ = 'AJ Purdy'
# Last updated: May 2017
# -------------------------------------------------------------
# FORCING DATASETS:
# SMAP DATA:
#  * SMAP L3 P36 [cm3/cm3]
#  * SMAP L3 P_E 9km [cm3/cm3]
#  * SMAP L4 soil properties: wilting point:
#  * SMAP L4 soil properties: porosity
# MERRA DATA:
#  * Radiation [W/m2]
#  * Air Temperature [K]
#  * Vapor pressure [Pa]
# MODIS DATA:
#  * NDVI-- MOD13 & MYD13 0.05o
# CANOPY HEIGHT:
#  * Simmard et al., 2011: [m]
# LANDCOVER:
#  * MODIS LANDCOVER PARAMETERS RESAMPLED FOR SMAP: SMAP ATBD
# 
# TO UPDATE THE ABOVE RUN THESE SCRIPTS
# -------------------------------------------------------------
from netCDF4 import Dataset
import xarray as xr
import json
from numpy import genfromtxt
import numpy as np
import os
import sys
import datetime
# -------------------------------------------------------------
# Move this into PTJPL_LIB
# os.chdir('/Volumes/AJ_RESEARCH/SMAP_ET_1_tower_eval/SCRIPTS/')
# Determine which functions I am grabbing from each of these:
# import FLUXNET_META
# import SM_ET_lib
# from fluxnet_diag_plots import *
from ptjpl_lib.ptjpl_smap_lib import *
# -------------------------------------------------------------
# set home directory to navigate from data folders to this folder:

home_dir = os.getcwd();
data_dir = '/data/data_global/'
data_out_dir = '/results/results_global/'

# -------------------------------------------------------------
# Get command line arguments
# -------------------------------------------------------------

IS_arg=len(sys.argv)
if IS_arg <4 or IS_arg>5 :
     print('ERROR - Expecting 3 inputs: \n\t1.\tResolution 3 9 or 36\n\t2.\tYear [format YYYY]\n\t3.\tStart DOY [format DOY]\n\t4.\End DOY [format DOY]')
     raise SystemExit(22) 
res_str=sys.argv[1]
res = int(res_str)
year = sys.argv[2]
start_doy = sys.argv[3]
end_doy   = sys.argv[4]
# -------------------------------------------------------------
# Check for resolution
# -------------------------------------------------------------
if res!=9 and res!=36:
    print('ERROR - resolution needs to be: 9, or 36')
    raise SystemExit(22)


# -------------------------------------------------------------
# Print input information
# -------------------------------------------------------------
print('\nCommand line inputs:')
print('\tresolution set to \t '+ str(res))
print('\tyear \t\t\t '         + str(year))
print('\tdoy start \t\t '      + start_doy)
print('\tdoy end \t\t '        + end_doy)

# -------------------------------------------------------------
with open('ptjpl_smap_forcing_meta_data_gp.json') as data_file:    
    meta_data = json.load(data_file)


# -------------------------------------------------------------    
# Data for original runs:
# for key in meta_data.keys(): print key

ndvi_path     = meta_data['NDVI']['path_'+str(res)+'km']
rnet_path     = meta_data['NET_RADIATION']['path_'+str(res)+'km']
tmax_path     = meta_data['AIR_TEMPERATURE']['path_'+str(res)+'km']
tmean_path    = meta_data['AIR_TEMPERATURE_MEAN']['path_'+str(res)+'km']
vp_path       = meta_data['VAPOR_PRESSURE']['path_'+str(res)+'km']
faparmax_path = meta_data['fAPAR_MAX']['path_'+str(res)+'km']
topt_path     = meta_data['T_OPT']['path_'+str(res)+'km']

# -------------------------------------------------------------
# New data added to algorithm:

smap_path     = meta_data['SOIL_MOISTURE_'+str(res)+'KM']['path_'+str(res)+'km']
ch_path       = meta_data['CANOPY_HEIGHT']['path_'+str(res)+'km']
wp_path       = meta_data['WILTING_POINT']['path_'+str(res)+'km']
por_path       = meta_data['POROSITY']['path_'+str(res)+'km']

# -------------------------------------------------------------
# New data added to algorithm: Add these to metadata file:
ndvi_days = ['001','009','017','25','033','41','49','57','065','73','081','89','097',
             '105','113','121','129','137','145','153','161','169','177','185','193',
             '201','209','217','225','233','241','249','257','265','273','281','289',
             '297','305','313','321','329','337','345','353','361']; 
            
ndvi_days = list(map(int, ndvi_days))

# -------------------------------------------------------------

# -------------------------------------------------------------

def read_soil_props():
    fc_file = data_dir+por_path+meta_data['POROSITY']['file_start'] + meta_data['POROSITY'][str(res)] + meta_data['POROSITY']['file_end']
    fh = Dataset(fc_file, mode='r')
    fc = fh.variables[meta_data['POROSITY']['varname']][:]
    fh.close()
    wp_file = data_dir+wp_path+meta_data['WILTING_POINT']['file_start'] + meta_data['WILTING_POINT'][str(res)] + meta_data['WILTING_POINT']['file_end']
    fh = Dataset(wp_file, mode='r')
    wp = fh.variables[meta_data['WILTING_POINT']['varname']][:]
    fh.close()
    return wp, fc

def read_CH():
    ch_file = data_dir+ch_path + meta_data['CANOPY_HEIGHT']['file_start'] + meta_data['CANOPY_HEIGHT'][str(res)] + meta_data['CANOPY_HEIGHT']['file_end']
    fh = Dataset(ch_file, mode='r')
    canopy_height = fh.variables[meta_data['CANOPY_HEIGHT']['varname']][:]
    fh.close()
    return canopy_height.T

def load_Topt_fAPARmax():
    # load optimum temperature (C)
    Topt_file = data_dir+topt_path+ meta_data['T_OPT']['file_start'] + meta_data['T_OPT'][str(res)] + meta_data['T_OPT']['file_end']
    fh = Dataset(Topt_file, mode='r')
    optimum_temperature = fh.variables[meta_data['T_OPT']['varname']][:]
    fh.close()
    
    fAPARmax_file = data_dir+faparmax_path+ meta_data['fAPAR_MAX']['file_start'] + meta_data['fAPAR_MAX'][str(res)] + meta_data['fAPAR_MAX']['file_end']
    fh = Dataset(fAPARmax_file, mode='r')
    fAPARmax = fh.variables[meta_data['fAPAR_MAX']['varname']][:]
    fh.close()
    return optimum_temperature, fAPARmax.T


def read_SMAP_data(smap_time,resolution): # change to L4_
    smap_key='SOIL_MOISTURE_'+str(resolution)+'KM'
    smap_file = data_dir+smap_path + meta_data[smap_key]['file_start'] + smap_time + meta_data[smap_key]['file_end']
    fh = Dataset(smap_file, mode='r')
    smap_sm = fh.variables[meta_data[smap_key]['varname']][:]
    fh.close()
    return smap_sm

def read_ndvi(time):
    ndvi_file = data_dir+ndvi_path +meta_data['NDVI']['file_start'] + time + meta_data['NDVI']['file_end_'+str(res)+'km']
    fh = Dataset(ndvi_file, mode='r')
    NDVI_array = fh.variables[meta_data['NDVI']['varname']][:]
    fh.close()
    return NDVI_array.T;

def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return array[idx]

def file_name_times(year,doy):
    ndvi_DOY = str(find_nearest(ndvi_days,doy)).zfill(3);
    ndvi_time = year + ndvi_DOY;
    DOY =  str(doy).zfill(3)
    date = year+'-'+DOY
    date = datetime.datetime.strptime(date,'%Y-%j')
    met_time = date.strftime('%Y%m%d')
    smap_time = date.strftime('%Y%m%d')
    return ndvi_time,met_time,smap_time

def read_met_forcing(time):
    '''
    input:   time in YYYY_DOY (%Y_%j)
    outputs: Tmax (K), Tmean (K), Vapor Pressure (Pa), Rnet (W/m2)
    '''
    ##read max temp data
    merra_file = data_dir+tmax_path+meta_data['AIR_TEMPERATURE']['file_start']+time+ meta_data['AIR_TEMPERATURE']['file_end']
    fh = Dataset(merra_file, mode='r')
    tmax_K = fh.variables[meta_data['AIR_TEMPERATURE']['varname']][:]
    ##read mean temp data 
    tmean_K = fh.variables[meta_data['AIR_TEMPERATURE_MEAN']['varname']][:]
    ##read mean temp data for 15 day avg
    tmax_K_15_day = fh.variables[meta_data['AIR_TEMPERATURE_15_day']['varname']][:]
    ##read radiation data
    rnet_Wm2 = fh.variables[meta_data['NET_RADIATION']['varname']][:]
    ##read vapor pressure
    vp_Pa = fh.variables[meta_data['VAPOR_PRESSURE']['varname']][:]
    ##read vapor pressure for 15 day avg
    vp_Pa_15_day = fh.variables[meta_data['VAPOR_PRESSURE_15_days']['varname']][:]
    
    fh.close()
    ## return values
    return tmax_K, tmean_K, vp_Pa, rnet_Wm2, tmax_K_15_day, vp_Pa_15_day

# -------------------------------------------------------------
# LOAD STATIC DATA
# -------------------------------------------------------------
print('\nStatic variables loaded:')
if res == 9:
    lon_fname = 'SMAP_L4_LON_1d_global.csv'
    lat_fname = 'SMAP_L4_LAT_1d_global.csv'
elif res == 36:
    lon_fname = 'SMAP_L3P_LON_1d_36km_global.csv'
    lat_fname = 'SMAP_L3P_LAT_1d_36km_global.csv'
else:
    print("please enter 3, 9 or 36")
# Coordinate path & filenames: #<--- consider adding these to json metadata file
lon_lat_path = data_dir+'LAT_LON/'
lons1d = genfromtxt(lon_lat_path+lon_fname, delimiter=',')
lats1d = genfromtxt(lon_lat_path+lat_fname, delimiter=',')
lon_0 = lons1d.mean();
lat_0 = lats1d.mean();
lon, lat = np.meshgrid(lons1d, lats1d)
print('\tlat lon')
# optimum temperature
optimum_temperature,fAPARmax = load_Topt_fAPARmax()
print('\toptimum_temperature')
print('\tfAPARmax')
#wilting_point & field capacity
wilting_point, field_capacity = read_soil_props()
print('\twilting_point')
print('\tporosity')
# canopy height
CHx = read_CH()
mask = (wilting_point>0.0);
canopy_height = CHx
canopy_height[~mask] =np.nan;
print('\tcanopy_height')

# -------------------------------------------------------------
print('\nstarting model run: start='+year+start_doy.zfill(3)+' end='+year+end_doy.zfill(3))
missing_data=[];
for doy in np.arange(int(start_doy),int(end_doy)): # 91 is the start date 4/1/2015
    ndvi_time, merra_time, smap_time = file_name_times(str(year),doy)
    try:
        print('\t'+year+'-'+str(doy).zfill(3))
        print('\t\tLOADING FORCING DATA:')
        air_temperature_K, air_temperature_mean_K, water_vapor_pressure_mean_Pa, net_radiation, air_temperature_max_K_15day,water_vapor_pressure_mean_Pa_15day = read_met_forcing(merra_time)
        print('\t\t\tmerra data loaded')
        ndvi_mean      = read_ndvi(ndvi_time);
        print('\t\t\tndvi data loaded')
        soil_moisture  = read_SMAP_data(smap_time,int(res))
        print('\t\t\tsmap data loaded')
        print('\t\tEXECUTING MODEL')
        results_i = ptjpl_smap(lat,
                               lon,
                               air_temperature_K,
                               air_temperature_mean_K,
                               water_vapor_pressure_mean_Pa,
                               air_temperature_max_K_15day,
                               water_vapor_pressure_mean_Pa_15day,
                               ndvi_mean,
                               optimum_temperature,
                               fAPARmax,
                               net_radiation,
                               wilting_point,
                               field_capacity,
                               canopy_height,
                               soil_moisture,
                               verbose=True,
                               );
    
        results_i.to_netcdf(data_out_dir+str(res)+'km/'+'PTJPL_SMAP_ET_'+smap_time+'_'+str(res)+'km.nc');
        print('\t\tFILE CREATED: PTJPL_SMAP_ET_'+smap_time+'_'+str(res)+'km.nc')
    except:
        print('\t**************************file missing for ' + year + '-' + str(doy).zfill(3)+'**************************')
        missing_data.append(smap_time)
        continue
print('ALL RESULTS SAVED TO '+data_out_dir+str(res)+'km/n')


