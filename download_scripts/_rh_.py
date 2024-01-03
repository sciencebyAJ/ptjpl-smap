import datetime
import glob
import json
from netCDF4 import Dataset
import numpy as np
from numpy import exp, log, genfromtxt
import os
import sys
import xarray as xr

home_dir = os.getcwd();
data_dir = '/data15/famiglietti2/ajpurdy/SMAP_ET_global/data/'
rh_out_dir_s = '/data15/famiglietti2/ajpurdy/SMAP_ET_global/data/'

IS_arg=len(sys.argv)
if IS_arg <4 or IS_arg>5 :
     print('ERROR - Expecting 3 inputs: \n\t1.\tResolution 3 9 or 36\n\t2.\tYear [format YYYY]\n\t3.\tStart DOY [format DOY]\n\t4.\End DOY [format DOY]')
     raise SystemExit(22)
res_str=sys.argv[1]
res = int(res_str)
year = sys.argv[2]
start_doy = sys.argv[3]
end_doy   = sys.argv[4]
if res!=3 and res!=9 and res!=36:
    print('ERROR - resolution needs to be: 9, or 36')
    raise SystemExit(22)

    print('\nCommand line inputs:')

    
print('\tresolution set to \t '+ str(res))
rh_out_dir = rh_out_dir_s+'EASE_'+str(res)+'/RH/'
print('\tyear \t\t\t '         + str(year))
print('\tdoy start \t\t '      + start_doy)
print('\tdoy end \t\t '        + end_doy)
with open('ptjpl_smap_forcing_meta_data_gp.json') as data_file:
    meta_data = json.load(data_file)

tmean_path    = meta_data['AIR_TEMPERATURE_MEAN']['path_'+str(res)+'km']
vp_path       = meta_data['VAPOR_PRESSURE']['path_'+str(res)+'km']


# -------------------------------------------------------------
def read_met_forcing(time):
    '''
    input:   time in YYYY_DOY (%Y_%j)
    outputs: Tmax (K), Tmean (K), Vapor Pressure (Pa), Rnet (W/m2)
    '''
    ##read max temp data
    merra_file = data_dir+tmean_path+meta_data['AIR_TEMPERATURE']['file_start']+time+ meta_data['AIR_TEMPERATURE']['file_end']
    #print(merra_file)
    fh = Dataset(merra_file, mode='r')
    ##read mean temp data 
    tmean_K = fh.variables[meta_data['AIR_TEMPERATURE_MEAN']['varname']][:]
    print('\t\t\ttmean')
    ##read vapor pressure
    vp_Pa = fh.variables[meta_data['VAPOR_PRESSURE']['varname']][:]
    print('\t\t\tvp')
    ##read vapor pressure for 15 day avg
    fh.close()
    return tmean_K, vp_Pa

# -------------------------------------------------------------
def EASE2_LON_LAT(res):
    directory_start = '/data15/famiglietti2/ajpurdy/SMAP_ET_global/data/LAT_LON/EASE2/';
    folder_path = 'gridloc_EASE2_M'+str(res).zfill(2)+'km';
    fnameLON = glob.glob(directory_start+folder_path+'/*lons*')[0]
    fnameLAT = glob.glob(directory_start+folder_path+'/*lats*')[0]
    scale = int(res/3.)
    ny3 = int(4872/scale)
    nx3 = int(11568/scale)
    LAT = np.fromfile(fnameLAT, dtype=np.double).reshape((ny3,nx3))
    LON = np.fromfile(fnameLON, dtype=np.double).reshape((ny3,nx3))
    return LON[0,:], LAT[:,0]
# -------------------------------------------------------------
def file_name_times(year,doy):
    DOY =  str(doy).zfill(3)
    date = year+'-'+DOY
    date = datetime.datetime.strptime(date,'%Y-%j')
    met_time = date.strftime('%Y%m%d')
    smap_time = date.strftime('%Y%m%d')
    return met_time,smap_time
# -------------------------------------------------------------
def saturation_vapor_pressure_from_air_temperature(air_temperature):
    SVP_BASE = 0.611
    SVP_MULT = 17.27
    SVP_ADD = 237.7
    svp = SVP_BASE * exp((air_temperature * SVP_MULT) / (air_temperature + SVP_ADD))
    return svp
# -------------------------------------------------------------
def rh_calculator(lat,lon, air_temperature_mean_K, water_vapor_pressure_mean_Pa):
    air_temperature_mean = air_temperature_mean_K - 273
    saturation_vapor_pressure = saturation_vapor_pressure_from_air_temperature(air_temperature_mean)
    saturation_vapor_pressure = np.clip(saturation_vapor_pressure, 1, None)
    water_vapor_pressure_mean = water_vapor_pressure_mean_Pa * 0.001
    vapor_pressure_deficit = saturation_vapor_pressure - water_vapor_pressure_mean
    relative_humidity = water_vapor_pressure_mean / saturation_vapor_pressure
    relative_humidity[(relative_humidity<0.)|(relative_humidity>1.)]=np.nan;
    results = xr.Dataset({'relative_humidity': (['x', 'y'],           relative_humidity)},
                         coords={'lon': (['x', 'y'], lon),'lat': (['x', 'y'], lat)})
    return results
# -------------------------------------------------------------
lon_lat_path = data_dir+'LAT_LON/'
lons1d, lats1d = EASE2_LON_LAT(res);
lon, lat = np.meshgrid(lons1d, lats1d)
# -------------------------------------------------------------
print('\ncalculating RH and saving new files for: start='+year+start_doy.zfill(3)+' end='+year+end_doy.zfill(3))
missing_data=[];

for doy in np.arange(int(start_doy),int(end_doy)): # 91 is the start date 4/1/2015
    merra_time, smap_time = file_name_times(str(year),doy)
    try:
        print('\tSTARTING TO RUN MODEL ON:'+str(year)+str(doy))
        print('\t\tLOADING FORCING DATA:')
        #print(merra_time)
        air_temperature_mean_K, water_vapor_pressure_mean_Pa= read_met_forcing(merra_time)
        print('\t\t\tmerra data loaded')
        print('\t\tMAKING RH FILE')
        results_i = rh_calculator(lat,lon,air_temperature_mean_K, water_vapor_pressure_mean_Pa);
        results_i.to_netcdf(rh_out_dir+'RH_'+smap_time+'_'+str(res)+'km.nc');
        print('\tFILE CREATED: RH_'+smap_time+'_'+str(res)+'km.nc')
    except:
        print('\t**************************file missing for ' + year + '-' + str(doy).zfill(3)+'**************************')
        missing_data.append(smap_time)
        continue
print('\nALL RESULTS SAVED TO '+rh_out_dir+' ')

