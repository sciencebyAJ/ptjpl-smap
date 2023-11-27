import numpy as np
import xarray as xr
import os
import sys
from scipy.interpolate import griddata

IS_arg=len(sys.argv)
if IS_arg <3 or IS_arg>4 :
    print('ERROR - Expecting 3 inputs: \n\t1.\tYear [format YYYY]\n\t2.\tStart Month [format mm]\n\t3.\tEnd Month [format mm]')
    raise SystemExit(22)
year_in = sys.argv[1]
start_month = int(sys.argv[2])
end_month   = int(sys.argv[3])

if start_month == end_month:
    end_month = end_month + 1

years = [str(year_in)]
#years = ['2017']#,'2015','2016']

month_arr = np.arange(start_month, end_month)
month_list = month_arr.tolist()
months=[str(item).zfill(2) for item in month_list]
print(months)

# for local use
path_root = '/Volumes/AJ_4TB/Coding/Research/'
# for h2observer use
path_root = '/raid/apurdy/projects/global_et/'


def create_MERRA_filenames(YYYY,MM,DD):
    '''creates paths to opendap access to MERRA2 datasets for temperature and radiation
    returns filenameTEMP for access to daily air temperature, surface pressure and specific humidity
    returns filenameRAD for access to daily short and longwave net radiation
    '''
#    TEMPfile_start = 'https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2I1NXASM.5.12.4/'
    TEMPfile_start = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I1NXASM.5.12.4/'
    #2016/12/MERRA2_400.inst1_2d_asm_Nx.20161201.nc4'

    TEMPfile_end   = '/MERRA2_400.inst1_2d_asm_Nx.'+YYYY+MM+DD+'.nc4'
    filenameTEMP   = TEMPfile_start+YYYY+'/'+MM+TEMPfile_end;

#    RADfile_start = 'https://goldsmr4.gesdisc.eosdis.nasa.gov:443/opendap/MERRA2/M2T1NXLND.5.12.4/'
    RADfile_start = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXLND.5.12.4/'
    RADfile_end   = '/MERRA2_400.tavg1_2d_lnd_Nx.'+YYYY+MM+DD+'.nc4'

    filenameRAD   = RADfile_start+YYYY+'/'+MM+RADfile_end;
    
    return(filenameTEMP, filenameRAD)

def specic_humidty_to_relative_humidity(Specific_H, T_air_K, Psurf_PA):
    Psurf_mb = 0.01*Psurf_PA;
    T_air_C = T_air_K-273.15;
    Esat_mb =  6.112 * np.exp((17.67 * T_air_C)/(T_air_C + 243.5)); 
    # source: https://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
    E_mb = Psurf_mb/(0.622/Specific_H + 0.378);
    RH = E_mb/Esat_mb;
    return RH

def specic_humidty_to_vapor_pressure(Specific_H, T_air_K, Psurf_PA):
    # source: https://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
    Psurf_mb = 0.01*Psurf_PA;
    T_air_C = T_air_K-273.15;
    E_mb = Psurf_mb/(0.622/Specific_H + 0.378);
    E_PA = E_mb/0.01
    return E_PA

def daily_mean(VAR_24hrs):
    VAR_day_mean = np.nanmean(VAR_24hrs,0)
    return VAR_day_mean

def daily_max(VAR_24hrs):
    VAR_day_max = np.nanmax(VAR_24hrs,0)
    return VAR_day_max

def check_directory(destination_directory):
    if not os.path.exists(destination_directory):
        print('creating directory: '+destination_directory)
        os.makedirs(destination_directory)
    else:
        print('directory exists: '+ destination_directory)
    return None

# DOES NOT WORK WITH NAN VALUES OVER THE OCEAN
# from scipy.interpolate import interp2d
#def regrid_MERRA_to_EASE_2_grid(data,data_lon,data_lat,EASE2_lon,EASE2_lat):
#    ip = interp2d(data_lon, data_lat, data)
#    EASE_2_data = ip(EASE2_lon,EASE2_lat)
#    return EASE_2_data

def regrid_MERRA_to_EASE_2_grid(data,data_lon,data_lat,EASE2_lon,EASE2_lat):
    x, y = np.meshgrid(data_lon,data_lat);
    x = x + 0.625/2.
    z = data;
    x=x.ravel()
    x=list(x[x!=np.isnan])
    y=y.ravel()
    y=list(y[y!=np.isnan])
    z=z.ravel()
    z=list(z[z!=np.isnan])
    EASE_2_data = griddata((x, y), z, (EASE2_lon[None,:], EASE2_lat[:,None]), method='linear')
    return EASE_2_data

def EASE2_36km_LON_LAT():
    path_to_lon = path_root+'PTJPL_model_Purdy_et_al_2017/data/LON_LAT/SMAP_L3P_LON_1d_36km_global.csv'
    path_to_lat = path_root+'PTJPL_model_Purdy_et_al_2017/data/LON_LAT/SMAP_L3P_LAT_1d_36km_global.csv'
    LON = np.genfromtxt(path_to_lon, dtype=float, delimiter=',') 
    LAT = np.genfromtxt(path_to_lat, dtype=float, delimiter=',') 
    return LON, LAT

def EASE2_9km_LON_LAT():
    path_to_lon = path_root+'PTJPL_model_Purdy_et_al_2017/data/LON_LAT/SMAP_L4_LON_1d_global.csv'
    path_to_lat = path_root+'PTJPL_model_Purdy_et_al_2017/data/LON_LAT/SMAP_L4_LAT_1d_global.csv'
    
    LON = np.genfromtxt(path_to_lon, dtype=float, delimiter=',')
    LAT = np.genfromtxt(path_to_lat, dtype=float, delimiter=',')
    return LON, LAT


days_num = np.arange(23,32);
days = list(map(str, days_num))

for y in np.arange(0,np.size(years)):
    for m in np.arange(0,np.size(months)):
        for d in np.arange(0,np.size(days)):
            print(years[y], months[m], days[d])
            filenameTEMP, filenameRAD = create_MERRA_filenames(years[y], months[m], days[d].zfill(2))
            print(filenameTEMP)
            fileNAMETEMP = filenameTEMP.split('/')[8];
            print(filenameRAD)
            fileNAMERAD = filenameRAD.split('/')[8];

            try:
                print('opening MERRA2 data for: '+years[y]+'-'+months[m]+'-'+days[d].zfill(2))
                os.system('wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies '+filenameTEMP)
                os.system('wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies '+filenameRAD)
                datasetT = xr.open_dataset(fileNAMETEMP);
                print('temperature and vapor data download complete')
                datasetR = xr.open_dataset(fileNAMERAD);
                print('radiation data download complete')
                os.system('rm ' +fileNAMETEMP)
                print('temperature and vapor data deleted from local system')
                os.system('rm ' +fileNAMERAD)
                print('radiation data deleted from local system')
                data_avail=True
            except:
                print('No available data for '+years[y]+'-'+months[m]+'-'+days[d].zfill(2))
                data_avail=False
            
            if data_avail:
                # create variables
                t2m_daily_mean  = daily_mean(datasetT['T2M'].data[:,:,:]);
                t2m_daily_max   = daily_max(datasetT['T2M'].data[:,:,:])
                ea2m_daily_mean = specic_humidty_to_vapor_pressure(daily_mean(datasetT['QV2M'].data[:,:,:]), t2m_daily_mean, daily_mean(datasetT['PS'].data[:,:,:]))
                rnet_daily_mean = daily_mean(datasetR['SWLAND'].data[:,:,:]) + daily_mean(datasetR['LWLAND'].data[:,:,:]);
                rnet_daily_max  = daily_max(datasetR['SWLAND'].data[:,:,:])  + daily_mean(datasetR['LWLAND'].data[:,:,:]); print(np.nanmax(np.nanmax(rnet_daily_max)))

                LAT=datasetT['lat'].data[:].astype(float)
                LON=datasetT['lon'].data[:].astype(float)
                # resample to 36km EASE2Grid
                for res in [36]:# 9,36
                    if res == 9:
                        ease2lon,ease2lat = EASE2_9km_LON_LAT();
                        out_dir = path_root+'PTJPL_model_Purdy_et_al_2017/data/data_global/merra_forcing/ease2_'+str(res)+'km/'
                    if res == 36:
                        ease2lon,ease2lat=EASE2_36km_LON_LAT()
                        out_dir = path_root+'PTJPL_model_Purdy_et_al_2017/data/data_global/merra_forcing/ease2_'+str(res)+'km/'

                    # resample to EASE2Grid
                    t2_daily_mean_36km   = regrid_MERRA_to_EASE_2_grid(t2m_daily_mean,LON,LAT,ease2lon,ease2lat);


                    t2_daily_max_36km    = regrid_MERRA_to_EASE_2_grid(t2m_daily_max,LON,LAT,ease2lon,ease2lat);
                    ea2m_daily_mean_36km = regrid_MERRA_to_EASE_2_grid(ea2m_daily_mean,LON,LAT,ease2lon,ease2lat);
                    rnet_daily_mean_36km = regrid_MERRA_to_EASE_2_grid(rnet_daily_mean,LON,LAT,ease2lon,ease2lat);
                    rnet_daily_max_36km  = regrid_MERRA_to_EASE_2_grid(rnet_daily_max,LON,LAT,ease2lon,ease2lat);
                    # Masks out ocean areas and assings data to null values
                    rnet_daily_mean_36km[rnet_daily_mean_36km>3000]      = np.nan;
                    rnet_daily_max_36km[rnet_daily_max_36km>3000]        = np.nan;
                    t2_daily_mean_36km[np.isnan(rnet_daily_mean_36km)]   = np.nan;
                    t2_daily_max_36km[np.isnan(rnet_daily_mean_36km)]    = np.nan;
                    ea2m_daily_mean_36km[np.isnan(rnet_daily_mean_36km)] = np.nan;
                    
                    # Masks out ocean areas and assigns null values
                    ds = xr.Dataset({'temperature_mean': (['lat', 'lon'], t2_daily_mean_36km),
                                     'temperature_max':  (['lat', 'lon'], t2_daily_max_36km),
                                     'ea2m_daily_mean':  (['lat', 'lon'], ea2m_daily_mean_36km),
                                     'rnet_daily_mean':  (['lat', 'lon'], rnet_daily_mean_36km),
                                     'rnet_daily_max':   (['lat', 'lon'], rnet_daily_max_36km)},
                                    coords={'lat': ease2lat,
                                            'lon': ease2lon})
                                            
                    # save dataset to netcdf file
                    try:
                        ds.to_netcdf(out_dir+'merra2_forcing_'+str(res)+'km_'+years[y]+months[m]+days[d].zfill(2)+'.nc')
                        print('MERRA2 resampled and saved for: '+years[y]+'-'+months[m]+'-'+days[d].zfill(2))
                        print('file saved to -> '+out_dir+'merra2_forcing_'+str(res)+'km_'+years[y]+months[m]+days[d]+'.nc')
                        print('\n')
                    except:
                        print('file already created \n\n')
                        pass
            
print('All MERRA2 forcing data resampled and saved to: '+out_dir)
