import numpy as np
import os
import xarray as xr
import sys
from pydap.client import open_url

res         = int(sys.argv[1])
year_in     = int(sys.argv[2])
month_start = int(sys.argv[3])
month_end   = int(sys.argv[4])

print '\n\nNow downloading data from SMAP\t'+str(year_in)+'-'+str(month_start)+'\t'+str(year_in)+'-'+str(month_end)


if res ==9:
    data_out = '/Volumes/AJ_4TB/Coding/Research/PTJPL_model_Purdy_et_al_2017/data/data_global/smap_sm/ease2_9km/'

if res ==36:
    data_out = '/Volumes/AJ_4TB/Coding/Research/PTJPL_model_Purdy_et_al_2017/data/data_global/smap_sm/ease2_36km/'

def check_directory(destination_directory):
    if not os.path.exists(destination_directory):
        print("output directory does not exist \nmaking new directory \ndata saved to:\t"+ destination_directory)
        os.makedirs(destination_directory)
    else:
        print("output directory exists \ndata saved to:\t"+ destination_directory)
    return None

def SMAP_L3_P_36km_Path(yyyy,mm,dd):#,PRODUCT,VERSION):
    fpath_start = 'https://n5eil01u.ecs.nsidc.org/SMAP/'
    product = 'SPL3SMP';
    version = '.004/'
    year = str(yyyy); month = str(mm).zfill(2); day = str(dd).zfill(2);
    yearmonthday= year+month+day
    product_long = '/SMAP_L3_SM_P_'
    fpath_end = '_R14010_001.h5'
    SMAP_DATA_PATH = fpath_start +product+version+year+'.'+month+'.'+day+product_long+yearmonthday+fpath_end
    fileName = SMAP_DATA_PATH.split('/')[6]
    return SMAP_DATA_PATH,fileName

def SMAP_L3_P_E_9km_Path(yyyy,mm,dd):
    #SMAP SM Passive Enhanced
    smap_path_start = 'https://n5eil01u.ecs.nsidc.org/SMAP/SPL3SMP_E.001/'
    year = str(yyyy)+'.'
    month = str(mm).zfill(2)+'.'
    day = str(dd).zfill(2)+'/'
    smap_product_ln = 'SMAP_L3_SM_P_E'+'_'
    yyyymmdd = str(yyyy)+str(mm).zfill(2)+str(dd).zfill(2)
    smap_path_end = '_R14010_001.h5'
    SMAP_DATA_PATH = smap_path_start+year+month+day+smap_product_ln+yyyymmdd+smap_path_end
    fileName = SMAP_DATA_PATH.split('/')[6]
    return SMAP_DATA_PATH, fileName

def fname_out(fileName):
    ''' returns filename to save SMAP subset
        '''
    len(fileName.split('_'))
    fname_out=''
    for l in np.arange(0,len(fileName.split('_'))-2):
        fname_out = fname_out+(fileName.split('_')[l])+'_'
    fname_out= fname_out[:-1]+'.nc'
    return fname_out

def subset_SMAP(fileName):
    f = h5py.File(fileName,'r')
    sm_dataAM     =  f['Soil_Moisture_Retrieval_Data_AM_soil_moisture'][:,:]
    qc_flag_am    =  f['Soil_Moisture_Retrieval_Data_AM_retrieval_qual_flag'][:,:]
    #     sm_dataAM[((qc_flag_am>>0)&1)==1]=np.nan;
    #     sm_dataAM[(sm_dataAM==-9999.0)]=np.nan;
    
    sm_dataPM  =  f['Soil_Moisture_Retrieval_Data_PM_soil_moisture_pm'][:,:]
    qc_flag_pm =  f['Soil_Moisture_Retrieval_Data_PM_retrieval_qual_flag_pm'][:,:]
    #     sm_dataPM[((qc_flag_pm>>0)&1)==1]=np.nan;
    #     sm_dataPM[sm_dataPM==-9999.0]=np.nan;
    sm_both=np.ones((np.shape(sm_dataAM)[0],np.shape(sm_dataAM)[1],2))
    sm_both[:,:,0]=sm_dataAM;
    sm_both[:,:,1]=sm_dataPM;
    sm_avg = np.nanmax(sm_both,axis=2)
    pathOUT = data_out
    sm_ds = xr.Dataset({'soil_moisture':(['x','y'],sm_avg),
                       'qc_flag_am':(['x','y'],qc_flag_am),
                       'qc_flag_pm':(['x','y'],qc_flag_pm)})
    fNameOUT = fname_out(fileName)
    sm_ds.to_netcdf(pathOUT +fNameOUT)
    print 'smap data saved to '+pathOUT+fNameOUT


WGET_COMMAND_START = "wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies "

home_dir = os.getcwd();
check_directory(data_out)

yyyy=year_in;
MM = np.arange(month_start,month_end);
DD = np.arange(1,2);
for mm in MM[:]:
    for dd in DD[:]:
        try:
            print 'opening SMAP data for: '+str(yyyy)+'-'+str(mm).zfill(2)+'-'+str(dd).zfill(2)
            if res == 36:
                full_path36,fileName36 = SMAP_L3_P_36km_Path(yyyy,mm,dd)
                WGET_COMMAND = WGET_COMMAND_START + full_path36
                print WGET_COMMAND
                os.system(WGET_COMMAND)
                subset_SMAP(fileName36)
                os.system('rm '+fileName36)
                print fileName36+' downloaded & subset'
                data_avail=True
            if res ==9:
                full_path9,fileName9 = SMAP_L3_P_E_9km_Path(yyyy,mm,dd)
                WGET_COMMAND = WGET_COMMAND_START + full_path9
                print full_path9
                    
                os.system(WGET_COMMAND)
                subset_SMAP(fileName9)
                os.system('rm '+fileName9)
                print full_path9
                print fileName9+' downloaded & subset'
                data_avail=True
        except:
            data_avail=False
            print 'No available data for '+str(yyyy)+'-'+str(mm).zfill(2)+'-'+str(dd).zfill(2)
            pass
        if data_avail:
            if res ==36:
                print '*** SM data saved from \t'+fileName36+'\t *** '
                print fileName36 + ' removed from scratch'
            if res ==9:
                print '*** SM data saved from \t'+fileName9+'\t *** '
                print fileName9 + ' removed from scratch'


print '\n'
print 'SCRIPT COMPLETE'