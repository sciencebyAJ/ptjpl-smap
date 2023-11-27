% find FAPAR MAX and save file
% find TOPT and save file
current = dir('/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/data_conus/modis_ndvi/ndvi*.nc');
% Use current directory:
current_list = {current.name}';
fAPAR_max = NaN(size(ncread(current_list{1},'ndvi')));
NDVI_all = NaN(size(ncread(current_list{1},'ndvi'),1),size(ncread(current_list{1},'ndvi'),2),size(current_list,1));
NDVI_yr = NaN(size(ncread(current_list{1},'ndvi'),1),size(ncread(current_list{1},'ndvi'),2),365);
for i = 1:size(current_list,1)
    filename_in  = current_list{i};
    NDVI_all(:,:,i) = ncread(filename_in,'ndvi');
    doy = str2num(filename_in(10:12));
    NDVI_yr(:,:,doy) =ncread(filename_in,'ndvi');
end
SAVI_MULT = 0.45;
SAVI_ADD = 0.132;
savi_all = NDVI_all.* SAVI_MULT + SAVI_ADD;
FAPAR_MULT = 1.3632;
FAPAR_ADD = -0.048;
fAPAR_all = savi_all.* FAPAR_MULT + FAPAR_ADD;
fAPAR_max = max(fAPAR_all,[],3);
clear savi_all fAPAR_all NDVI_all

%% radiation all
cd '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/data_conus/daymet_rnet/'
current = dir('/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/data_conus/daymet_rnet/*.nc');
current_list = {current.name}';
RNET_yr = NaN(size(ncread(current_list{1},'srad'),1),size(ncread(current_list{1},'srad'),2),365);
for i = 1:size(current_list,1)
    filename_in  = current_list{i};
    doy = str2num(filename_in(21:23));
    RNET_yr(:,:,doy) =ncread(filename_in,'srad');
end

%% temperature all
cd '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/data_conus/temp/'
current = dir('/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/data_conus/temp/*tmax*.nc');
current_list = {current.name}';
TMAX_yr = NaN(size(ncread(current_list{1},'tmax'),1),size(ncread(current_list{1},'tmax'),2),365);
for i = 1:size(current_list,1)
    filename_in  = current_list{i};
    doy = str2num(filename_in(21:23));
    TMAX_yr(:,:,doy) =ncread(filename_in,'tmax');
end

%% vapor pressure all
cd '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/data_conus/temp/'
current = dir('/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/data_conus/temp/*vp*.nc');
current_list = {current.name}';
VP_yr = NaN(size(ncread(current_list{1},'vp'),1),size(ncread(current_list{1},'vp'),2),365);
for i = 1:size(current_list,1)
    filename_in  = current_list{i};
    doy = str2num(filename_in(19:21));
    VP_yr(:,:,doy) =ncread(filename_in,'vp');
end

%% saturated vapor pressure 
SVP_BASE = 0.611;
SVP_MULT = 17.27;
SVP_ADD = 237.7;
SVP_yr = SVP_BASE*exp((TMAX_yr.* SVP_MULT)./(TMAX_yr + SVP_ADD));

%% Relative Humidity
RH_yr = VP_yr*0.001./SVP_yr;
RH_yr(RH_yr>1.0) = NaN;
%% Vapor Pressure Deficit

VPD_yr = SVP_yr-RH_yr.*SVP_yr;
%% FPAR yr
FPAR_yr = (NDVI_yr.* SAVI_MULT + SAVI_ADD).* FAPAR_MULT + FAPAR_ADD;

%% Find Maximum index for each value:
PHEN_ind = NaN(size(fAPAR_max,1),size(fAPAR_max,2),336);
TMAX_vals = NaN(size(fAPAR_max,1),size(fAPAR_max,2),336);
for i = 15:350
    FPAR_use = nanmean(FPAR_yr(:,:,i-14:i+15),3);
    RNET_use = nanmean(RNET_yr(:,:,i-14:i+15),3);
    TMAX_use = nanmean(TMAX_yr(:,:,i-14:i+15),3);
    VPD_use  = nanmean(VPD_yr(:,:,i-14:i+15), 3);
    TMAX_vals(:,:,i) = TMAX_use;
    PHEN_ind(:,:,i) = FPAR_use.*RNET_use.*TMAX_use./VPD_use;
end
%%
Topt = NaN(size(fAPAR_max,1),size(fAPAR_max,2));
for j = 1:size(fAPAR_max,1)
    for k = 1:size(fAPAR_max,2)
        [Val, loc] = max(PHEN_ind(j,k,:));
        Topt(j,k)  = TMAX_vals(j,k,loc);    
    end
end
Topt(Topt<0)=NaN;

%% Write OUTPUT to Files
% outfilefAPARfilename = '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/data_conus/fapar_max/faparMAX.nc';
% nccreate(outfilefAPARfilename,'fapar_max','Dimensions',{'lat',length(lat_conus),'lon',length(lon_conus)},'Format','classic');
% ncwrite(outfilefAPARfilename,'fapar_max',fAPAR_max);

outfilefAPARfilename = '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/data_conus/t_opt/topt.nc';
nccreate(outfilefAPARfilename,'topt','Dimensions',{'lat',size(Topt,1),'lon',size(Topt,2)},'Format','classic');
ncwrite(outfilefAPARfilename,'topt',Topt);

%% Topt
filename_LON = '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/SMAP_L4_LON_1d_CONUS.csv';
filename_LAT = '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/SMAP_L4_LAT_1d_conus.csv';
lon_1d = csvread(filename_LON);
lat_1d = csvread(filename_LAT);

[LONconus,LATconus]=meshgrid(lon_1d,lat_1d);
%%
figure(); pcolor(LONconus,LATconus,RNET_use);shading flat