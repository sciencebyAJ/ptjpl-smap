%% REGRID TEST FOR DAYMET:
clear
clc
filename_LON      = '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/SMAP_L4_LON_1d_CONUS.csv';
filename_LAT      = '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/SMAP_L4_LAT_1d_CONUS.csv';
filename_NDVI     = '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/modis_download/ndvi_2015_wgs84/ndvi_2015017.nc';
filename_NDVI_LON = '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/modis_download/ndvi_2015_wgs84/MODIS_LON_1d_CONUS.csv';
filename_NDVI_LAT = '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/modis_download/ndvi_2015_wgs84/MODIS_LAT_1d_CONUS.csv';
% filename_DAYMET = '/Volumes/AJ_RESEARCH/RAW_DATA/DAYMET/daymet_v3_vp_2015_066.nc';
% lonDM = ncread(filename_DAYMET,'lon'); lonDM1d = reshape(lonDM,[],1);   clear lonDM
% latDM = ncread(filename_DAYMET,'lat'); latDM1d = reshape(latDM,[],1);   clear latDM
% parts = strsplit(filename_DAYMET,'_'); varname = parts{5};
% dataDM = ncread(filename_DAYMET,varname); dataDM1d = reshape(dataDM,[],1); clear dataDM

%% LOAD NDVI LAT LON OF DATASET for CONUS
lonNDVI = csvread(filename_NDVI_LON);     
latNDVI = csvread(filename_NDVI_LAT);     
[lonNDVI2d,latNDVI2d] = meshgrid(lonNDVI,latNDVI);
% lonNDVI1d  = reshape(lonNDVI2d,[],1);  %clear lonNDVI2d lonNDVI
% latNDVI1d  = reshape(latNDVI2d,[],1);  %clear latNDVI2d latNDVI
%% LOAD SMAP 9km LAT LON for CONUS
lon_1d = csvread(filename_LON);
lat_1d = csvread(filename_LAT);
[lonSMAP,latSMAP] = meshgrid(lon_1d,lat_1d);clear latSMAP
%% For loop to put around below
ndvi_scale_factor = 0.0001;
lon_dif = abs(lon_1d(1)-lon_1d(2))/2.0;

fileStart = '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/modis_download/ndvi_2015_wgs84/';
fileList = dir('/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/modis_download/ndvi_2015_wgs84/*.nc');
for k = 1:size(fileList,1)
    filename_NDVI = strcat(fileStart,fileList(k).name)
    dataNDVI = ncread(filename_NDVI, 'ndvi'); 
    dataNDVI2d = dataNDVI*ndvi_scale_factor;clear dataNDVI
    dataNDVI = dataNDVI2d; clear dataNDVI2d
    data_out = NaN(size(lonSMAP)); 
    for i = 1:size(lat_1d,1)-1
        if i ==1
            lat_N = lat_1d(i)+abs((lat_1d(i+1)-lat_1d(i))/2.0);
            lat_S = lat_1d(i)-abs((lat_1d(i+1)-lat_1d(i))/2.0);
        elseif i == size(lat_1d)
            lat_N = lat_1d(i)+abs((lat_1d(i-1)-lat_1d(i))/2.0);
            lat_S = lat_1d(i)-abs((lat_1d(i-1)-lat_1d(i))/2.0);
        else
            lat_N = lat_1d(i)+abs((lat_1d(i-1)-lat_1d(i))/2.0);
            lat_S = lat_1d(i)-abs((lat_1d(i+1)-lat_1d(i))/2.0);
        end
        i
        lat_ind=find(latNDVI2d<lat_N & latNDVI2d>lat_S);
        parfor j = 1:size(lon_1d,1)-1
            lon_W = lon_1d(j)-lon_dif;
            lon_E = lon_1d(j)+lon_dif;
            lon_ind = find(lonNDVI2d>lon_W & lonNDVI2d<lon_E);
            lat_lon_ind = intersect(lat_ind,lon_ind);
            if size(lat_lon_ind,1)==0
                'not yet'
            else
                data_out(i,j) = nanmean(dataNDVI(lat_lon_ind));
            end
        end
    end
    tempname = strsplit(fileList(k).name,'.');
    nccreate(strcat(tempname{1},'_smap_l4_res.nc'),'ndvi','Dimensions',{'lat',length(lat_1d),'lon',length(lon_1d)},'Format','classic')
    ncwrite(strcat(tempname{1},'_smap_l4_res.nc'),'ndvi',data_out)
end
%%
% [lonSMAP,latSMAP] = meshgrid(lon_1d,lat_1d);
% figure(); pcolor(lonSMAP,latSMAP,data_out);shading flat
%
% figure(); pcolor(lonNDVI2d,latNDVI2d,dataNDVI2d');shading flat

%%
% t = cputime;
% 
% for i = 1:size(lat_1d,1)-1
%     if i ==1
%         lat_N = lat_1d(i)+abs((lat_1d(i+1)-lat_1d(i))/2.0);
%         lat_S = lat_1d(i)-abs((lat_1d(i+1)-lat_1d(i))/2.0);
%     elseif i == size(lat_1d)
%         lat_N = lat_1d(i)+abs((lat_1d(i-1)-lat_1d(i))/2.0);
%         lat_S = lat_1d(i)-abs((lat_1d(i-1)-lat_1d(i))/2.0);
%     else
%         lat_N = lat_1d(i)+abs((lat_1d(i-1)-lat_1d(i))/2.0);
%         lat_S = lat_1d(i)-abs((lat_1d(i+1)-lat_1d(i))/2.0);
%     end
%     i
%     lat_ind=find(latNDVI1d<lat_N & latNDVI1d>lat_S);
%     parfor j = 1:size(lon_1d,1)-1
%         lon_W = lon_1d(j)-lon_dif;
%         lon_E = lon_1d(j)+lon_dif;
%         lon_ind = find(lonNDVI1d>lon_W & lonNDVI1d<lon_E);
%         lat_lon_ind = intersect(lat_ind,lon_ind);
%         if size(lat_lon_ind,1)==0
%         else
%             data_out(i,j) = nanmean(dataNDVI1d(lat_lon_ind));
%         end
%         
%     end
% end
% e = cputime-t
