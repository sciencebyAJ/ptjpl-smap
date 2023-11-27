%% get soil properties to 9 and 36km resolutions:
% only put wilting point into a netcdf
% only put porosity into a netcdf file
clear 
clc

%% LIST FILES IN THE CURRENT DIRECTORY
rel_path ='../soil_props/';
runsDir = dir(strcat(rel_path,'S*.h5')); %# Get the data for the current directory
dirIndex = [runsDir.isdir];  %# Find the index for directories
fileList = {runsDir(~dirIndex).name}';
%% Display meta data
h5disp(strcat(rel_path,fileList{1}))

%% Names of data for L3 AP dataset
% Define the HDF5 Group where to access the data
grp_id = '/Land-Model-Constants_Data/';
wp_id  = 'clsm_wp';
por_id = 'clsm_poros';

%% Quick function to grab data from group based on dataset name:
% filenum is the file name; 
% file_id is the dataset name
% rel_path is path to the file; 
% grp_id is the hdf5 group; 
getdata  = @(file_num, file_id) h5read(strcat(rel_path,file_num), strcat(grp_id, file_id));

%% Grab porosity and wilting point datasets
file1  = fileList{1};
cell_time_i = strsplit(file1,'_');
ymd_i = cell_time_i{5};
% Extract data and filter for nan values;
wilt_point = getdata(file1,wp_id);   wilt_point(wilt_point==-9999)=nan;
porosity   = getdata(file1,por_id);  porosity(porosity==-9999)=nan;

%% Create netcdf files to use
path_out = '../soil_props/';
outfile_filename_por = strcat(path_out,'porosity_9km_grid.nc');
nccreate(outfile_filename_por,'porosity','Dimensions',{'lon',size(porosity,1),'lat',size(porosity,2)},'Format','classic');
ncwrite(outfile_filename_por,'porosity',porosity);

outfile_filename_por = strcat(path_out,'wilting_point_9km_grid.nc');
nccreate(outfile_filename_por,'wp','Dimensions',{'lon',size(wilt_point,1),'lat',size(wilt_point,2)},'Format','classic');
ncwrite(outfile_filename_por,'wp',wilt_point);

%% Load geospatial information to regrid datasets
% read in latitude and longitude data
filename_LON9 = '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/LON_LAT/SMAP_L4_LON_1d_global.csv';
filename_LAT9 = '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/LON_LAT/SMAP_L4_LAT_1d_global.csv';
filename_LON36 = '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/LON_LAT/SMAP_L3P_LON_1d_36km_global.csv';
filename_LAT36 = '/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/LON_LAT/SMAP_L3P_LAT_1d_36km_global.csv';

lon_1d9 = csvread(filename_LON9);
lat_1d9 = csvread(filename_LAT9);
[LON_9,LAT_9]=meshgrid(lon_1d9,lat_1d9);
clear lon_1d9 lat_1d9 filename_LAT9 filename_LON9
lon_1d36 = csvread(filename_LON36);
lat_1d36 = csvread(filename_LAT36);
[LON_36,LAT_36]=meshgrid(lon_1d36,lat_1d36);
clear lon_1d36 lat_1d36 filename_LAT36 filename_LON36

%% Check geospatial information
figure(); pcolor(LON_9,LAT_9,porosity');   shading flat   % note need to transpose porosity data
figure(); pcolor(LON_9,LAT_9,wilt_point'); shading flat % note need to transpose wilting point data

%% Regrid data to EASE2 grid
% disp('Regridding data to EASE2 grid');
porT = double(porosity');
wpT  = double(wilt_point');

porosity_res36     = griddata(LON_9,LAT_9, porT, LON_36, LAT_36);

wiltingpoint_res36 = griddata(LON_9,LAT_9, wpT,  LON_36, LAT_36);
%%
por = porosity_res36';
wp = wiltingpoint_res36';
size(por)
%% Check geospatial information
figure(); pcolor(LON_36,LAT_36,por');   shading flat   % note need to transpose porosity data
figure(); pcolor(LON_36,LAT_36,wp' ); shading flat % note need to transpose wilting point data

%% Create netcdf files to use
path_out = '../soil_props/';
outfile_filename_por = strcat(path_out,'porosity_36km_grid.nc');
nccreate(outfile_filename_por,'porosity','Dimensions',{'lon',size(por,1),'lat',size(por,2)},'Format','classic');
ncwrite(outfile_filename_por,'porosity',por);

outfile_filename_por = strcat(path_out,'wilting_point_36km_grid.nc');
nccreate(outfile_filename_por,'wp','Dimensions',{'lon',size(wp,1),'lat',size(wp,2)},'Format','classic');
ncwrite(outfile_filename_por,'wp',wp);


