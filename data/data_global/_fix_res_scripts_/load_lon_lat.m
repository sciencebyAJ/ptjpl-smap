function [l on1d, lat1d] = load_lon_lat()

filePath = '/data15/famiglietti2/ajpurdy/NASA_NESSF_dir/data/LON_LAT/';
lonName = 'SMAP_L3P_LAT_1d_36km_global.csv';
latName = 'SMAP_L3P_LON_1d_36km_global.csv';

lon1d = csvread(strcat(filePath,lonName));
lat1d = csvread(strcat(filePath,latName));
end
