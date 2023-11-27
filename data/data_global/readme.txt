PT-JPL SMAP ET Data Description:

The PT-JPL model with SMAP data requires Inputs of:
	air_temperature_K,
	air_temperature_mean_K,
	water_vapor_pressure_mean_Pa,
	ndvi_mean,
	optimum_temperature,
	fAPARmax,
	NET RADIATION,
	WP, # < new for SMAP ET
	FC, # < new for SMAP ET
	CH, # < new for SMAP ET
	SM, # < new for SMAP ET):


DATA DESCRIPTION:

NET RADIATION
Script to regrid to EASE2grid: 
	scripts/merra2_to_ease2grid.py
Folders:
	merra_forcing/ease2_36km/
	merra_forcing/ease2_9km/
Filename:
	merra2_forcing_36km_YYYYMMDD.nc
Variable Name:
	rnet_daily_mean
	Size:       964x406
	Size:       3856x1624
	Units: W/m2

AIR TEMPERATURE
Script to regrid to EASE2grid: 
	scripts/merra2_to_ease2grid.py
Folders:
	merra_forcing/ease2_36km/
	merra_forcing/ease2_9km/
Filename:
	merra2_forcing_36km_YYYYMMDD.nc
Variable Name:
	temperature_max
	Size:       964x406
	Size:       3856x1624
	Units: K

AIR TEMPERATURE MEAN 
Script to regrid to EASE2grid: 
	scripts/merra2_to_ease2grid.py
Folders:
	merra_forcing/ease2_36km/
	merra_forcing/ease2_9km/
Filename:
	merra2_forcing_36km_YYYYMMDD.nc
Variable Name:
	temperature_mean
	Size:       964x406
	Size:       3856x1624
	Units: K

Vapor Pressure
Script to regrid to EASE2grid: 
	/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/scripts/merra2_to_ease2grid.py
Folders:
	merra_forcing/ease2_36km/
	merra_forcing/ease2_9km/
Filename:
	merra2_forcing_36km_YYYYMMDD.nc
Variable Name:
	ea2m_daily_mean
	Size:       964x406
	Size:       3856x1624
	Units: Pa

NDVI
Script to regrid to EASE2grid: 
	data/RAW/MOD13C1_to_EASE_grid.m
Folders:
	modis_ndvi/ease2_36km/
	modis_ndvi/ease2_9km/
Filename:
	MOD13C1.AYYYYDOY_36km_grid.nc
	MOD13C1.AYYYYDOY_9km_grid.nc
Variable Name:
	NDVI
	Units: NA 0.0 - 1.0

Canopy Height
Script to regrid to EASE2grid: 
	data/RAW/rescale_ch_to_ease2grid.m
Folders:
	canopy_height/
Filename:
	canopy_height_36km_grid.nc
	canopy_height_9km_grid.nc
Variable Name:
	CH
	Units: m
	description: mean canopy height for each ease2 grid cell

fAPAR maximum value
Script to regrid to coarsen to 36km: 
	scripts/find_faparmax.m
Folders:
	fapar_max/
Filename:
	fAPAR_max_36km_grid.nc
	fAPAR_max_9km_grid.nc
Variable Name:
	fAPAR_max
	Units: % [data are valid from 0-1]

Soil Properties
Script to obtain Soil Properties data
	data_global/_fix_res_scripts_/fix_res_soil_props.m
Links to data description:
	SMAP L4 Data Soil Propeties ATBD by Narandera Das
Folders:
	soil_props/	
Filename:
	porosity_9km_grid.nc
	porosity_36km_grid.nc
	wilting_point_9km_grid.nc
	wilting_point_36km_grid.nc
Variable Name:
	porosity
	wp
	Units: % [data are valid from 0-1]


Soil Moisture
Script to obtain SMAP data
	/Volumes/AJ_RESEARCH/SMAP_ET_2_global_application/scripts/get_SMAP.py
Links to data description:
	SMAP Passive L3 Data:	https://nsidc.org/data/SPL3SMP/versions/4
	SMAP Passive L4 Data:	https://nsidc.org/data/SPL3SMP_E/versions/1
Folders:
	NA data obtained through opendaap
Filename:
	SMAP L3 36km Passive
	http://n5eil01u.ecs.nsidc.org:80/opendap/SMAP/SPL3SMP.003/2015.04.06/SMAP_L3_SM_P_20150406_R13080_001.h5
	http://n5eil01u.ecs.nsidc.org:80/opendap/SMAP/SPL3SMP.004/2017.01.21/SMAP_L3_SM_P_20170121_R14010_001.h5
	Enhanced SMAP-SENTINEL 1 data
	http://n5eil01u.ecs.nsidc.org:80/opendap/SMAP/SPL3SMP_E.001/2016.03.02/SMAP_L3_SM_P_E_20160302_R14010_001.h5
Variable Name:
	Soil_Moisture_Retrieval_Data_AM_soil_moisture
	Soil_Moisture_Retrieval_Data_PM_soil_moisture_pm
#	Size:       964x406   #<-- NEED TO DOUBLE CHECK THESE
#	Size:       3856x1624 #<-- NEED TO DOUBLE CHECK THESE
	Units: K

	Citation: 
	O'Neill, P. E., S. Chan, E. G. Njoku, T. Jackson, and R. Bindlish. 2016. SMAP L3 Radiometer Global Daily 36 km EASE-Grid Soil Moisture, Version 4. [Indicate subset used]. Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center. doi: http://dx.doi.org/10.5067/OBBHQ5W22HME. [Date Accessed].
