# Module that contains the model for a FLUXNET site and a pandas dataframe 2-d array for evapotranspiration used in smap analysis
#

import numpy as np
from numpy import exp, log
import xarray as xr

DEFAULT_AVERAGING_PERIOD = 30
DEFAULT_FLOOR_SATURATION_VAPOR_PRESSURE = True

# Priestley-Taylor coefficient alpha
PRIESTLEY_TAYLOR_ALPHA = 1.26

BETA = 1.0

# psychrometric constant gamma in pascals over kelvin
# same as value for ventilated (Asmann type) psychrometers, with an air movement of some 5 m/s
# http://www.fao.org/docrep/x0490e/x0490e07.htm
PSYCHROMETRIC_GAMMA = 0.0662 # Pa/K
KRN = 0.6
KPAR = 0.5

def Gfrac_calc(CH,LAI):
    '''
        inputs:  CH, LAI
        outputs: Gfrac
    '''
    # create empty array for Gfrac
    g_frac = np.ones(shape=np.shape(LAI));
    g_frac[:,:]=np.nan;
    # Paramters for tall vegetation cover
    A_tall = 0.087; B_tall = 0.15; C_tall = 0.88; D_tall = 0.05;
    g_frac[(CH>1.0)&(LAI>0.5)] = A_tall + B_tall * exp(-1.0*C_tall*LAI[(CH>1.0)&(LAI>0.5)]);
    g_frac[(CH>1.0)&(LAI<0.5)] = D_tall;
    # Paramters for short vegetation cover
    A_short = 0.019; B_short = 0.079; C_short = 0.44; D_short = 0.20;
    g_frac[(CH<1.0)&(LAI>0.5)] = A_short + B_short * exp(-1.0*C_short*LAI[(CH<1.0)&(LAI>0.5)]);
    g_frac[(CH<1.0)&(LAI<0.5)] = D_short;
    # Return Gfrac values
    return g_frac



# calculate Soil-Adjusted Vegetation Index from Normalized Difference Vegetation Index
# using linear relationship
def savi_from_ndvi(ndvi):
    SAVI_MULT = 0.45
    SAVI_ADD = 0.132
    savi = ndvi * SAVI_MULT + SAVI_ADD
    
    return savi


# calculate fraction of absorbed photosynthetically active radiation
# from Soil-Adjusted Vegetation Index using linear relationship
def fAPAR_from_savi(savi):
    FAPAR_MULT = 1.3632
    FAPAR_ADD = -0.048    
    fAPAR = savi * FAPAR_MULT + FAPAR_ADD
    return fAPAR

# calculate fraction of absorbed photosynthetically active radiation
# from Normalized Difference Vegetation Index
def fAPAR_from_ndvi(ndvi):
    savi = savi_from_ndvi(ndvi)
    fAPAR = fAPAR_from_savi(savi)
    return fAPAR

# calculate fraction of intercepted photosynthetically active radiation
# from Normalized Difference Vegetation Index
def fIPAR_from_ndvi(ndvi):
    FIPAR_ADD = -0.05
    fIPAR = ndvi + FIPAR_ADD    
    return fIPAR

# saturation vapor pressure in kPa from air temperature in celsius
def saturation_vapor_pressure_from_air_temperature(air_temperature):
    SVP_BASE = 0.611
    SVP_MULT = 17.27
    SVP_ADD = 237.7
    svp = SVP_BASE * exp((air_temperature * SVP_MULT) / (air_temperature + SVP_ADD))
    
    return svp


# calculate slope of saturation vapor pressure to air temperature
# in pascals over kelvin
def delta_from_air_temperature(air_temperature):
    return 4098 * (0.6108 * exp(17.27 * air_temperature / (237.7 + air_temperature))) / (air_temperature + 237.3) ** 2


def enforce_boundaries(matrix, lower_bound, upper_bound):
    matrix[matrix < lower_bound] = lower_bound
    matrix[matrix > upper_bound] = upper_bound
    
    return matrix
################################################################################################
############################## ALL CODE BELOW IS NEWLY MODIFIED ################################
################################################################################################
def fTREW_CR_fun(VWC, WP, FC, CH):
    '''RETURNS: Plant relative extractable water
        INPUTS: soil moisture in VWC, wilting point, field capactiy'''
    P1 = 1.0
    P2 = 2.0
    CHsqrt = np.sqrt(CH);
    CHsqrt[CHsqrt<1.0]=1.0;
    CHsqrt[CHsqrt>4.0]=4.0;
    WP = WP/CHsqrt;
    CR = WP + (FC-WP)*0.7;
    VWC = np.minimum(VWC, CR);
    return 1-((((CR)-VWC)/((CR)-(WP/P1)))**P2)

def LE_2_ETcm(LE_Wm2):
    '''
        This tool converts Latent Energy to Evapotranspiration
        INPUT DATA:  LE_2_ET (W/m2)
        OUTPUT DATA: ET_cm (cm/day)
        '''
    lambda_e = 2.460*10**6      # J kg^-1
    roe_w = 1000                # kg m^-3
    m_2_cm = 100                # convert m to mm
    s_2_30m = 60*30*48          # multiply by s in 30 min multiply by 48 to get cm/day
    mask = ~np.isnan(LE_Wm2)
    ET_cm = np.empty(LE_Wm2.shape)
    ET_cm[:] = np.NAN
    ET_cm[mask] = (LE_Wm2[mask]*(m_2_cm*s_2_30m)/(lambda_e*roe_w))
    return ET_cm


def fTREW_ch_pet_fun(theta_fc, theta_wp, pet, CH, theta_obs):
    '''
        INPUTS:
        theta_fc  = field capacity               [cm3/cm3]
        theta_wp  = wilting point                [cm3/cm3]
        pet       = potential evapotranspiration [cm/day]
        CH        = canopy height                [m]
        theta_obs = observed soil moisture       [cm3/cm3]
        OUTPUT:
        sens: plant sensitivity to soil moisture
        The above includes impacts of reducing the wilting point with increasing canopy height & reducing the point of inflection of soil moisture critical by canopy heigh AND potential evapotranspiration.
        
        '''
    alpha = 0.76# van diepen et al., 1988
    beta  = 1.5 # van diepen et al., 1988
    c     = 0.5 # parm to be optimized
    CHfactor = 5-c*CH;
    CHfactor[CHfactor<0.]=0.0;
    p   = (1.0/(alpha+beta*pet))-0.1*CHfactor; # van diepen et al., 1988
    theta_cr = (1-p)*(theta_fc-theta_wp)+theta_wp
    theta_cr = np.minimum(theta_cr,theta_fc)
    CH_max_scale = 4.0; # Cut off to avoid 'no-sensitivity to surface soil moisture'
    CH_min_scale = 2.0; # Cut off to ensure not double accounting of stress from original stress functions
    CH_scale = CH;
    CH_scale[CH_scale>CH_max_scale]=CH_max_scale;
    CH_scale[CH_scale<CH_min_scale]=CH_min_scale;
    CHsqrt = np.sqrt(CH);
    CHsqrt[CHsqrt<1.0]=1.0;
    CHsqrt[CHsqrt>4.0]=4.0;
    theta_wp = theta_wp/CHsqrt;
    theta_cr = np.maximum(theta_cr,theta_wp);
    theta_obs = np.maximum(theta_obs,theta_wp);
    theta_obs = np.minimum(theta_obs,theta_cr);
    sens = 1.0-((theta_cr-theta_obs)/(theta_cr-theta_wp))**(CH_scale);
    sens[sens<0.]=0.0
    sens[sens>1.0]=1.0
    return sens

def ptjpl_smap(lat,
               lon,
               air_temperature_K,
               air_temperature_mean_K,
               water_vapor_pressure_mean_Pa,
               air_temp_moving_window,
               water_vapor_moving_window,
               ndvi_mean,
               optimum_temperature,
               fAPARmax,
               net_radiation,
               WP, # < new for SMAP ET
               FC, # < new for SMAP ET
               CH, # < new for SMAP ET
               SM, # < new for SMAP ET
               verbose=True,
               floor_saturation_vapor_pressure=DEFAULT_FLOOR_SATURATION_VAPOR_PRESSURE):
    """
        :param air_temperature_K:
        Numpy matrix of air temperature near the surface in kelvin
        :param air_temperature_mean_K:
        Numpy matrix of average air temperature near the surface in kelvin
        :param water_vapor_pressure_mean_Pa:
        Numpy matrix of average vapor pressure in pascals
        :param air_temp_moving_window,
        Numpy matrix of 15 day average vapor air temperature in kelvin
        :param water_vapor_moving_window,
        Numpy matrix of 15 day average vapor pressure in pascals
        :param ndvi_mean:
        Numpy matrix of average Normalized Difference Vegetation Index
        :param optimum_temperature:
        Numpy matrix of phenologically optimum temperature in celsius
        :param fAPARmax:
        Numpy matrix of maximum fraction of photosynthetically active radiation
        :param net_radiation:
        Numpy matrix of instantaneous net radiation in watts per square meter
        :param WP, 
        Numpy matrix of Wiling Point
        :param FC,
        Numpy matrix of Field Capcity
        :param CH, 
        Numpy matrix of Canopy Height
        :param SM,
        Numpy matrix of SMAP soil moisture
        :param daily_radiation:
        Numpy matrix of daily net radiation in watts per square meter
        :param verbose:
        Flag to output activity to console
        :param floor_saturation_vapor_pressure:
        Option to floor calculation of saturation vapor pressure at 1 to avoid anomalous output
        :return:
        Named tuple of type PTJPLResults including the attributes:
        evapotranspiration,
        potential_evapotranspiration,
        daily_evapotranspiration,
        soil_evaporation,
        canopy_transpiration,
        interception_evaporation,
        relative_surface_wetness,
        green_canopy_fraction,
        plant_moisture_constraint,
        fractional_vegetation_cover,
        plant_temperature_constraint,
        epsilon,
        soil_net_radiation,
        canopy_net_radiation
        """
    #    RH =  clean_RH(np.array(AA.RH),0.0,1.0)
    # convert temperatures from kelvin to celsius
    air_temperature = air_temperature_K - 273
    air_temperature_mean = air_temperature_mean_K - 273
    air_temp_15day =air_temp_moving_window-273
    
    # scale water vapor pressure from Pa to kPa
    water_vapor_pressure_mean = water_vapor_pressure_mean_Pa * 0.001
    water_vapor_moving_window_kPa = water_vapor_moving_window * 0.001
    # calculate surface wetness values
    
    if verbose:
        print('\t\t\tcalculating surface wetness values')
    
    # calculate saturation vapor pressure in kPa from air temperature in celcius
    saturation_vapor_pressure = saturation_vapor_pressure_from_air_temperature(air_temperature_mean)

    saturation_vapor_pressure_15_day = saturation_vapor_pressure_from_air_temperature(air_temp_15day)

    # floor saturation vapor pressure at 1
    if floor_saturation_vapor_pressure:
        # saturation_vapor_pressure[saturation_vapor_pressure < 1] = 1
        saturation_vapor_pressure = np.clip(saturation_vapor_pressure, 1, None)
        saturation_vapor_pressure_15_day = np.clip(saturation_vapor_pressure_15_day, 1, None)

    # calculate vapor pressure deficit from water vapor pressure
    vapor_pressure_deficit = saturation_vapor_pressure - water_vapor_pressure_mean
    vapor_pressure_deficit_15_day = saturation_vapor_pressure_15_day - saturation_vapor_pressure_15_day
    
    # lower bound of vapor pressure deficit is zero, negative values replaced with nodata
    # vapor_pressure_deficit[vapor_pressure_deficit < 0] = np.nan
    vapor_pressure_deficit = np.where(vapor_pressure_deficit < 0, np.nan, vapor_pressure_deficit)
    vapor_pressure_deficit_15_day = np.where(vapor_pressure_deficit_15_day<0, np.nan, vapor_pressure_deficit_15_day);


    # calculate relative humidity from water vapor pressure and saturation vapor pressure
    relative_humidity = water_vapor_pressure_mean / saturation_vapor_pressure
    relative_humidity_15_day = water_vapor_moving_window_kPa / saturation_vapor_pressure_15_day
    
    # upper bound of relative humidity is one, results higher than one are capped at one
    #     relative_humidity[relative_humidity > 1] = 1
    relative_humidity[(relative_humidity<0.)|(relative_humidity>1.)]=np.nan;
    relative_humidity_15_day[(relative_humidity_15_day<0.)|(relative_humidity_15_day>1.)]=np.nan;
    #    relative_humidity = np.clip(relative_humidity, None, 1)
    #    relative_humidity_15_day = np.clip(relative_humidity_15_day, None, 1)

    # calculate relative surface wetness from relative humidity #<--- Can change this to bi-weekly if need be
    relative_surface_wetness = relative_humidity ** 4
    # reduce as in Mu et al., 2011
    #    relative_surface_wetness[relative_surface_wetness<0.7]=0.; This ommits unrealistic areas globally due to coarse resolution model is run at
    # limit cold temperatures to 0 as water would be frozen in these circumstances
    relative_surface_wetness[air_temperature_mean<0.]=0.;

    # calculate slope of saturation to vapor pressure curve Pa/K
    delta = delta_from_air_temperature(air_temperature)
    
    # calculate vegetation values
    
    if verbose:
        print('\t\t\tcalculating vegetation values')

    # calculate fAPAR from NDVI mean
    fAPAR = fAPAR_from_ndvi(ndvi_mean)

    # calculate fIPAR from NDVI mean
    fIPAR = fIPAR_from_ndvi(ndvi_mean)
    
    # calculate green canopy fraction (fg) from fAPAR and fIPAR, constrained between zero and one
    green_canopy_fraction = np.clip(fAPAR / fIPAR, 0, 1)
    
    # calculate plant moisture constraint (fM) from fraction of photosynthetically active radiation,
    # constrained between zero and one
    plant_moisture_constraint = np.clip(fAPAR / fAPARmax, 0, 1)


    # calculate soil moisture constraint from relative humidity and vapor pressure deficit,
    # reflects time-integrated component
    # constrained between zero and one
    # soil_moisture_constraint = np.clip(relative_humidity ** (vapor_pressure_deficit / BETA), 0, 1)
    soil_moisture_constraint = np.clip(relative_humidity_15_day ** (vapor_pressure_deficit_15_day / BETA), 0, 1)

    # take fractional vegetation cover from fraction of photosynthetically active radiation
    fractional_vegetation_cover = np.clip(fIPAR, 0, 1)
    
    # calculate plant temperature constraint (fT) from optimal phenology
    plant_temperature_constraint = exp(-(((air_temperature_mean - optimum_temperature) / optimum_temperature) ** 2))
    plant_temperature_constraint[air_temperature_mean<-5.] = 0.05; #< cold temperature constraint

    # calculate leaf area index
    leaf_area_index = -log(1 - fIPAR) * (1 / KPAR)
    
    # calculate delta / (delta + gamma)
    epsilon = delta / (delta + PSYCHROMETRIC_GAMMA)
    
    # soil evaporation
    
    if verbose:
        print('\t\t\tcalculating soil evaporation')

    # caluclate net radiation of the soil from leaf area index
    soil_net_radiation = net_radiation * exp(-KRN * leaf_area_index)


    # calculate instantaneous soil heat flux from net radiation and fractional vegetation cover
    #    soil_heat_flux = net_radiation * (0.05 + (1 - fractional_vegetation_cover) * 0.265)

    # validate
    # calculate the fraction of ground heat flux from net radiation
    gfrac = Gfrac_calc(CH, leaf_area_index);
    soil_heat_flux = net_radiation*gfrac;
    # floor soil heat flux at zero
    soil_heat_flux = np.clip(soil_heat_flux, 0, None)
    maximum_soil_heat_flux = 0.35 * soil_net_radiation
    # soil_heat_flux[soil_heat_flux > maximum_soil_heat_flux] = 0.35 * soil_net_radiation[soil_heat_flux > maximum_soil_heat_flux]
    soil_heat_flux = np.where(soil_heat_flux > maximum_soil_heat_flux, maximum_soil_heat_flux, soil_heat_flux)
    
    # calculate soil evaporation (LEs) from relative surface wetness, soil moisture constraint,
    # priestley taylor coefficient, epsilon = delta / (delta + gamma), net radiation of the soil,
    # and soil heat flux
    soil_evaporation = (relative_surface_wetness + soil_moisture_constraint * (1 - relative_surface_wetness)) * \
        PRIESTLEY_TAYLOR_ALPHA * epsilon * (soil_net_radiation - soil_heat_flux)
    # # replace missing soil evaporation with zero
    # soil_evaporation[np.isnan(soil_evaporation)] = 0
    soil_evaporation = np.where(np.isnan(soil_evaporation), 0, soil_evaporation)
    soil_evaporation[air_temperature_mean<0.]=0.
    # # floor soil evaporation at zero
    soil_evaporation = np.clip(soil_evaporation, 0, None)

    # canopy transpiration
    if verbose:
        print('\t\t\tcalculating canopy transpiration')
    # calculate net radiation of the canopy from net radiation of the soil
    canopy_net_radiation = net_radiation - soil_net_radiation
    # calculate canopy transpiration (LEc) from priestley taylor, relative surface wetness,
    # green canopy fraction, plant temperature constraint, plant moisture constraint,
    # epsilon = delta / (delta + gamma), and net radiation of the canopy
    canopy_transpiration = PRIESTLEY_TAYLOR_ALPHA *(1 - relative_surface_wetness) * green_canopy_fraction * plant_temperature_constraint *plant_moisture_constraint * epsilon *canopy_net_radiation
    # replace missing canopy transpiration with zero
    # canopy_transpiration[np.isnan(canopy_transpiration)] = 0
    canopy_transpiration = np.where(np.isnan(canopy_transpiration), 0, canopy_transpiration)
    # floor canopy transpiration at zero
    canopy_transpiration[canopy_transpiration < 0] = 0
    canopy_transpiration = np.clip(canopy_transpiration, 0, None)
    #***********************************************************************************************************
    #***********************************************************************************************************
    #***********************************************************************************************************
    #***********************************************************************************************************
    # NEW EVAPORATION COMPONENT
    rew_constraint = enforce_boundaries((SM-WP)/((FC*0.7)-WP),0,1);
    soil_evaporation_rew  = (relative_surface_wetness + rew_constraint * (1 - relative_surface_wetness)) * PRIESTLEY_TAYLOR_ALPHA * epsilon * (soil_net_radiation - soil_heat_flux)
    soil_evaporation_rew[soil_evaporation_rew < 0.] = 0.
    soil_evaporation_rew = np.where(np.isnan(soil_evaporation_rew), 0, soil_evaporation_rew)
    soil_evaporation_rew = np.clip(soil_evaporation_rew, 0, None)
    # New TREW CANOPY TRANSPIRATION
    potential_transpiration = PRIESTLEY_TAYLOR_ALPHA * epsilon * canopy_net_radiation
    rew_constraint_T_CH_PET = enforce_boundaries(fTREW_ch_pet_fun(FC,WP,LE_2_ETcm(potential_transpiration),CH,SM),0.0,1.0)
    canopy_transpiration_rew = PRIESTLEY_TAYLOR_ALPHA * (1 - relative_surface_wetness) * np.sqrt(green_canopy_fraction * plant_moisture_constraint * plant_temperature_constraint) * (rew_constraint_T_CH_PET)* epsilon * canopy_net_radiation
    canopy_transpiration_rew[canopy_transpiration_rew < 0] = 0
    canopy_transpiration_rew = np.where(np.isnan(canopy_transpiration_rew), 0, canopy_transpiration_rew)
    canopy_transpiration_rew = np.clip(canopy_transpiration_rew, 0, None)
    canopy_transpiration_rew[np.isnan(canopy_transpiration_rew)] = 0
    #***********************************************************************************************************
    #***********************************************************************************************************
    #***********************************************************************************************************
    #***********************************************************************************************************
    
    # interception evaporation
    
    if verbose:
        print('\t\t\tcalculating interception evaporation')

    # calculate interception evaporation (LEi) from relative surface wetness, priestley taylor coefficient,
    # epsilon = delta / (delta + gamma), and net radiation of the canopy
    interception_evaporation = relative_surface_wetness * PRIESTLEY_TAYLOR_ALPHA * epsilon * canopy_net_radiation

    # replace missing interception evaporation with zero
    # interception_evaporation[np.isnan(interception_evaporation)] = 0
    interception_evaporation = np.where(np.isnan(interception_evaporation), 0, interception_evaporation)
    
    # floor interception evaporation at zero
    # interception_evaporation[interception_evaporation < 0] = 0
    interception_evaporation = np.clip(interception_evaporation, 0, None)
    
    # combined evapotranspiration
    
    if verbose:
        print('\t\t\tcombining evapotranspiration')

    # combine soil evaporation (LEs), canopy transpiration (LEc), and interception evaporation (LEi)
    # into instantaneous evapotranspiration (LE)
    evapotranspiration = soil_evaporation + canopy_transpiration + interception_evaporation
    # limit evapotranspiration to net radiation
    # evapotranspiration[evapotranspiration > net_radiation] = net_radiation[evapotranspiration > net_radiation]
    evapotranspiration = np.clip(evapotranspiration, None, net_radiation)

    # replace infinite evapotranspiration with with nodata
    # evapotranspiration[np.isinf(evapotranspiration)] = np.nan
    evapotranspiration = np.where(np.isinf(evapotranspiration), np.nan, evapotranspiration)
    
    # remove negative and zero values from evapotranspiration
    # evapotranspiration[evapotranspiration <= 0] = np.nan
    evapotranspiration = np.where(evapotranspiration <= 0, np.nan, evapotranspiration)
    
    # daily evapotranspiration
    
    if verbose:
        print('\t\t\tcalculating daily evapotranspiration')

    # calculate evaporative fraction (EF) from evapotranspiration, net radiation, and soil heat flux
    evaporative_fraction = evapotranspiration / (net_radiation - soil_heat_flux)

    # calculate daily evapotranspiration from daily net radiation and evaporative fraction
    #     daily_evapotranspiration = daily_radiation * evaporative_fraction

    #***********************************************************************************************************
    #***********************************************************************************************************
    #***********************************************************************************************************
    #***********************************************************************************************************
    SMAPsoil_evaporation = soil_evaporation_rew;
    SMAPcanopy_transpiration = canopy_transpiration_rew
    SMAPevapotranspiration = SMAPsoil_evaporation + SMAPcanopy_transpiration + interception_evaporation
    # Mask to land area only
    mask = np.isnan(FC)
    SMAPevapotranspiration[mask]=np.nan
    SMAPsoil_evaporation[mask]=np.nan
    SMAPcanopy_transpiration[mask]=np.nan
    evapotranspiration[mask]=np.nan
    soil_evaporation[mask]=np.nan
    canopy_transpiration[mask]=np.nan
    interception_evaporation[mask]=np.nan
    #***********************************************************************************************************
    #***********************************************************************************************************
    #***********************************************************************************************************
    
    if verbose:
        print('\t\t\tcalculating potential evapotranspiration')

    # calculate potential evapotranspiration (pET) from priestley taylor coefficient,
    # epsilon = delta / (delta + gamma), net radiation, and soil heat flux
    potential_evapotranspiration = PRIESTLEY_TAYLOR_ALPHA * epsilon * (net_radiation - soil_heat_flux)

    ## Create xarray dataset to export and save as NETCDF File
    results = xr.Dataset({
                         'evapotranspiration': (['x', 'y'],           evapotranspiration),
                         'soil_evaporation': (['x', 'y'],             soil_evaporation),
                         'canopy_transpiration': (['x', 'y'],         canopy_transpiration),
                         'interception_evaporation':  (['x', 'y'],    interception_evaporation),
                         'potential_evapotranspiration': (['x', 'y'], potential_evapotranspiration),
                         'potential_transpiration':(['x', 'y'],       potential_transpiration),
                         'SMAPevapotranspiration': (['x', 'y'],       SMAPevapotranspiration),
                         'SMAPsoil_evaporation': (['x', 'y'],         SMAPsoil_evaporation),
                         'SMAPcanopy_transpiration': (['x', 'y'],     SMAPcanopy_transpiration)#,
#                         'soil_heat_flux':(['x','y'],                 soil_heat_flux)#,
#                         'relative_surface_wetness':  (['x', 'y'],    relative_surface_wetness),
#                         'green_canopy_fraction':  (['x', 'y'],       green_canopy_fraction),
#                         'plant_moisture_constraint':  (['x', 'y'],   plant_moisture_constraint),
#                         'fractional_vegetation_cover': (['x', 'y'],  fractional_vegetation_cover),
#                         'plant_temperature_constraint':  (['x', 'y'],plant_temperature_constraint),
#                         'epsilon':  (['x', 'y'],                     epsilon),
#                         'soil_net_radiation': (['x', 'y'],           soil_net_radiation),
#                         'canopy_net_radiation': (['x', 'y'],         canopy_net_radiation),
#                         'Soil_Moisture': (['x', 'y'],                SM)
                         },
                         coords={'lon': (['x', 'y'], lon),'lat': (['x', 'y'], lat)})
    return results
