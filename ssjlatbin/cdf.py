# (C) 2021 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam Kilcommons - University of Colorado, Boulder - Colorado Center for Astrodynamics Research
# Mar 2021
import numpy as np
import pandas as pd
import datetime

from ssjlatbin.io import ssjcdffn
from ssjlatbin.fluxcalculations import integrate_flux
from ssjlatbin.tools import median_date,derivative

from functools import partial

from pycdflib.cdf import ReadOnlyCDF
from geospacepy.satplottools import simple_passes
from geospacepy.sun import solar_zenith_angle
from geospacepy.special_datetime import datetimearr2jd

def _define_ssj_dataframe_contents(config):
    #Define variables to load from CDF into dataframe which are already 1D
    dataframevar_to_cdfvar = config['dataframevar_to_cdfvar']

    soft_channels = config['soft_channels']
    hard_channels = config['hard_channels']
    all_channels = config['all_channels']

    #Define functions for getting total particle flux
    #from channel-by-channel flux
    totaling_functions = {}
    for channels,key in [(soft_channels,'soft'),(hard_channels,'hard'),(all_channels,'total')]:
        for fluxtype in ['energy','number']:
            totaling_functions[key+'_'+fluxtype] = partial(integrate_flux,
                                                            channels=channels,
                                                            energy_or_number=fluxtype)

    #Define variables which are not 1D in the CDF, along with a function that makes them 1D
    #so they can be added to the dataframe
    #(in this case the variables are [n_times x 19] and the function integrates along some or all of the 19 columns)
    for diff_flux_var,key in [('ELE_DIFF_ENERGY_FLUX','ele'),('ION_DIFF_ENERGY_FLUX','ion')]:
        for func_key,totalling_func in totaling_functions.items(): 
            dataframevar_to_cdfvar[key+'_'+func_key] = (diff_flux_var,totalling_func) 
    return dataframevar_to_cdfvar

def _read_ssj_cdf(ssjcdffn,config):
    """Read one spacecraft day of DMSP SSJ data into a dataframe"""
    dataframevar_to_cdfvar = _define_ssj_dataframe_contents(config)
    
    if 'uncertainty_tolerance' not in config['calculation']:
        uncertainty_tolerance = None
        print(config['calculation'])
    else:
        uncertainty_tolerance = config['calculation']['uncertainty_tolerance']

    cdf = ReadOnlyCDF(ssjcdffn)
    #read timestamps
    dts=cdf['Epoch']
    #read variables
    data = {}
    for dfvar,cdfvar_or_tup in dataframevar_to_cdfvar.items():
        if isinstance(cdfvar_or_tup,tuple):
            cdfvar = cdfvar_or_tup[0]
            to_1D_func = cdfvar_or_tup[1]
            if uncertainty_tolerance is None: #No uncertainty filtering
                data[dfvar]=to_1D_func(cdf[cdfvar])
            else:
                data[dfvar]=to_1D_func(cdf[cdfvar],
                                        diff_flux_rel_uncert=cdf[cdfvar+'_STD'],
                                        uncertainty_tolerance=uncertainty_tolerance)
        else:
            cdfvar = cdfvar_or_tup
            data[dfvar]=cdf[cdfvar]

    data['time']=dts
    data['solar_zenith_angle']=np.degrees(solar_zenith_angle(datetimearr2jd(dts),
                                                    data['glats'],
                                                    data['glons']))
    ssjcdfdf = pd.DataFrame(data,index=dts)
    return ssjcdfdf

def _read_current_previous_next_ssj_cdfs(prevcdffn,currcdffn,nextcdffn,config):
    """Read three consecutive spacecraft-days of data"""

    prevdf = _read_ssj_cdf(prevcdffn,config)
    currdf = _read_ssj_cdf(currcdffn,config)
    nextdf = _read_ssj_cdf(nextcdffn,config)
        
    if median_date(prevdf) != median_date(currdf)-datetime.timedelta(days=1):
        raise ValueError('Date of {} != date of {} - 1 day'.format(prevcdffn,
                                                                   currcdffn))
    if median_date(nextdf) != median_date(currdf)+datetime.timedelta(days=1):
        raise ValueError('Date of {} != date of {} + 1 day'.format(nextcdffn,
                                                                   currcdffn))

    ssjcdfdf = pd.concat([prevdf,currdf,nextdf])
    return ssjcdfdf


def _number_orbits(df,reference_date,latvar):
    """Find all equator crossings using latitudes in dataframe column latvar,
    marking each orbit with a integer, with orbit 0 being the first orbit
    ending after 00:00 UT of date refrence_date."""
    lats = df[latvar].values
    eqxinginds = simple_passes(lats,half_or_full_orbit='full')
    if isinstance(reference_date,datetime.datetime):
        refdate = reference_date.date()
    else:
        refdate = reference_date

    xing_zero = None
    for ixing,(ind_st,ind_ed) in enumerate(zip(eqxinginds[:-1],eqxinginds[1:])):
        xing_dates = df.iloc[ind_st:ind_ed].index.date
        if np.any([xing_date==refdate for xing_date in xing_dates]):
            xing_zero = ixing
            break
    if xing_zero is None:
        raise RuntimeError(('No data on {}'.format(reference_date)
                            +'in dataframe with dates {}'.format(df.index)))
    
    orbit_number = np.full_like(lats,np.nan)
    for ixing,(ind_st,ind_ed) in enumerate(zip(eqxinginds[:-1],eqxinginds[1:])):
        orbit_number[ind_st:ind_ed] = ixing-xing_zero

    return orbit_number

def _orbit_start_time(df):
    """Get a Series of the same length as the original data frame df.
    which has the first value from the time index for each orbit"""
    return df.groupby('orbit_number')['time'].transform('min')

def get_orbit_numbered_ssj_dataframe(dmsp_number,dt,config):
    """Get a dataframe of the SSJ data for one day, but ensuring the full orbit's data
    from the first and last orbits of the day is present from the previous and next days' data"""
    prevfn = ssjcdffn(dmsp_number,dt-datetime.timedelta(days=1),config)
    currfn = ssjcdffn(dmsp_number,dt,config)
    nextfn = ssjcdffn(dmsp_number,dt+datetime.timedelta(days=1),config)
    df = _read_current_previous_next_ssj_cdfs(prevfn,currfn,nextfn,config)
    df['orbit_number'] = _number_orbits(df,dt,'glats')
    df['orbit_start_time'] = _orbit_start_time(df)
    df['dglats'] = derivative(df['glats'].values)
    
    on_day = df.index.date==dt.date()
    day_orbits = np.unique(df[on_day]['orbit_number'].values)
    on_day_whole_orbits = np.logical_and(df['orbit_number']>=np.nanmin(day_orbits),
                                         df['orbit_number']<=np.nanmax(day_orbits))
    return df[on_day_whole_orbits].dropna().sort_index()

def get_orbit_numbered_ssj_range_dataframe(dmsp_number,dt_start,dt_end,config):
    """Get a dataframe of SSJ data for an arbitrary continuous time range"""
    dt = dt_start-datetime.timedelta(days=1)
    dfs = []
    while dt<dt_end:
        try:
            dfs.append(_read_ssj_cdf(ssjcdffn(dmsp_number,dt,config),config))
        except IOError:
            print(f'File not found for date {dt}')
        dt+=datetime.timedelta(days=1)
        
    df = pd.concat(dfs)
    df['orbit_number'] = _number_orbits(df,dt_start,'glats')
    df['orbit_start_time'] = _orbit_start_time(df)
    df['dglats'] = derivative(df['glats'].values)
    return df.dropna().sort_index()

