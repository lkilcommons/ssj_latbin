# (C) 2021 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam Kilcommons - University of Colorado, Boulder - Colorado Center for Astrodynamics Research
# Mar 2021
import numpy as np
import pandas as pd
from dateutil.relativedelta import *
def _dawn_dusk(hemi,asc_desc):
    """Determine if the spacecraft was in the dawn or
    dusk sector from which hemisphere it is in an whether
    it is ascending (going north) or descending (going south)
    Returns 1 for dawn, -1 for dusk"""
    if hemi=='N' and asc_desc=='asc':
        return 1
    elif hemi=='N' and asc_desc=='desc':
        return -1
    elif hemi=='S' and asc_desc=='desc':
        return -1
    elif hemi=='S' and asc_desc=='asc':
        return 1
    else:
        raise ValueError(('Invalid hemisphere {} (N or S) or'.format(hemi)
                         +'ascending/descending {} (asc or desc)'.format(asc_desc)))
        
def _ascending_descending_masks(glats):
    """Calculate boolean arrays (masks) for when the spacecraft is ascending (going northward),
    or descending (going southward) from the timeseries of it's geographic latitudes"""
    dlat = np.diff(glats)
    dlat = np.concatenate([dlat,np.array([dlat[0]])],axis=0) #make same length as timeseries
    masks = {}
    masks['asc'] = dlat>0 #spacecraft is ascending (going north)
    masks['desc'] = dlat<0 #spacecraft is descending (going south)
    return masks
        
def bin_by_latitude(orbit_numbered_ssj_dataframe,var_to_bin,latvar='glats',delta_lat=5,max_lat=80):
    """Extract each orbit as one row in an array, binning the data
    into latitude bins of width delta_lat degrees to get a constant
    number of columns for each orbit"""
    df = orbit_numbered_ssj_dataframe
    glats = df['glats'].values
    asc_desc_masks = _ascending_descending_masks(glats)
    lats = df[latvar].values
    y = df[var_to_bin].values
    if np.mod(max_lat*4./delta_lat,delta_lat)!=0:
        raise ValueError('Non-integer number of latitude bins with maxlat {}, delta lat {}'.format(max_lat,delta_lat))
    n_bins = int(np.round(max_lat*4./delta_lat))
    orbitnums = np.unique(df['orbit_number'].values)
    n_orbits = orbitnums.size
    binned_y = np.full((n_orbits,n_bins),np.nan)
    orbit_ts = np.full((n_orbits,),np.nan,dtype=object)
    bin_lats = np.full((n_bins,),np.nan)
    dawn_dusk_flag = np.full((n_bins),np.nan)
    lat_bin_edges = {'N':{},'S':{}}
    lat_bin_edges['N']['asc'] = np.arange(0,max_lat+delta_lat,delta_lat)
    lat_bin_edges['N']['desc'] = np.arange(max_lat,0-delta_lat,-delta_lat)
    lat_bin_edges['S']['desc'] = np.arange(0,-1*max_lat-delta_lat,-delta_lat)
    lat_bin_edges['S']['asc'] = np.arange(-1*max_lat,0+delta_lat,delta_lat)
    for i_orbit,orbitnum in enumerate(orbitnums):
        orbit_mask = df['orbit_number']==orbitnum
        orbit_ts[i_orbit]=pd.to_datetime(df[orbit_mask].index[0])
        print("Orbit {}: {}".format(orbitnum,orbit_ts[i_orbit]))
        i_lat_bin = 0
        for hemi,asc_desc in [('N','asc'),('N','desc'),('S','desc'),('S','asc')]:
            #Treat each quarter of the orbit individually
            asc_desc_mask = asc_desc_masks[asc_desc]
            for lat_start,lat_end in zip(lat_bin_edges[hemi][asc_desc][:-1],lat_bin_edges[hemi][asc_desc][1:]):
                if i_orbit == 0: #Only need to do this once
                    bin_lats[i_lat_bin]=(lat_start+lat_end)/2
                    dawn_dusk_flag[i_lat_bin]=_dawn_dusk(hemi,asc_desc)

                if lat_start < lat_end:
                    lat_mask = np.logical_and(lats>lat_start,lats<=lat_end)
                else:
                    lat_mask = np.logical_and(lats<lat_start,lats>=lat_end)
                mask = (orbit_mask & asc_desc_mask & lat_mask)
                if np.count_nonzero(mask)==0:
                    print('No data {} {}-{}'.format(asc_desc,lat_start,lat_end))
                
                binned_y[i_orbit,i_lat_bin] = np.nanmean(y[mask])
                i_lat_bin+=1
                
    return orbit_ts,bin_lats,dawn_dusk_flag,binned_y
