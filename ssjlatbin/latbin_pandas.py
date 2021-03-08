import numpy as np
import pandas as pd

def _dawn_dusk(hemi,asc_desc):
    """Determine if the spacecraft was in the dawn or
    dusk sector from which hemisphere it is in an whether
    it is ascending (going north) or descending (going south)
    Returns 1 for dawn, -1 for dusk"""
    if hemi=='N' and asc_desc=='asc':
        return 'dawn'
    elif hemi=='N' and asc_desc=='desc':
        return 'dusk'
    elif hemi=='S' and asc_desc=='desc':
        return 'dusk'
    elif hemi=='S' and asc_desc=='asc':
        return 'dawn'
    else:
        raise ValueError(('Invalid hemisphere {} (N or S) or'.format(hemi)
                         +'ascending/descending {} (asc or desc)'.format(asc_desc)))

def define_latbins(delta_lat,max_lat):
    lat_bin_edges = {'N':{},'S':{}}
    lat_bin_edges['N']['asc'] = np.arange(0,max_lat+delta_lat,delta_lat)
    lat_bin_edges['N']['desc'] = np.arange(max_lat,0-delta_lat,-delta_lat)
    lat_bin_edges['S']['desc'] = np.arange(0,-1*max_lat-delta_lat,-delta_lat)
    lat_bin_edges['S']['asc'] = np.arange(-1*max_lat,0+delta_lat,delta_lat)
    return lat_bin_edges

def latbin_label(lat1,lat2,hemi,asc_desc):
    return _dawn_dusk(hemi,asc_desc)+'_{:.1f}'.format((lat1+lat2)/2)

def bin_by_latitude(orbit_numbered_ssj_dataframe,config,latvar='glats'):
    """Extract each orbit as one row in an array, binning the data
    into latitude bins of width delta_lat degrees to get a constant
    number of columns for each orbit"""
    delta_lat=config['latbin']['delta_lat']
    max_lat=config['latbin']['max_lat']
    
    df = orbit_numbered_ssj_dataframe
    df['hemi']=pd.cut(df['glats'],[-91.,0.,91.],labels=['S','N'])
    df['asc_desc']=pd.cut(df['dglats'],[-np.inf,0,np.inf],labels=['desc','asc'])

    categories = []
    conditions = []
    lat_bin_edges = define_latbins(delta_lat,max_lat)
    lat_bin_centers = []
    for hemi,asc_desc in [('N','asc'),('N','desc'),('S','desc'),('S','asc')]:
        edges = lat_bin_edges[hemi][asc_desc]
        for lat1,lat2 in zip(edges[:-1],edges[1:]):
            lat_bin_centers.append((lat1+lat2)/2)
            categories.append(latbin_label(lat1,lat2,hemi,asc_desc))
            slat = np.nanmin([lat1,lat2])
            elat = np.nanmax([lat1,lat2])
            bin_mask = np.logical_and(df[latvar]>=slat,df[latvar]<elat)
            hemi_mask = df['hemi']==hemi
            asc_mask = df['asc_desc']==asc_desc
            conditions.append((bin_mask & hemi_mask & asc_mask))

    #Returns an array with length the same number of rows as df
    #of strings which are the labels for latitude bins
    latbinned = np.select(conditions,categories)

    #The bins do have an order because we want them to plot in a particular
    #order
    df['latbin']=pd.Categorical(latbinned,categories=categories,ordered=True)
    
    #Store the time as datetime64 for the averaging operation   
    df['time']=df.index.values.astype(np.int64)
    
    binneddf = df.groupby(['orbit_number','latbin']).mean()
    binneddf['time'] = pd.to_datetime(binneddf['time'])
    return binneddf

