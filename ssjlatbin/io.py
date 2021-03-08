# (C) 2021 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam Kilcommons - University of Colorado, Boulder - Colorado Center for Astrodynamics Research
# Mar 2021
import datetime
import glob
import toml
import numpy as np

def read_config(tomlfn):
    """Read the configuration file used by the rest of the application"""
    with open(tomlfn,'r') as f:
        tomldict = toml.loads(f.read())
    return tomldict

def ssjcdffn(dmsp_number,dt,config):
    """Find the full path to DMSP SSJ CDF file by looking
    recursively in directory ssj_cdf_root_dir set in config file
    """
    version = config['io']['ssj_cdf_version']
    cdfdir = config['io']['ssj_cdf_root_dir']
    fname = 'dmsp-f{:02d}_ssj_precipitating-electrons-ions_{}{:02d}{:02d}_v{}.cdf'.format(dmsp_number,dt.year,dt.month,dt.day,version)
    fullpath = glob.glob(cdfdir+'/**/'+fname,recursive=True)
    if len(fullpath)==0:
        raise IOError('No DMSP file {} found in {}'.format(fname,cdfdir))
    elif len(fullpath)>1:
        raise RuntimeError('Multiple matches {} this should not happen!'.format(fullpath))
    return fullpath[0]

def latbinned_flux_to_dataframe(t,lats,lats_dawn_dusk_flag,fluxes):
    """Columnify a fluxes array which originally had shape (len(t),len(lats))
    into a dataframe
    
    Parameters
    ----------
    t - list
        List of datetime.datetime objects, the times for each row
        of the fluxes array
    lats - np.array
        Array of center latitudes for each latitude bin
    lats_dawn_dusk_flag - np.array
        Array of same size as lats with value 1 if the corresponding lat
        is on the dawn side of the orbit, -1 if it is on the dusk side
    fluxes - np.array
        Array of electron or ion flux values
    
    Returns
    -------
    df - pd.DataFrame
        Pandas dataframe, indexed by time, with each column being
        the flux for a particular latitude bin. Column labels
        are either 'dusk' or 'dawn', then an underscore, then 
        the latitude as a string with 1 decimal place precision
    """
    datadict = {}
    for i_column in range(lats.size):
        if lats_dawn_dusk_flag[i_column] == 1:
            dawn_dusk = 'dawn'
        elif lats_dawn_dusk_flag[i_column] == -1:
            dawn_dusk = 'dusk'
        else:
            raise ValueError('Unexpected lats_dawn_dusk_flag value {}'.format(lats_dawn_dusk_flag[i_column]))
        col_label = '{}_{:.1f}'.format(dawn_dusk,lats[i_column])
        datadict[col_label] = fluxes[:,i_column]
    return pd.DataFrame(datadict,index=t)

def dataframe_to_latbinned_flux(binneddf,fluxvar):
    """Extract time and latitude 1D arrays and a flux 2D array for 
    use in plotting.
    
    Parameters
    ----------
    binneddf - pd.DataFrame
        Pandas dataframe with an orbit_number / latitude bin MultiIndex
        (as returned by ssjlatbin.latbin_pandas.bin_by_latitude).
    fluxvar - str
        Type of flux to extract (e.g. total_ele_number) 
        Must be a valid column name for binneddf
        
    Returns
    -------
    t - list
        List of datetime.datetime objects, the times for each row
        of the fluxes array
    lats - np.array
        Array of center latitudes for each latitude bin
    lats_dawn_dusk_flag - np.array
        Array of same size as lats with value 1 if the corresponding lat
        is on the dawn side of the orbit, -1 if it is on the dusk side
    fluxes - np.array
        Array of electron or ion flux values
    """
    t = [timestamp.to_pydatetime() for timestamp in binneddf.groupby('orbit_number')['time'].first()]
    df2d = binneddf[fluxvar].unstack()
    lats = np.full((len(df2d.columns),),np.nan)
    lats_dawn_dusk_flag = np.full((len(df2d.columns),),np.nan)
    fluxes = np.full((len(t),lats.size),np.nan)
    for icol,colname in enumerate(df2d.columns):
        dawn_dusk,latstr = colname.split('_')
        if dawn_dusk == 'dawn':
            lats_dawn_dusk_flag[icol]=1
        elif dawn_dusk == 'dusk':
            lats_dawn_dusk_flag[icol]=-1
        else:
            raise ValueError('Unexpected column name prefix {}'.format(dawn_dusk))
        lats[icol]=np.float(latstr)
        fluxes[:,icol] = df2d[colname].values
    return t,lats,lats_dawn_dusk_flag,fluxes