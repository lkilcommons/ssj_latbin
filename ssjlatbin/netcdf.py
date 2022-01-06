import datetime,os
from collections.abc import Mapping #base class used by ReadOnlyCDF

import numpy as np
from netCDF4 import Dataset
from dateutil.relativedelta import relativedelta

class ReadOnlyConvertedNC(Mapping):
    """Class providing dict-like reading syntax for 
    NetCDF files which were created from CDF files 
    using the NASA cdf_to_netcdf tool"""

    def __init__(self,fn):
        self.fn=fn
        self.ds = Dataset(fn,'r')
        self._keys = list(self.ds.variables.keys())
        self.epocharr_to_datetime = np.vectorize(self._ms_since_0AD_to_datetime)
            
    @staticmethod
    def _ms_since_0AD_to_datetime(ms_since_0AD):
        """CDF files which have been converted to netCDF using 
        the NASA cdf_to_netcdf tool appear to convert EPOCH CDF type
        into a float64 variable which is milliseconds since 0 AD"""
        
        dt_1AD = datetime.datetime(1,1,1,0,0) #0 AD is not defined as a Python datetime
        s_since_0AD=ms_since_0AD/1000.

        return dt_1AD+datetime.timedelta(seconds=s_since_0AD)-relativedelta(years=1)

    def __str__(self):
        keystr = '\n'.join([key for key in self])
        return f'Converted CDF NetCDF File {self.fn} with variables:\n{keystr}'

    def __len__(self):
        return self._keys.__len__()

    def __contains__(self,key):
        return self._keys.__contains__(key)

    def __iter__(self):
        return self._keys.__iter__()

    def __getitem__(self,key):
        """netCDF4 Datasets return masked arrays, so the .data is required 
        to get at the bare numpy array"""
        if key=='Epoch':
            return self.epocharr_to_datetime(self.ds['Epoch'][:].data)    
        else:
            return self.ds[key][:].data