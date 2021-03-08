# (C) 2021 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam Kilcommons - University of Colorado, Boulder - Colorado Center for Astrodynamics Research
# Mar 2021
import numpy as np
import pandas as pd

def median_date(df):
    """Determine median date from a possibly unsorted dataframe with a datetime index"""
    return pd.to_datetime(df.sort_index().index[len(df.index)//2]).date()

def derivative(y):
    """Returns 1st order finite difference with the first value
    duplicated to keep the array the same size"""
    dy = np.diff(y)
    dy = np.concatenate([dy,np.array([dy[0]])],axis=0) #make same length as timeseries
    return dy