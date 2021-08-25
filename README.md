# ssj_latbin

A library which creates a machine-learning-ready Defense Meterology Satellite Program auroral particle precipitation (SSJ) dataset.

Data is reduced from raw data as follows:

1. The 19 SSJ energy channels are grouped in 'bands', by default "hard" (particle kinetic energy >1keV) and "soft" (particle kinetic energy <= 1keV) 
2. Uncertain (low detector count) measurements are discarded using flexible threshold (settable in configuration file)
3. The remaining measurements from all channels in a given band are summed 
4. Data is grouped by orbit (in the ML-ready dataset these would be the 'samples')
5. Data from each orbit is divided into latitude bins and the average flux in each bin is calculated (each bin is a 'feature' in the ML-ready dataset)
 
Takes Level 1B NASA Common Data Format (CDF) SSJ files (Redmon et al. 2017) from [NASA CDAWeb](https://cdaweb.gsfc.nasa.gov/pub/data/dmsp/dmspf13/ssj/precipitating-electrons-ions/) and produces the latitude-and-orbit binned data in Apache Parquet format.

## References

Redmon, R. J., Denig, W. F., Kilcommons, L. M., and Knipp, D. J. (2017), New DMSP database of precipitating auroral electrons and ions, J. Geophys. Res. Space Physics, 122, 9056– 9067, doi:10.1002/2016JA023339. 

