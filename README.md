# ssj_latbin

A library which creates a machine-learning-ready Defense Meterology Satellite Program auroral particle precipitation (SSJ) dataset.

Data is reduced from raw data as follows:


Takes Level 1B NASA Common Data Format (CDF) SSJ files [[1]](#1) from [NASA CDAWeb](https://cdaweb.gsfc.nasa.gov/pub/data/dmsp/dmspf13/ssj/precipitating-electrons-ions/) and produces the latitude-and-orbit binned data in Apache Parquet format.

# References 
<a id="1">[1]</a>
Redmon, R. J., Denig, W. F., Kilcommons, L. M., and Knipp, D. J. (2017), New DMSP database of precipitating auroral electrons and ions, J. Geophys. Res. Space Physics, 122, 9056– 9067, doi:10.1002/2016JA023339. 
