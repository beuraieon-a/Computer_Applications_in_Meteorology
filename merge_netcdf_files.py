# Program for merging multiple NetCDF files into one file

import xarray

ds = xarray.open_mfdataset('C:/Users/Brian/Desktop/Future thesis works/Shearline/For merging/733967.Z.e5*.nc',
                           combine='nested', concat_dim='time')
ds.to_netcdf('C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/geopotential.nc')

print('\nDONE!')

# End of program