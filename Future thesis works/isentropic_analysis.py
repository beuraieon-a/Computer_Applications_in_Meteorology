# Program for implementing isentropic analysis

import pandas as pd
import numpy as np
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

###################
#### FUNCTIONS ####
###################

# FUNCTION 1: Data extraction
def extract_data(nc_file, variable_name, time):
    dataset = xr.open_dataset(nc_file)
    
    if variable_name == 'VAR_100U' or variable_name == 'VAR_100V':
        extracted_data = dataset[variable_name].sel(time=time)
        dataset.close()
    else:
        extracted_data = dataset[variable_name].sel(time=time, level=slice(100, 1000))
        dataset.close()

    return extracted_data

#####################################
#### MAIN PROGRAM IMPLEMENTATION ####
#####################################

# Setting input file paths
temp_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/temperature.nc'
uwind_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/u-wind.nc'
vwind_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/v-wind.nc'
geopot_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/geopotential.nc'
spechum_path = 'C:/Users/Brian/Desktop/Future thesis works/Shearline/Pressure level/specific_humidity.nc'

# isentropic_level_list = [290.0, 296.0, 300.0, 302.0, 304.0, 306.0] * units.kelvin
isentropic_level_list = [306.0] * units.kelvin # Best choices for tropics
time_range = pd.date_range(start='2022-12-21T00:00:00', end='2022-12-27T00:00:00', freq='h')

for level in isentropic_level_list:
    for time in time_range:
        time = str(time)
        
        # Data extraction
        temp_data = extract_data(temp_path, 'T', time)
        uwind_data = extract_data(uwind_path, 'U', time)
        vwind_data = extract_data(vwind_path, 'V', time)
        geopot_data = extract_data(geopot_path, 'Z', time)
        spechum_data = extract_data(spechum_path, 'Q', time)
        
        # Deriving geopotential height values from the geopotential data
        geopot_hgt_data = geopot_data / 9.80665
        
        # Combining all variable data into a single DataArray
        data = xr.Dataset({
            'temp': temp_data,
            'uwind': uwind_data,
            'vwind': vwind_data,
            'geopot_hgt': geopot_hgt_data,
            'spechum': spechum_data,
            })
        
        data = data.metpy.parse_cf().squeeze()
        
        isentropic_data = mpcalc.isentropic_interpolation_as_dataset(
            isentropic_level_list,
            data['temp'],
            data['uwind'],
            data['vwind'],
            data['spechum'],
            data['geopot_hgt']
        )
        
        isentropic_data['relhum'] = mpcalc.relative_humidity_from_specific_humidity(
            isentropic_data['pressure'],
            isentropic_data['temperature'],
            isentropic_data['spechum']
        ).metpy.convert_units('percent')
        
        isentropic_data['uwind'] = isentropic_data['uwind'].metpy.convert_units('kt')
        isentropic_data['vwind'] = isentropic_data['vwind'].metpy.convert_units('kt')
        isentropic_data['spechum'] = (isentropic_data['spechum']) * 1000
        
        crs = ccrs.PlateCarree()
        lat = isentropic_data['pressure'].metpy.latitude
        lon = isentropic_data['pressure'].metpy.longitude
        
        bounds = [(100, 160, 0, 40)]
        
        fig = plt.figure(figsize=(17., 12.), dpi=250)
        ax = fig.add_subplot(1, 1, 1, projection=crs)
        ax.set_extent(*bounds, crs=ccrs.PlateCarree())
        
        """
        spechum_contourf = ax.contourf(lon, lat, isentropic_data['spechum'].sel(isentropic_level=level),
                 levels=np.arange(0, 20, 0.001), cmap='gist_earth_r', transform=ccrs.PlateCarree())
        """
        
        relhum_contourf = ax.contourf(lon, lat, isentropic_data['relhum'].sel(isentropic_level=level),
                 levels=np.arange(0, 130, 0.01), cmap='gist_earth_r', transform=ccrs.PlateCarree())
        
        ax.coastlines('10m', edgecolor='black', linewidth=2.5)
        
        isentropic_contour_levels = np.arange(0, 1000, 25)
        isentropic_contour = ax.contour(lon, lat, isentropic_data['pressure'].sel(isentropic_level=level),
                isentropic_contour_levels, colors='k', linewidths=1.0, linestyles='solid',
                transform=ccrs.PlateCarree())
        isentropic_contour.clabel(fontsize=10, inline=1, inline_spacing=7, fmt='%i', rightside_up=True,
                  use_clabeltext=True)
        
        """
        spechum_colorbar = fig.colorbar(spechum_contourf, orientation='horizontal', aspect=65, shrink=0.5, pad=0.05,
                                        extendrect='False', ticks=np.arange(0, 20, 2))
        spechum_colorbar.set_label('Specific Humidity (g/kg)', size='x-large')
        """
        
        relhum_colorbar = fig.colorbar(relhum_contourf, orientation='horizontal', aspect=65, shrink=0.5, pad=0.05,
                                        extendrect='False', ticks=np.arange(0, 130, 5))
        relhum_colorbar.set_label('Relative Humidity (%)', size='x-large')
        
        # Plotting wind barbs
        ax.barbs(lon.values, lat.values, isentropic_data['uwind'].sel(isentropic_level=level).values,
                 isentropic_data['vwind'].sel(isentropic_level=level).values, length=6,
                 regrid_shape=20, transform=ccrs.PlateCarree())

        # Setting the plot titles
        """
        ax.set_title(f'{level:~.0f} Isentropic Pressure (hPa), Wind (kt) and Specific Humidity (g/kg) | '
                     f'Valid: {isentropic_data["time"].dt.strftime("%Y-%m-%d %H:%MZ").item()}\n', loc='left')
        """
        ax.set_title(f'{level:~.0f} Isentropic Pressure (hPa), Wind (kt) and Relative Humidity (%) | '
                     f'Valid: {isentropic_data["time"].dt.strftime("%Y-%m-%d %H:%MZ").item()}\n', loc='left')
        fig.tight_layout()
        
        plt.show()
        
        print('\n', level, 'isentropic analysis @' , time,'>> DONE!\n')

# End of program