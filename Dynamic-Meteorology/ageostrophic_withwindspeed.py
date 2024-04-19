# Program for visualizing the actual, geostrophic and ageostrophic winds at a given pressure level
# Ver. 1 (actual wind speed in shaded contours in the background)

import numpy as np
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
import pandas as pd

# FUNCTION: Data extraction
def extract_data(nc_file, variable_name, pressure_level, lat_range, lon_range, time):
    dataset = xr.open_dataset(nc_file)
    
    variable_data = dataset[variable_name].sel(level=pressure_level, time=time)
    
    lats = dataset['lat'].values
    lons = dataset['lon'].values

    # Find indices corresponding to the specified lat/lon range
    lat_indices = np.where((lats >= lat_range[0]) & (lats <= lat_range[1]))[0]
    lon_indices = np.where((lons >= lon_range[0]) & (lons <= lon_range[1]))[0]

    extracted_data = variable_data[lat_indices, lon_indices]
    
    dataset.close()

    return extracted_data

# FUNCTION: Computing wind speed (magnitude of the wind velocity vectors) in km/h
def wind_speed_calc(u, v):
    return np.sqrt((u*3.6)**2 + (v*3.6)**2)

# Opening files
uwind_file = 'C:/Users/Brian/Desktop/Mete134 (Dynamic Meteorology II)/uwnd.2012.nc'
vwind_file = 'C:/Users/Brian/Desktop/Mete134 (Dynamic Meteorology II)/vwnd.2012.nc'
geopot_hgt_file = 'C:/Users/Brian/Desktop/Mete134 (Dynamic Meteorology II)/hgt.2012.nc'

# Preparing inputs to the function "extract_data": wind variable names and area of study (lat-lon range)
uwind_variable_name = 'uwnd'
vwind_variable_name = 'vwnd'
geopot_hgt_variable_name = 'hgt'
press_level = 300
lat_range = (20, 70)
lon_range = (230, 295)

# For today's demo, start at 2012-03-29T00:00:00 and end at 2012-04-08T00:00:00

start_time = pd.to_datetime('2012-03-29T00:00:00')
end_time = pd.to_datetime('2012-04-08T00:00:00')
times = pd.date_range(start_time, end_time, freq='6H')

for time in times:
    time_current = str(time)[:10] + 'T' + str(time)[11:]
    time_prev = str(time - pd.Timedelta(hours=6))[:10] + 'T' + str(time - pd.Timedelta(hours=6))[11:]
    
    # Extracting current data
    uwind_data = extract_data(uwind_file, uwind_variable_name, press_level, lat_range, lon_range, time_current)
    vwind_data = extract_data(vwind_file, vwind_variable_name, press_level, lat_range, lon_range, time_current)
    geopot_hgt_data = extract_data(geopot_hgt_file, geopot_hgt_variable_name, press_level, lat_range, lon_range, time_current)
    
    # Extracting data from 6 hours ago
    uwind_data_prev = extract_data(uwind_file, uwind_variable_name, press_level, lat_range, lon_range, time_prev)
    vwind_data_prev = extract_data(vwind_file, vwind_variable_name, press_level, lat_range, lon_range, time_prev)
    geopot_hgt_data_prev = extract_data(geopot_hgt_file, geopot_hgt_variable_name, press_level, lat_range, lon_range, time_prev)
    
    # Initializing the map
    x, y = np.meshgrid(geopot_hgt_data['lon'], geopot_hgt_data['lat'])
    fig = plt.figure(figsize=(25., 25.), dpi=250)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([230, 295, 20, 70], ccrs.PlateCarree())
    ax.coastlines('10m', edgecolor='black', linewidth=3)
    projection = ccrs.PlateCarree()
    
    # Add gridlines
    gridlines = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 181, 10), ylocs=np.arange(-90, 91, 10), color='gray', linestyle='--')
    gridlines.top_labels = False
    gridlines.right_labels = False
    
    # Visualizing 6-hour negative geopotential height tendencies (negative isallohypses),
    # i.e. amount of decrease in geopotential height
    negative_isallohypse_range = np.arange(-100, 0, 2)
    negative_isallohypse_visuals =  ax.contour(x, y, ((geopot_hgt_data/10)-(geopot_hgt_data_prev/10)), negative_isallohypse_range, colors='red', linewidths=2,
                                               linestyles='dashdot', transform=projection)
    kw_clabels = {'fontsize': 18, 'inline': True, 'inline_spacing': 5, 'fmt': '%i',
                  'rightside_up': True, 'use_clabeltext': True}
    plt.clabel(negative_isallohypse_visuals, **kw_clabels)
    
    # Visualizing 6-hour positive geopotential height tendencies (positive isallohypses),
    # i.e. amount of increase in geopotential height
    positive_isallohypse_range = np.arange(0, 100, 2)
    positive_isallohypse_visuals =  ax.contour(x, y, ((geopot_hgt_data/10)-(geopot_hgt_data_prev/10)), positive_isallohypse_range, colors='blue', linewidths=2,
                                               linestyles='dashdot', transform=projection)
    kw_clabels = {'fontsize': 18, 'inline': True, 'inline_spacing': 5, 'fmt': '%i',
                  'rightside_up': True, 'use_clabeltext': True}
    plt.clabel(positive_isallohypse_visuals, **kw_clabels)
    
    # Visualizing geopotential height
    geopot_hgt_range = np.arange(0, 1800, 5)
    geopot_hgt_visuals =  ax.contour(x, y, (geopot_hgt_data/10), geopot_hgt_range, colors='k', linewidths=2,
                                     linestyles='solid', transform=projection)
    kw_clabels = {'fontsize': 18, 'inline': True, 'inline_spacing': 5, 'fmt': '%i',
                  'rightside_up': True, 'use_clabeltext': True}
    plt.clabel(geopot_hgt_visuals, **kw_clabels)
    
    # Visualizing actual wind
    actual_wind_speed = wind_speed_calc(uwind_data, vwind_data)
    wind_speed_visuals = ax.contourf(x, y, actual_wind_speed, transform=ccrs.PlateCarree(), levels=np.arange(0, 400, 5), cmap='BuPu')
    
    # Calculating geostrophic wind
    geostrophic_wind_u, geostrophic_wind_v = mpcalc.geostrophic_wind((geopot_hgt_data) * units.m)
    
    # Calculating ageostrophic wind
    ageostrophic_wind_u, ageostrophic_wind_v = mpcalc.ageostrophic_wind((geopot_hgt_data) * units.m, uwind_data, vwind_data)
    
    # Set up parameters for quiver plot. The slices below are used to subset the data (here
    # taking every 4th point in x and y). The quiver_kwargs are parameters to control the
    # appearance of the quiver so that they stay consistent between the calls.
    quiver_slices = (slice(None, None, 2), slice(None, None, 3))
    quiver_kwargs = {'headlength': 4, 'headwidth': 3, 'angles': 'uv', 'scale_units': 'xy',
                     'scale': 6, 'width': 0.005}
    
    # Plot the wind vectors
    actual_wind = ax.quiver(x[quiver_slices], y[quiver_slices],
                            uwind_data[quiver_slices], vwind_data[quiver_slices],
                            color='lime', **quiver_kwargs)
    
    geostrophic_component = ax.quiver(x[quiver_slices], y[quiver_slices],
                                      geostrophic_wind_u[quiver_slices], geostrophic_wind_v[quiver_slices],
                                      color='red', **quiver_kwargs)
    
    ageostrophic_component = ax.quiver(x[quiver_slices], y[quiver_slices],
                                       ageostrophic_wind_u[quiver_slices], ageostrophic_wind_v[quiver_slices],
                                       color='blue', **quiver_kwargs)
    
    # Adding colorbar
    cax = fig.add_axes([0.92, 0.25, 0.02, 0.5])
    norm = Normalize(vmin=0, vmax=400)
    colorbar = ColorbarBase(cax, cmap='BuPu', norm=norm, extend='max', ticks=np.arange(0, 401, 10))
    colorbar.set_label('Actual wind speed (km/h)')
    
    # Adding plot title and captions
    title = 'Wind flow at the ' + str(press_level) + ' hPa level, ' + time_current[:10] + ' ' + time_current[11:16] + ' UTC'
    caption = 'Geopotential height (dam, black solid isohypses), 6-hour geopotential height change (dam, red/blue dashed-dotted isallohypses), actual wind speed (km/h, shading),\nactual wind velocity (m/s, lime-green arrow), geostrophic wind (m/s, red arrow), ageostrophic wind (m/s, blue arrow)'
    plt.suptitle(title, x=0.125, y=0.839, fontsize=35, ha='left', va='top')
    plt.figtext(0.125, 0.820, caption, fontsize=15, ha='left', va='top', style='italic')
    
    plt.show()

    # Printing the status
    print('\n', time,'>> DONE!')

# End of program
