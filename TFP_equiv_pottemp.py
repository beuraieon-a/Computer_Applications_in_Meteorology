# Program for calculating and visualizing equivalent potential temperature

import datetime
import numpy as np
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.ndimage import gaussian_filter

###################
#### FUNCTIONS ####
###################

# FUNCTION 1: Data extraction
def extract_data(nc_file, variable_name, time, level):
    dataset = xr.open_dataset(nc_file)
    
    extracted_data = dataset[variable_name].sel(time=time, level=level)
    dataset.close()

    return extracted_data

# FUNCTION 2: Reading color values from a non-RGB colormap expressed in the range [1,0]
def read_colormap_file(colormap_file_path):
    with open(colormap_file_path, 'r') as open_file:
        color_lines = open_file.readlines()

    colors = []
    for color_line in color_lines:
        # Split the line into components
        r, g, b = color_line.split()
        colors.append((r, g, b))
    
    colors.reverse()
    
    return colors

#####################################
#### MAIN PROGRAM IMPLEMENTATION ####
#####################################

# Setting input file paths
temp_path = 'C:/Users/Brian/Desktop/RESEARCH/0_data/New folder/temp.nc'
spechum_path = 'C:/Users/Brian/Desktop/RESEARCH/0_data/New folder/spechum.nc'
uwind_path = 'C:/Users/Brian/Desktop/RESEARCH/0_data/New folder/uwind.nc'
vwind_path = 'C:/Users/Brian/Desktop/RESEARCH/0_data/New folder/vwind.nc'
colorfile_path = 'C:/Users/Brian/Desktop/RESEARCH/RdBl_diverge.rgb'

# Preparing the colormap
TFPcolors = read_colormap_file(colorfile_path)
TFPcolormap = mcolors.ListedColormap(TFPcolors)
TFPcolormap = TFPcolormap.reversed()

# Data extraction
time = '2022-12-18T00:00:00'
level = 925

temp = extract_data(temp_path, 'T', time, level)
spechum = extract_data(spechum_path, 'Q', time, level)
uwind = extract_data(uwind_path, 'U', time, level)
vwind = extract_data(vwind_path, 'V', time, level)

# Calculating dew-point temperature
dewtemp = mpcalc.dewpoint_from_specific_humidity(temp['level']*units.hPa,
                                                 temp*units.kelvin,
                                                 spechum*units('kg/kg'))
dewtemp = dewtemp.metpy.convert_units('kelvin')

# Calculating equivalent potential temperature
equiv_pottemp = mpcalc.equivalent_potential_temperature(temp['level']*units.hPa,
                                                        temp*units.kelvin,
                                                        dewtemp)

# Calculating gradient of equivalent potential temperature
grad_equiv_pottemp = mpcalc.gradient(equiv_pottemp, axes=(0,1)) # axes=(0,1) stands for latitude (0) and longitude (1)

# Calculating gradient magnitude of equivalent potential temperature
mag_grad_equiv_pottemp = np.sqrt((grad_equiv_pottemp[0] ** 2) + (grad_equiv_pottemp[1] ** 2))

# Calculating the negative gradient of the gradient magnitude of equivalent potential temperature
grad2_equiv_pottemp = mpcalc.gradient(mag_grad_equiv_pottemp, axes=(0,1))
neggrad2_equiv_pottemp = (grad2_equiv_pottemp[0] * (-1),
                          grad2_equiv_pottemp[1] * (-1))

# Calculating the unit vector in the direction of the gradient of equivalent potential temperature
unit_vector = (grad_equiv_pottemp[0] / mag_grad_equiv_pottemp,
               grad_equiv_pottemp[1] / mag_grad_equiv_pottemp)

# Calculating thermal front parameter (TFP) of equivalent potential temperature
# scaling the values in terms of x10^(-9)
TFP_equiv_pottemp = (neggrad2_equiv_pottemp[0] * unit_vector[0]) + (neggrad2_equiv_pottemp[1] * unit_vector[1])
TFP_equiv_pottemp = TFP_equiv_pottemp / (10**(-9))

# Initializing the map
x, y = np.meshgrid(temp['longitude'], temp['latitude'])
fig = plt.figure(figsize=(25., 25.), dpi=250)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([125, 175, 10, 55], ccrs.PlateCarree())

# Add gridlines
gridlines = ax.gridlines(draw_labels=False, xlocs=np.arange(-180, 181, 5),
                         ylocs=np.arange(-90, 91, 5), color='gray', linestyle='--')
xticks = np.arange(125, 176, 5)
yticks = np.arange(10, 56, 5)
ax.set_xticks(xticks, crs=ccrs.PlateCarree())
ax.set_yticks(yticks, crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)

# Add coastlines
ax.coastlines('10m', edgecolor='black', linewidth=2.5)

# Data smoothing
# TFP_equiv_pottemp = gaussian_filter(TFP_equiv_pottemp, sigma=2)

# Visualizing the thermal front parameter (TFP) of equivalent potential temperature
TFP_equiv_pottemp_visuals = ax.contourf(x, y, TFP_equiv_pottemp,
                                        transform=ccrs.PlateCarree(),
                                        levels=np.arange(-5, 5, 0.01),
                                        cmap=TFPcolormap,
                                        extend='both')

# Plotting wind velocity
ax.barbs(x, y, uwind, vwind, length=10, linewidth=1.5, regrid_shape=20,
         transform=ccrs.PlateCarree())

# Adding a colorbar
cax = fig.add_axes([0.92, 0.25, 0.02, 0.5])
norm = Normalize(vmin=-5, vmax=5)
colorbar = ColorbarBase(cax, cmap=TFPcolormap,
                        norm=norm, extend='both',
                        ticks=np.arange(-5, 6, 0.5))
colorbar.set_label('Thermal front parameter of θ$_e$ (K/m$^2$, ×10$^{-9}$)', fontsize=25)
colorbar.ax.tick_params(labelsize=20)

# Adding plot title and captions
time_truncated = time[:26]
time_obj = datetime.datetime.strptime(time_truncated, '%Y-%m-%dT%H:%M:%S')
formatted_date = time_obj.strftime('%d %B %Y')
formatted_time = time_obj.strftime('%H:%M')

title = f'{level}-hPa winds and thermal front parameter of equivalent potential temperature (θ$_e$) | ERA5'
caption = f'{formatted_time} UTC, {formatted_date}'
plt.suptitle(title, x=0.125, y=0.882, fontsize=32, ha='left', va='top')
plt.figtext(0.125, 0.863, caption, fontsize=25, ha='left', va='top')

plt.show()

# End of program