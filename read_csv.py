# Program for reading data from a comma-separated values (CSV) file, e.g. IBTrACS TC data

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

IBTrACS = pd.read_csv('C:/Users/Brian/Desktop/Mete134 (Dynamic Meteorology II)/Meteorological data analysis and research/ibtracs.WP.list.v04r00.csv',
                      usecols=['SEASON', 'NAME', 'ISO_TIME', 'LAT', 'LON', 'TOKYO_GRADE', 'TOKYO_WIND', 'TOKYO_PRES'])

TCs_2013 = IBTrACS[IBTrACS['SEASON'] == 2013]

haiyan_track = TCs_2013[TCs_2013['NAME'] == 'HAIYAN']

fig = plt.figure(figsize=(25., 25.), dpi=250)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([90, 180, -10, 50], ccrs.PlateCarree())

gridlines = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 181, 10), ylocs=np.arange(-90, 91, 10), color='gray', linestyle='--')
gridlines.top_labels = False
gridlines.right_labels = False

ax.coastlines('10m', edgecolor='black', linewidth=2.5)
ax.plot(haiyan_track['LON'], haiyan_track['LAT'], 'ro', transform=ccrs.PlateCarree(), markersize=15)