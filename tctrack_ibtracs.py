# Program for extracting TC track data from a comma-separated values (CSV) file
# from the IBTrACS archive

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

IBTrACS = pd.read_csv('C:/Users/Brian/Desktop/Mete134 (Dynamic Meteorology II)/Meteorological data analysis and research/ibtracs.WP.list.v04r00.csv',
                      usecols=['SEASON', 'NAME', 'ISO_TIME', 'LAT', 'LON', 'TOKYO_GRADE', 'TOKYO_WIND', 'TOKYO_PRES'])

TCs_2013 = IBTrACS[IBTrACS['SEASON'] == 2013]

haiyan_track = TCs_2013[TCs_2013['NAME'] == 'HAIYAN']
haiyan_track['TOKYO_WIND'] = haiyan_track['TOKYO_WIND'].replace(' ', '0')
haiyan_track['TOKYO_WIND'] = haiyan_track['TOKYO_WIND'].astype(int)

haiyan_LPA_track = haiyan_track[haiyan_track['TOKYO_GRADE'] == ' ']
haiyan_TD_track = haiyan_track[haiyan_track['TOKYO_GRADE'] == '2']
haiyan_TS_track = haiyan_track[haiyan_track['TOKYO_GRADE'] == '3']
haiyan_STS_track = haiyan_track[haiyan_track['TOKYO_GRADE'] == '4']
haiyan_TY_track = haiyan_track[(haiyan_track['TOKYO_GRADE'] == '5') & (haiyan_track['TOKYO_WIND'] < 100)]
haiyan_STY_track = haiyan_track[(haiyan_track['TOKYO_GRADE'] == '5') & (haiyan_track['TOKYO_WIND'] >= 100)]
haiyan_exTC_track = haiyan_track[haiyan_track['TOKYO_GRADE'] == '6']

fig = plt.figure(figsize=(25., 25.), dpi=250)
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([100, 170, 0, 30], ccrs.PlateCarree())

gridlines = ax.gridlines(draw_labels=True, xlocs=np.arange(-180, 181, 5), ylocs=np.arange(-90, 91, 5), color='gray', linestyle='--')
gridlines.top_labels = False
gridlines.right_labels = False

ax.coastlines('10m', edgecolor='black', linewidth=2.5)
ax.plot(haiyan_track['LON'], haiyan_track['LAT'], '-k', transform=ccrs.PlateCarree(), markersize=15)
ax.plot(haiyan_LPA_track['LON'], haiyan_LPA_track['LAT'], 'X', color='skyblue', transform=ccrs.PlateCarree(), markersize=10)
ax.plot(haiyan_TD_track['LON'], haiyan_TD_track['LAT'], 'o', color='blue', transform=ccrs.PlateCarree(), markersize=10)
ax.plot(haiyan_TS_track['LON'], haiyan_TS_track['LAT'], 'o', color='green', transform=ccrs.PlateCarree(), markersize=10)
ax.plot(haiyan_STS_track['LON'], haiyan_STS_track['LAT'], 'o', color='orange', transform=ccrs.PlateCarree(), markersize=10)
ax.plot(haiyan_TY_track['LON'], haiyan_TY_track['LAT'], 'o', color='red', transform=ccrs.PlateCarree(), markersize=10)
ax.plot(haiyan_STY_track['LON'], haiyan_STY_track['LAT'], 'o', color='purple', transform=ccrs.PlateCarree(), markersize=10)
ax.plot(haiyan_exTC_track['LON'], haiyan_exTC_track['LAT'], 's', color='gray', transform=ccrs.PlateCarree(), markersize=10)

plt.show()
