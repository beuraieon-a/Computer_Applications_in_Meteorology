# Program for implementing time series analysis of a meteorological parameter
# for a specific point location extracted from a gridded dataset

import numpy as np
import xarray as xr
from statsmodels.tsa.stattools import acf as autocorr
from statsmodels.graphics.tsaplots import plot_acf as plot_autocorr
from statsmodels.tsa.seasonal import seasonal_decompose
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime

###################
#### FUNCTIONS ####
###################

# FUNCTION 1: Data extraction
def extract_data(nc_file, variable_name, lat, lon):
    dataset = xr.open_dataset(nc_file)
    
    extracted_data = dataset[variable_name].sel(latitude=lat, longitude=lon, method='nearest')
    dataset.close()

    return extracted_data

#####################################
#### MAIN PROGRAM IMPLEMENTATION ####
#####################################

# Opening files
temp_file = 'C:/Users/Brian/Desktop/CSci135 (ComSci in Meteorology)/For practice/2mtemp_guiuan_1993_2023.nc'
mslp_file = 'C:/Users/Brian/Desktop/CSci135 (ComSci in Meteorology)/For practice/mslp_guiuan_1993_2023.nc'
tcrw_file = 'C:/Users/Brian/Desktop/CSci135 (ComSci in Meteorology)/For practice/totalcolumn_rainwater_guiuan_1993_2023.nc'

# Preparing wind variable names as inputs to FUNCTION 1 (extract_data)
temp_variable_name = 'VAR_2T'
mslp_variable_name = 'MSL'
tcrw_variable_name = 'TCRW'
lat = 11.033232682104137
lon = 125.7252134467548

# Data extraction
temp_data = extract_data(temp_file, temp_variable_name, lat, lon)
mslp_data = extract_data(mslp_file, mslp_variable_name, lat, lon)
tcrw_data = extract_data(tcrw_file, tcrw_variable_name, lat, lon)

# Unit conversion
temp_data = temp_data - 273.15 # Temperature, from Kelvin to Celsius
mslp_data = mslp_data/100 # Atmospheric pressure, pascal to hectopascal

# Time series visualization: Air temperature
plt.figure(figsize=(15., 10.), dpi=250)
temp_data.plot(color='red')
plt.title('Air temperature in Guiuan, Eastern Samar (1993–2023)', fontsize=20)
plt.xlabel('Time', fontsize=16)
plt.ylabel('Temperature (°C)', fontsize=16)
plt.xlim([datetime.date(1993, 1, 1), datetime.date(2023, 1, 1)])
plt.gca().xaxis.set_major_locator(mdates.YearLocator(2))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
plt.tick_params(axis='both', labelsize=15)
plt.show()

# Time series visualization: MSL atmospheric pressure
plt.figure(figsize=(15., 10.), dpi=250)
mslp_data.plot(color='green')
plt.title('Atmospheric pressure at mean sea level in Guiuan, Eastern Samar (1993–2023)', fontsize=20)
plt.xlabel('Time', fontsize=16)
plt.ylabel('Air pressure (hPa)', fontsize=16)
plt.xlim([datetime.date(1993, 1, 1), datetime.date(2023, 1, 1)])
plt.gca().xaxis.set_major_locator(mdates.YearLocator(2))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
plt.tick_params(axis='both', labelsize=15)
plt.show()

# Time series visualization: Total column rainwater
plt.figure(figsize=(15., 10.), dpi=250)
tcrw_data.plot(color='blue')
plt.title('Total column rainwater in Guiuan, Eastern Samar (1993–2023)', fontsize=20)
plt.xlabel('Time', fontsize=16)
plt.ylabel('Amount of rainwater (mm)', fontsize=16)
plt.xlim([datetime.date(1993, 1, 1), datetime.date(2023, 1, 1)])
plt.gca().xaxis.set_major_locator(mdates.YearLocator(2))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
plt.tick_params(axis='both', labelsize=15)
plt.show()

# Converting the meteorological data from xarray DataArray to pandas Series
temp_data_series = temp_data.to_series()
mslp_data_series = mslp_data.to_series()
tcrw_data_series = tcrw_data.to_series()

# Tests for autocorrelation
autocorr_temp = autocorr(temp_data_series)
autocorr_mslp = autocorr(mslp_data_series)
autocorr_tcrw = autocorr(tcrw_data_series)


# Visualization of autocorrelation
for data_series, title in zip([temp_data_series, mslp_data_series, tcrw_data_series], 
                              ['Air temperature', 'Atmospheric pressure at mean sea level', 'Total column rainwater']):
    fig, ax = plt.subplots(figsize=(15., 10.), dpi=250)
    plot_autocorr(data_series, ax=ax, lags=365) # Added lags=365 (annual, ~365 days)
    ax.set_title(f'Autocorrelation function:\n{title} in Guiuan, Eastern Samar (1993–2023)', fontsize=20)
    ax.set_xlabel('Lag', fontsize=16)
    ax.set_ylabel('Autocorrelation', fontsize=16)
    ax.tick_params(axis='both', labelsize=15)
    plt.show()

# Time series decomposition: Calculation and visualization
for data_series, title in zip([temp_data_series, mslp_data_series, tcrw_data_series], 
                              ['Air temperature', 'Atmospheric pressure at mean sea level', 'Total column rainwater']):
    tsDecompose = seasonal_decompose(data_series, model='additive', period=8760)
    # Added annual period=8760 since the data is hourly and there is 8760 hours in a year
    # (i.e. there are 24 hours in a day and ~365 days in a year -> 24*365)
    fig, axes = plt.subplots(4, 1, figsize=(15., 10.), dpi=250)
    tsDecompose.observed.plot(ax=axes[0], legend=False)
    axes[0].set_ylabel('Observed', fontsize=16)
    tsDecompose.trend.plot(ax=axes[1], legend=False)
    axes[1].set_ylabel('Trend', fontsize=16)
    tsDecompose.seasonal.plot(ax=axes[2], legend=False)
    axes[2].set_ylabel('Seasonal', fontsize=16)
    tsDecompose.resid.plot(ax=axes[3], legend=False)
    axes[3].set_ylabel('Residual', fontsize=16)
    axes[3].set_xlabel('Time', fontsize=16)
    fig.suptitle(f'Time series decomposition:\n{title} in Guiuan, Eastern Samar (1993–2023)', fontsize=20)
    plt.xticks([])
    plt.tick_params(axis='both', labelsize=15)
    plt.show()

# 1-D discrete fast Fourier transform (FFT) of the time series
for data_series, title in zip([temp_data_series, mslp_data_series, tcrw_data_series], 
                              ['Air temperature', 'Atmospheric pressure at mean sea level', 'Total column rainwater']):
    N = 365
    T = 1.0 / 365.0
    x = np.linspace(0.0, N*T, N)
    y = data_series
    xf = fftfreq(N, T)[1:N//2]
    yf = fft(y.values)
    fig, ax = plt.subplots(figsize=(15., 10.), dpi=250)
    plt.plot(xf, 2.0/N * np.abs(yf[1:N//2]))
    plt.grid()
    ax.set_title(f'1-D discrete Fourier transform:\n{title} in Guiuan, Eastern Samar (1993–2023)', fontsize=20)
    ax.set_xlabel('Frequency', fontsize=16)
    ax.set_ylabel('Amplitude', fontsize=16)
    ax.tick_params(axis='both', labelsize=15)
    plt.tight_layout()
    plt.show()

# End of program