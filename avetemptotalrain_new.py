# Program for calculating average temperature and accumulated rainfall with application of NumPy, pandas and matplotlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime

print('\n<< Average Temperature and Accumulated Rainfall Calculator >>\nThis program calculates the monthly average temperature and monthly accumulated rainfall using daily temperature and rainfall amount data for a given month.\n\nNOTE: This program requires the input data to be stored in a plain text file with three columns (date, temperature, rainfall amount) separated with a single space.')

# Preparation for printing user input and output: Units and date

# Temprature unit
temp_unit = input('\nChoose unit for temperature (c | f | k): ')

while(temp_unit != 'c' and temp_unit != 'C' and temp_unit != 'f' and temp_unit != 'F' and temp_unit != 'k' and temp_unit != 'K'):
    temp_unit = input('\nINVALID INPUT! Choose unit for temperature (c | f | k): ')

if(temp_unit == 'c' or temp_unit == 'C' or temp_unit == 'f' or temp_unit == 'F'):
    temp_unit = '°' + temp_unit.upper()
else:
    temp_unit = temp_unit.upper()

# Rainfall unit
rain_unit = input('\nChoose unit for rainfall amount (mm | in): ')

while(rain_unit != 'mm' and rain_unit != 'in'):
    rain_unit = input('\nINVALID INPUT! Choose unit for rainfall amount (mm | in): ')

# Month and year
year = int(input('\nEnter the year (YYYY): '))
month = int(input('\n1 - January     2 - February     3 - March\n4 - April     5 - May     6 - June\n7 - July     8 - August     9 - September\n10 - October     11 - November     12 - December\nEnter the month: '))

if(month == 1):
    month_name = 'January'
elif(month == 2):
    month_name = 'February'
elif(month == 3):
    month_name = 'March'
elif(month == 4):
    month_name = 'April'
elif(month == 5):
    month_name = 'May'
elif(month == 6):
    month_name = 'June'
elif(month == 7):
    month_name = 'July'
elif(month == 8):
    month_name = 'August'
elif(month == 9):
    month_name = 'September'
elif(month == 10):
    month_name = 'October'
elif(month == 11):
    month_name = 'November'
else:
    month_name = 'December'

# Input data

dates = []
temp = []
rain = []

input_data_file_path = input('\nEnter the absolute path to the TXT file containing daily temperature and rainfall data:\n')

with open(input_data_file_path, 'r') as data_file:
    for data_line in data_file.readlines():
        date_data, temp_data, rain_data = data_line.strip().split()
        dates.append(date_data)
        temp.append(float(temp_data))
        rain.append(float(rain_data))

dates_array = np.array(dates)
temp_array = np.array(temp)
rain_array = np.array(rain)

print('\n----------------\n\nINPUT DATA (temperature & rainfall)\n')

for x in range(len(dates)):
    print(dates_array[x] + '   ' + str(temp_array[x]) + temp_unit + '    ' + str(rain_array[x]) + ' ' + rain_unit + ' ')

input('\nPress ENTER to proceed to calculation.')

# Calculation

temp_sum = np.sum(temp_array)
temp_ave = np.mean(temp_array)
temp_stdev = np.std(temp_array)
temp_var = np.var(temp_array)

rain_sum = np.sum(rain_array)
rain_ave = np.mean(rain_array)
rain_stdev = np.std(rain_array)
rain_var = np.var(rain_array)

print('\n----------------\n\nOUTPUT\n\n>> Temperature in ' + month_name + ' ' + str(year))
print('[Sum: ' + str(temp_sum) + ']\nMonthly average temperature: ' + str(round(temp_ave, 2)) + temp_unit + '\nStandard deviation: ' + str(round(temp_stdev, 2)) + temp_unit + '\nVariance: ' + str(round(temp_var, 2)) + temp_unit)
print('\n\n>> Rainfall in ' + month_name + ' ' + str(year))
print('Monthly accumulated rainfall: ' + str(rain_sum) + ' ' + rain_unit + '\nMonthly average rainfall: ' + str(round(rain_ave, 2)) + ' ' + rain_unit + '\nStandard deviation: ' + str(round(rain_stdev, 2)) + ' ' + rain_unit + '\nVariance: ' + str(round(rain_var, 2)) + ' ' + rain_unit)

input('\nPress ENTER to proceed.')

data = {
    'Date': dates,
    'Temperature (' + temp_unit + ')' : temp,
    'Rainfall (' + rain_unit + ')': rain
}

print('\n----------------\n\nINPUT DATA (temperature & rainfall)\n')

data_DataFrame = pd.DataFrame(data)

print(data_DataFrame)

# Preparations for visualization

datesNum_list = [mdates.date2num(datetime.strptime(date, '%m/%d/%Y')) for date in dates_array]

plt.rcParams['font.size'] = 16

title_temp = 'Temperature in ' + month_name + ' ' + str(year)
title_rain = 'Rainfall in ' + month_name + ' ' + str(year)
label_temp = 'Temperature (' + temp_unit + ')'
label_rain = 'Rainfall (' + rain_unit + ')'

# Visualization: Temperature line graph

plt.figure(figsize=(15., 10.), dpi=250)
plt.plot_date(datesNum_list, temp_array, 'r-', xdate=True)
plt.title(title_temp)
plt.ylabel(label_temp)
plt.xlabel('Date')
plt.gca().xaxis.set_major_locator(mdates.DayLocator())
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
plt.gcf().autofmt_xdate()

plt.show()

# Visualization: Rainfall bar graph

plt.figure(figsize=(15., 10.), dpi=250)
plt.bar(datesNum_list, rain_array, color='b')
plt.title(title_rain)
plt.ylabel(label_rain)
plt.xlabel('Date')
plt.gca().xaxis.set_major_locator(mdates.DayLocator())
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
plt.gcf().autofmt_xdate()

plt.show()

# Visualization: Box plots

fig, ax = plt.subplots(1, 2, figsize=(15., 10.), dpi=250)

ax[0].boxplot(temp_array)
ax[0].set_title(title_temp)
ax[0].set_ylabel(label_temp)

ax[1].boxplot(rain_array)
ax[1].set_title(title_rain)
ax[1].set_ylabel(label_rain)

plt.show()

# End of program