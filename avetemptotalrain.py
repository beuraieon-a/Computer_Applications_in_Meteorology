# Program for calculating average temperature and accumulated rainfall

print('\n<< Average Temperature and Accumulated Rainfall Calculator >>\nThis program calculates the monthly average temperature and monthly accumulated rainfall using daily temperature and rainfall amount data for a given month.\n\nNOTE: This program requires the input data to be stored in two separate plain text files (one for each of the two variables), wherein a data value is indicated per line.')

# Preparation for printing user input and output: Units and date

# Temprature unit
temp_unit = input('\nChoose unit for temperature (c | f | k): ')

while(temp_unit != 'c' and temp_unit != 'C' and temp_unit != 'f' and temp_unit != 'F' and temp_unit != 'k' and temp_unit != 'K'):
    temp_unit = input('\nINVALID INPUT! Choose unit for temperature (c | f | k): ')

if(temp_unit == 'c' or temp_unit == 'C'):
    temp_unit = temp_unit.upper()
    temp_unit_name = 'Celsius'
elif(temp_unit == 'f' or temp_unit == 'F'):
    temp_unit = temp_unit.upper()
    temp_unit_name = 'Fahrenheit'
else:
    temp_unit = temp_unit.upper()
    temp_unit_name = 'Kelvin'

# Precipitation unit
prcp_unit = input('\nChoose unit for precipitation amount (mm | in): ')

while(prcp_unit != 'mm' and prcp_unit != 'in'):
    prcp_unit = input('\nINVALID INPUT! Choose unit for precipitation amount (mm | in): ')

if(prcp_unit == 'mm'):
    prcp_unit_name = 'millimeters'
else:
    prcp_unit_name = 'inches'

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

# User input

temp = []
prcp = []

temp_file_directory = input('\nEnter the absolute path to the TXT file containing daily temperature data:\n')
prcp_file_directory = input('\nEnter the absolute path to the TXT file containing daily rainfall amount data:\n')

with open(temp_file_directory, 'r') as temp_file:
        for temp_line in temp_file.readlines():
            temp.append(float(temp_line.strip()))

with open(prcp_file_directory, 'r') as prcp_file:
        for prcp_line in prcp_file.readlines():
            prcp.append(float(prcp_line.strip()))

if(len(temp) != len(prcp)):
    print('\nWARNING: Input files do not have the same number of data values!')

print('\n----------------\n\nINPUT DATA for', month_name, year)
print('(temperature in ' + temp_unit_name + ', precipitation amount in ' + prcp_unit_name + ')')

count = int(0)

for x in range(len(temp)):
    count = count + 1
    print(year, month_name, count, '     Temp:', temp[x], '     Rainfall:', prcp[x])

input('\nPress ENTER to proceed to calculation.')

# Calculation proper and output display

count = 0
temp_sum = float(0)
prcp_sum = float(0)

for temp_val in temp:
    count = count + 1
    temp_sum = temp_sum + temp_val

temp_ave = temp_sum / count

for prcp_val in prcp:
    prcp_sum = prcp_sum + prcp_val

print('\n----------------\n\nOUTPUT for', month_name, year)
print('Monthly mean temperature: ' + str(round(temp_ave,1)) + temp_unit + '\nMonthly accumulated precipitation: ' + str(round(prcp_sum,1)) + ' ' + prcp_unit)

# End of program