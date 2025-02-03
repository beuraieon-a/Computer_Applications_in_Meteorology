import os
from datetime import datetime, timedelta

# Define the start date and time
start_datetime = datetime(2025, 2, 3, 0, 0)
time_interval = timedelta(minutes=10)

# Get the list of files in the directory
directory = 'C:/Users/charl/Pictures/sat_img/true_color/WPAC/'
# directory = 'C:/Users/charl/Pictures/sat_img/true_color/PH/'
# directory = 'C:/Users/charl/Pictures/sat_img/eir/WPAC/'
# directory = 'C:/Users/charl/Pictures/sat_img/eir/PH/'
files = sorted([f for f in os.listdir(directory) if f.startswith('Image') and f.endswith(".png")])

# Loop through each file and rename it
for index, file in enumerate(files):
    # Calculate the new date and time
    new_datetime = start_datetime + index * time_interval
    
    new_filename = f"sat_truecolor_WPAC_{new_datetime.strftime('%Y%m%d_%H%M')}.png"
    # new_filename = f"sat_truecolor_PH_{new_datetime.strftime('%Y%m%d_%H%M')}.png"
    # new_filename = f"sat_EIR_WPAC_{new_datetime.strftime('%Y%m%d_%H%M')}.png"
    # new_filename = f"sat_EIR_PH_{new_datetime.strftime('%Y%m%d_%H%M')}.png"
    
    # Construct full file paths
    old_file_path = os.path.join(directory, file)
    new_file_path = os.path.join(directory, new_filename)
    
    # Rename the file
    os.rename(old_file_path, new_file_path)
    print(f"Renamed '{file}' to '{new_filename}'")

print("Renaming completed.")