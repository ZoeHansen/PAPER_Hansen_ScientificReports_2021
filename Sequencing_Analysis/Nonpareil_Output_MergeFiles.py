# -*- coding: utf-8 -*-
#################################################

# Merging Output from Nonpareil for Analysis

#################################################

### This code is designed to use nonpareil_output.npo files and merge them into a single file for plotting in R.

# Load necessary packages
import pandas as pd
import numpy as np
import os
import csv

#Set root directory
rootdir1=r'J://HPCC/SeqRun_1/'
rootdir2=r'J://HPCC/SeqRuns_234/'
os.chdir(rootdir1)

# Designate our sample numbers we wish to iterate through
samples = os.listdir(rootdir1)

# Make a table containing our file name, sample ID, and color designations

data = []

for i in samples:
    sample_dir = rootdir1 + ''.join(i) + '/nonpareil/'
    for file in os.listdir(sample_dir):
        if file == 'nonpareil_output.npo':
            dir_path = sample_dir+file
            print(dir_path)
            data.append((dir_path, i))
            
df = pd.DataFrame(data, columns = ['File', 'ID'])
print(df)


R_vector = np.random.randint(low = 0, high = 256, size = 79)
G_vector = np.random.randint(low = 0, high = 256, size = 79)
B_vector = np.random.randint(low = 0, high = 256, size = 79)

df['R'] = R_vector
df['G'] = G_vector
df['B'] = B_vector

print(df)


# Now, make a second file for the SeqRuns_234; after we construct this, we will join the two dataframes

rootdir2=r'J://HPCC/SeqRuns_234/'
os.chdir(rootdir2)

samples2 = os.listdir(rootdir2)

data2 = []

for j in samples2:
    sample_dir2 = rootdir2 + ''.join(j) + '/nonpareil/'
    for file2 in os.listdir(sample_dir2):
        if file2 == 'nonpareil_output.npo':
            dir_path2 = sample_dir2+file2
            print(dir_path2)            
            data2.append((dir_path2, j))
            
df2 = pd.DataFrame(data2, columns = ['File', 'ID'])
print(df2)

R_vector2 = np.random.randint(low = 0, high = 256, size = 197)
G_vector2 = np.random.randint(low = 0, high = 256, size = 197)
B_vector2 = np.random.randint(low = 0, high = 256, size = 197)

df2['R'] = R_vector2
df2['G'] = G_vector2
df2['B'] = B_vector2

print(df2)


# Append df2 to the first dataframe for a comprehensive dataframe. 

df_all = df.append(df2, ignore_index = True)
df_all['ID']=df_all['ID'].astype(int)
print(df_all)



# Now, we will take information from our filtered subset of samples to exclude those that we 
# are not interested in. 

filtered = pd.read_csv('I://Manning_ER1N/CampylobacterSubset_AIM_ONE/Second_Analysis/Data_files_Hansen_2020/campylobacter_metadata_Hansen_2020.csv', 
                       sep = ',', header = 0)
f = pd.DataFrame(filtered[['ID','Case.status']])

# Merge our df_all dataframe with the subset of the filtered output to include health status
# as well. 

final_df = pd.merge(f, df_all, on = 'ID', how = 'left')
print(final_df)

# Write to a .csv file for later implementation in R
final_df.to_csv('J://HPCC/Nonpareil_output_all_FINAL2.csv', sep = ',', index = False )
