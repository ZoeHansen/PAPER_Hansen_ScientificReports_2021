# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 14:20:54 2019

@author: hansenzo
"""

#Convert the .tsv files to .csv files (try to remove columns, rename columns in the process)
import pandas as pd
import os
import numpy as np

rootdir = r'C://Users/Zoe/Manning_ER1N/Campylobacter_Resistomes/'
os.chdir(rootdir)
    
samples = os.listdir(rootdir)

for i in samples:
    sample_dir = r'I://Manning_ER1N/Campylobacter_Resistomes/' + "".join(i)
    os.chdir(sample_dir)
    sample_files = os.listdir(sample_dir)
    
    for file in sample_files:
        if file.startswith('gene'):
            
            gene_file = pd.read_csv(file, sep = "\t", header = None)
            gene_files = pd.DataFrame(gene_file)
            gene_files.columns = gene_files.iloc[0]
            gene_files.drop(gene_files.index[0], inplace = True)
            gene_files.drop("Gene Fraction", axis = 1, inplace = True)
            gene_files.drop("Sample", axis = 1, inplace = True)
            gene_files.drop(gene_files.columns[[0]], axis = 1, inplace = True)
            gene_files.replace({'Hits':'ERIN_'+"".join(i)}, regex = True, inplace = True)
            gene_files.to_csv('gene_' + ''.join(i) + '.csv', sep = ',', index = True)
            
        print(gene_files)    
        
        
#Merge the resistome gene data (can use this code as a template for the mechanism, group, class data)
# 'ccZID' refers to cases and controls that were included in the final analysis of this manuscript 

ccZID = ['27','29','52','63','68','73','86','88','110','117','121','133','138','141','146','150','162','165','171','174',
        '179','182','183','268','271','273','37','38','39','40','42','44','65','67','70','71','72','75','76','78','79',
        '80','81','197','198','199','200','224','225','226','228','229','230','250','257','259','260','269','270','277',
        '278','279','280','281','282']

# Create a "starter" for the file merge
merged = pd.DataFrame()
c = ['Gene', 'ERIN_27','ERIN_29']
subZ = ccZID[2:]

os.chdir(rootdir + ''.join(ccZID[0]))
main = pd.read_csv('gene27.csv', sep = ',', header = 0) #names = ["Gene", "ERIN_27"])

os.chdir(rootdir + ''.join(ccZID[1]))
main2 = pd.read_csv('gene29.csv', sep = ',', header = 0) #names = ["Gene", "ERIN_29"])

merged = pd.merge(main, main2, on = 'Gene', how = 'outer')


#Merge the rest of the samples:

for i in subZ:
    os.chdir(rootdir + ''.join(i))
    main3 = pd.read_csv('gene' + ''.join(i) + '.csv', sep = ',', header = 0)
    
    merged = pd.merge(merged, main3, on = 'Gene', how = 'outer')
    merged.replace(np.nan, 0, inplace = True) #Replace NaN with zeroes
    
print(merged)

os.chdir("I://Manning_ER1N/Campylobacter_Resistome_Data/")
merged.to_csv('GENE_merged_clean_CC.csv', sep = ',', index = False)


#### Repeat this script for the Group, Mechanism, and Class level output for all samples







