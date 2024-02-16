"""Post-imputational filtering: implementation of quality control measures.   

This script treats raw data from Impute2 and filters out imputed SNPs with MAF < 0.01 and info score < 0.2 
(thresholds may be adjusted).
"""

import os
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import glob

# load data
print("Read info data")
#paths_info = glob.glob("/Volumes/My Passport for Mac/GWAS/AllImputation_chr16/imputed/chr16_chunk_*_info")
paths_info = glob.glob("/Volumes/My Passport for Mac/GWAS/chr13/CHR13/chr13_chunk_*_info")

def treat_data_chunk(info_path):
    """Add MAF column."""
    pd_info_chunk = pd.read_csv(info_path, sep=" ")
    pd_info_chunk["exp_freq_a0"] = 1-pd_info_chunk["exp_freq_a1"]
    pd_info_chunk["MAF"] = pd_info_chunk[['exp_freq_a0','exp_freq_a1']].min(axis=1)

    return pd_info_chunk

def filter_info_MAF(info_path, info_limit):
    """Filter dataframe of imputation info by MAF (>0.01) and "info" scores."""
    info_chunk = treat_data_chunk(info_path)
    info_chunk_filtered_MAF = info_chunk[info_chunk["MAF"]>0.01]
    info_chunk_filtered = info_chunk_filtered_MAF[info_chunk_filtered_MAF["info"]>info_limit]
    info_chunk_clean = info_chunk_filtered[info_chunk_filtered["a0"].str.len()==1]
    info_chunk_clean = info_chunk_clean[info_chunk_clean["a1"].str.len()==1]
    return info_chunk_clean, len(info_chunk["MAF"])

def filter_MAF_zero(info_path):
    """Filter dataframe of imputation info, exclude MAF = 0."""
    info_chunk = treat_data_chunk(info_path)
    info_chunk_fltered_light = info_chunk[info_chunk["MAF"]>0]
    return info_chunk_fltered_light

def filter_MAF(info_path):
    """Filter dataframe of imputation info, exclude MAF < 0.01."""
    info_chunk = treat_data_chunk(info_path)
    info_chunk_fltered = info_chunk[info_chunk["MAF"]>0.1]
    return info_chunk_fltered

# filter all info files and store positions of filtered SNPs
filtered_positions = []
n_lines = 0
for info_path in paths_info:
    info_chunk_filtered, n = filter_info_MAF(info_path, info_limit=0.2)
    filtered_positions.extend(info_chunk_filtered["position"])
    n_lines += n
    
print("Filtered positions stored: ", len(filtered_positions), "positions")
print("Number of lines to treat: ", n_lines)

# Plot post-imputational quality metrics
info_list = []
maf_list = []

for info_path in paths_info:
    info_chunk = treat_data_chunk(info_path)
    info_list.extend(info_chunk["info"])
    maf_list.extend(info_chunk["MAF"])

# Info and MAF before filtering
fig = plt.figure(figsize=(8,5))
fig.patch.set_facecolor('white')
plt.title("Info scores distribution before filtering",fontsize=16)
plt.hist(info_list, bins=100)
plt.xlabel("Info",fontsize=14)
fig.savefig('figures/Info_raw.png')
plt.close(fig)

fig = plt.figure(figsize=(8,5))
plt.title("MAF distribution before filtering",fontsize=16)
plt.hist(maf_list, bins=100)
plt.xlabel("MAF",fontsize=14)
fig.savefig('figures/MAF_raw.jpg')
plt.close(fig)

print("Total number of imputations before filtering: ", len(info_list))
print("Number of imputations before filtering with info > 0.2: ", sum(np.array(info_list) > 0.2))

# Info after excluding MAF = 0
info_list_maf_light = []
maf_list_maf_light = []

info_list_maf_filtered = []
maf_list_maf_filtered = []

for info_path in paths_info:
    info_chunk_filtered_zero = filter_MAF_zero(info_path)
    info_chunk_filtered = filter_MAF(info_path)
    
    info_list_maf_light.extend(info_chunk_filtered_zero["info"])
    maf_list_maf_light.extend(info_chunk_filtered_zero["MAF"])
    
    info_list_maf_filtered.extend(info_chunk_filtered["info"])
    maf_list_maf_filtered.extend(info_chunk_filtered["MAF"])
    
fig = plt.figure(figsize=(8,5))
plt.title("Info scores distribution for MAF>0", fontsize=16)
plt.hist(info_list_maf_light, bins=200)
plt.xlabel("Info", fontsize=14)
fig.savefig('figures/Info_MAFpos.jpg')
plt.close(fig)

print("Total number of imputations with MAF > 0: ", len(info_list_maf_light))
print("Number of imputations with MAF > 0 and info > 0.2: ", sum(np.array(info_list_maf_light) > 0.2))

# Info and MAF after filtering
fig = plt.figure(figsize=(8,5))
plt.title("Info scores distribution for MAF>0.01", fontsize=16)
plt.hist(info_list_maf_filtered, bins=200)
plt.xlabel("Info", fontsize=14)
fig.savefig('figures/Info_MAFfilt.jpg')
plt.close(fig)

fig = plt.figure(figsize=(8,5))
plt.title("MAF distribution with MAF > 0.01",fontsize=16)
plt.hist(maf_list_maf_filtered, bins=100)
plt.xlabel("MAF",fontsize=14)
fig.savefig('figures/MAF_filt.jpg')
plt.close(fig)

print("Total number of imputations with MAF > 0.01: ", len(info_list_maf_filtered))
print("Number of imputations with  MAF > 0.01 and info > 0.2: ", sum(np.array(info_list_maf_filtered) > 0.2))

print("Figures stored")

######## filter chunk of haplotypes    
impute_chunk = "/Volumes/My Passport for Mac/GWAS/chr13/CHR13/chr13_concat"
#impute_chunk = "/Users/yuliya/Documents/Master2023/cours/GENOM/imputation/concatenated chunks/chr16_chunkAll3.impute2"

out_file = open("CHR13_filtered_new", 'a')

"""Iterate line-by-line over raw data of imputed SNPs and write in a new file only SNPs in filtered_positions."""
s = 0
# catch lines of non-standard format
position_problem = 0
# iterate over input file
with open(impute_chunk) as f:
    for line in f:
        s+=1
        line_elements = line.split(" ")
        if len(line_elements) > 2:
            position = line.split(" ")[2]
        else: 
            position_problem += 1
            print(line_elements[:50])
            continue
        if position.isdigit():
            SNP_position = int(position)
            if SNP_position in filtered_positions:
                out_file.write(line) 
        else:
           position_problem += 1
           print(line_elements[:50])
           print("Problem with BP position in line ", s) 
        if s%10000 == 0: print("Lines treated: ", s)
            
out_file.close()

print("N lines with wrong positions: ", position_problem)
print("All done!")