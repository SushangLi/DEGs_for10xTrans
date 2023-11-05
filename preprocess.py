# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 14:55:49 2023

@author: LiuJiangyuan
"""
## Import dependecies
import os
import gzip
import scanpy as sc
import pandas as pd
from ppfun import ungz, split_region, summary_data, sum_duplicated, mean_duplicated

## Set the working directory
os.chdir("../data")

## Set the paths of raw barcodes, features, matrix files
# here are 4 groups of mice, 2 of young vs 2 of old, each group has 4 mice as biological repeats. 
g1_1 = '/work/path/WTY2-A_result/outs/filtered_feature_bc_matrix'
g1_2 = '/work/path/WTY2-B_result/outs/filtered_feature_bc_matrix'
g1_3 = '/work/path/WTY2-C_result/outs/filtered_feature_bc_matrix'
g1_4 = '/work/path/WTY2-D_result/outs/filtered_feature_bc_matrix'

g2_1 = '/work/path/P14-A_result/outs/filtered_feature_bc_matrix'
g2_2 = '/work/path/P14-B_result/outs/filtered_feature_bc_matrix'
g2_3 = '/work/path/P14-C_result/outs/filtered_feature_bc_matrix'
g2_4 = '/work/path/P14-D_result/outs/filtered_feature_bc_matrix'

g3_1 = '/work/path/LLS20201222A_result/outs/filtered_feature_bc_matrix'
g3_2 = '/work/path/LLS20201222B_result/outs/filtered_feature_bc_matrix'
g3_3 = '/work/path/LLS20201222C_result/outs/filtered_feature_bc_matrix'
g3_4 = '/work/path/LLS20201222D_result/outs/filtered_feature_bc_matrix'

g4_1 = '/work/path/total_Slice-MEC/data/LLS-1/A'
g4_2 = '/work/path/total_Slice-MEC/data/LLS-1/B'
g4_3 = '/work/path/total_Slice-MEC/data/LLS-1/C'
g4_4 = '/work/path/total_Slice-MEC/data/LLS-1/D'

g5_1 = '/work/path/total_Slice-MEC/data/LLS-2/A'
g5_2 = '/work/path/total_Slice-MEC/data/LLS-2/B'
g5_3 = '/work/path/total_Slice-MEC/data/LLS-2/C'
g5_4 = '/work/path/total_Slice-MEC/data/LLS-2/D'



# concat and name the list of groups
sample_list = [g1_1, g1_2, g1_3, g1_4, 
	       g2_1, g2_2, g2_3, g2_4, 
	       g3_1, g3_2, g3_3, g3_4, 
	       g4_1, g4_2, g4_3, g4_4, 
	       g5_1, g5_2, g5_3, g5_4]
sample_list_name = ['g1_1', 'g1_2', 'g1_3', 'g1_4', 
                   	  'g2_1', 'g2_2', 'g2_3', 'g2_4', 
                  	  'g3_1', 'g3_2', 'g3_3', 'g3_4', 
                   	  'g4_1', 'g4_2', 'g4_3', 'g4_4', 
                   	  'g5_1', 'g5_2', 'g5_3', 'g5_4']

# unzip files in the group list
if 1:
    for path in sample_list:
        ungz(path)
    print('Ungz Finished!')
    
## Get the total counts of selected spots (with the individual or mean value for each probe) and then export
# split the MEC labels to seperate files
region_name = ['layer1','layer2','layer3','deep_layer']
split_region(sample_list, general_label_name = 'MEClabel', region_name = region_name)

# summary the regions (MEC layer, MEC total, slice total) in each sample and all the samples
# the MEC layer
for region in region_name:
    summary_data(sample_list, sample_list_name, use_label = True, general_label_name = region, output_name = 'MEC_' + region + '_summary')
# the MEC total
summary_data(sample_list, sample_list_name, use_label = True, general_label_name = 'MEClabel', output_name = 'MEC_total_summary')
# the slice total
summary_data(sample_list, sample_list_name, use_label = False, general_label_name = '', output_name = 'Slice_total_summary')


# compute the sum and mean of the summary data (which locates at the current working directory)
# the MEC layer
for region in region_name:
    sum_table_MEC_layer = sum_duplicated(summary_file_name = 'MEC_' + region + '_summary', sample_list_name = sample_list_name)
    mean_table_MEC_layer = mean_duplicated(sum_table = sum_table_MEC_layer, sample_list_name = sample_list_name, output_name = 'MEC_' + region + '_mean')

# the MEC total 
sum_table_MEC = sum_duplicated(summary_file_name = 'MEC_total_summary', sample_list_name = sample_list_name)
mean_table_MEC = mean_duplicated(sum_table = sum_table_MEC, sample_list_name = sample_list_name, output_name = 'MEC_total_mean')

# the slice total
sum_table_slice = sum_duplicated(summary_file_name = 'Slice_total_summary', sample_list_name = sample_list_name)
mean_table_slice = mean_duplicated(sum_table = sum_table_slice, sample_list_name = sample_list_name, output_name = 'Slice_total_mean')










