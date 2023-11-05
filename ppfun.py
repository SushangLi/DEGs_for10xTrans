# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 18:29:47 2021

@author: LiuJiangyuan
"""

## Import dependencies
import os
import gzip
import scanpy as sc
import pandas as pd
import numpy as np

## Functions
# unzip .gz files
def ungz(gz_path):
 	for f in  os.listdir(gz_path):
         if ".gz" in f:
             g = gzip.GzipFile(mode="rb", fileobj=open(gz_path+"\\"+f, 'rb'))
             open(gz_path+"\\"+f.replace(".gz",""), "wb").write(g.read())

# split the tissue regions in selected spots to different files
def split_region(sample_list, general_label_name, region_name):
    for sample in sample_list:
        # read in the user-defined barcode file
        selected_labels = pd.read_csv(sample + '/' + general_label_name + '.csv', sep=',')
        for region in region_name:
            selected_labels[selected_labels.iloc[:,1] == region].to_csv(sample + '/' + region + ".csv", sep=',', index=False, header=True)


# compute the total count of all selected spots for each gene probe, 
# return the summarized table for each sample and an integrated table for the input list.
def summary_data(sample_list, sample_list_name, use_label, general_label_name, output_name):
    colnames = ['gene_name'] + sample_list_name
    summary = pd.DataFrame(columns=colnames)

    # summarize all samples in the sample_list list
    index = 0
    for sample in sample_list:
        print('Current sample: ' + sample)
        # read in the raw data of the current sample
        features = pd.read_csv(sample + '/features.tsv', sep='\t', names=['id','name','type'])
        barcodes = pd.read_csv(sample + '/barcodes.tsv', sep='\t', names=['barcodes'])
        matrix = sc.read(sample + '/matrix.mtx').X.todense()

        # extract some spots from the raw data or not 
        if use_label == True:
            # read in the user-defined barcode file
            selected_labels = pd.read_csv(sample+'/'+ general_label_name +'.csv', sep=',')
    
            # extract the barcodes as variables
            sample_barcodes = barcodes.barcodes
            selected_barcodes = selected_labels.Barcode
    
            # extract barcode matched spots to a new dataframe
            match_df = pd.DataFrame(index = range(0, len(matrix)))
            for i in range(0, len(sample_barcodes)):
                for j in range(0, len(selected_labels)):
                    if selected_barcodes[j] == sample_barcodes[i]:
                        tmpdf = pd.DataFrame(matrix[:,i], columns = [sample_barcodes[i]])
                        match_df = match_df.join(tmpdf)
        else:
            sample_barcodes = barcodes.barcodes
            match_df = pd.DataFrame(matrix, columns = [sample_barcodes])

        # summarize the (selected) spots
        geneName = features.name   
        geneExp = match_df.sum(axis=1)     # total count of all selected spots for each gene probe
        summary['gene_name'] = geneName
        summary[sample_list_name[index]] = geneExp
        
        # export up to current sample summary table to the sample directory
        summary.to_csv(sample + '/' + output_name + '.csv', sep=',', index=False, header=True)

        # move to the next sample id
        index = index + 1
    
    # export the integrated summary table to current working directory
    summary.to_csv('./' + output_name + '.csv', sep=',', index=False, header=True)
    print('Summary Finished!')

# append frequency & add same gene data for sum_table
def sum_duplicated(summary_file_name, sample_list_name):
    summary = pd.read_csv(summary_file_name + '.csv')
    colnames = ['gene_name'] + sample_list_name + ['frequency']
    sample_size = len(sample_list_name)
    sum_table = pd.DataFrame(columns = colnames)    
    gene_list = list()
    
    unique_index = 0     # record the unique number of genes
    for _, probe in summary.iterrows():
        # mark the first appearing genes and extract the same rows from summary data to the sum_table
        if not probe['gene_name'] in gene_list:
            gene_list.append(probe['gene_name'])
            sum_table.loc[unique_index] = probe
            sum_table.at[unique_index, 'frequency'] = 1
            unique_index += 1
        # add the counts to existed rows for genes have multiple probes, recording the number of probes
        else:
            former_index = gene_list.index(probe['gene_name'])     # extract the unique index of this probe in former recored gene list
            for sample_index in range(1, sample_size + 1):
                sum_table.iloc[former_index, sample_index] += probe[sample_index]
            sum_table.at[former_index, 'frequency'] += 1
    return sum_table

# get the mean value of gene counts for multiple probes and remove the column 'frequency'
def mean_duplicated(sum_table, sample_list_name, output_name):
    sample_size = len(sample_list_name)
    
    gene_name = sum_table.iloc[:, 0]
    frequency_vec = sum_table.iloc[:, -1]
    sum_table_count =  sum_table.iloc[:, 1:sample_size + 1]
    
    mean_table_count = sum_table_count.div(frequency_vec, axis = 'rows')
    
    mean_table = pd.concat([gene_name, mean_table_count], axis = 1)

    mean_table.to_csv('./' + output_name + '.csv', sep=',', index=False, header=True)
    print(output_name + ' Finished!')
    
    return mean_table

