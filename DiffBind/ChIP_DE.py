#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 14:43:13 2021

@author: weiweijin
"""

# %% import lib & def functions
import pandas as pd
import os

def find_csv_filenames( path_to_dir, suffix=".csv" ):
    filenames = os.listdir(path_to_dir)
    return [ filename for filename in filenames if filename.endswith( suffix ) ]


# %% Import data
input_de = 'input_DE'
input_chip = 'input_ChIP_RTKI'
output_dir = 'output_RTKI'

# de = pd.read_excel(os.path.join(input_de,'full_de_syngeneic_only.xlsx'),index_col=0).dropna(subset=['Gene Symbol']).drop(['consistent'], axis=1)

de = pd.read_csv(os.path.join(input_de,'RNA_RTKI.csv')).dropna(subset=['Names']).rename(columns={"Names":"Gene Symbol"})


# Get the current working directory
cwd = os.getcwd()
all_csv_files = find_csv_filenames(os.path.join(cwd,input_chip), suffix=".csv" )

# %% data processing 
# list of combinations should have positive logFC
pos_list = ['K27ac-K4me3_gain','K27ac+K4me3_gain','K27ac+K36me3_gain','K27me3-K27ac_switch_gain','K27me3-K4me3_loss','K27me3-K36me3_loss','K27me3+K4me3_loss']
neg_list = ['K27ac-K4me3_loss','K27ac+K4me3_loss','K27ac+K36me3_loss','K27me3-K4me3_gain','K27me3-K36me3_gain','K27me3+K4me3_gain']

for ff in all_csv_files:
    tmp = pd.read_csv(os.path.join(cwd,input_chip,ff), error_bad_lines=False).rename(columns={"geneName":"Gene Symbol"})
    if tmp['Gene Symbol'].isnull().values.all() == False:
        tmp2 = pd.merge(tmp,de,how="inner",on="Gene Symbol").drop_duplicates(subset=['Gene Symbol']).set_index(['Gene Symbol'])
        tmp2.to_csv(os.path.join(output_dir, ff[:-16]+'with_DE.csv'))
        # if it is in the pos_list, only keep the logFC are positive
        if ff[5:-17] in pos_list:
            tmp2[tmp2.logFC>0].to_csv(os.path.join(output_dir, ff[:-16]+'with_DE_filt.csv'))
            # if it is in the neg_list, only keep the logFC are negative
        else:
        # elif ff[:-17] in neg_list:
            tmp2[tmp2.logFC<0].to_csv(os.path.join(output_dir, ff[:-16]+'with_DE_filt.csv'))