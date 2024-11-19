#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 12:41:50 2021

@author: weiweijin
"""

# %% import lib & def functions
import pandas as pd
from glob import glob
from pathlib import Path

# %% import data
input_dir = 'input'
output_dir = 'output'

input_files = glob(input_dir+ '/*.csv')

# %% find genes are close to enhancer regions (5kb)
th = 5000
for ff in input_files:
    temp = pd.read_csv(ff)
    temp[temp['shortestDistance']<th].to_csv(Path(output_dir)/(ff[6:-4]+"filt.csv"))