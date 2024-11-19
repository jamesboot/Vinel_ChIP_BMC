#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 10:04:25 2021

@author: weiweijin
"""

# %% import lib & def functions
from pathlib import Path
import glob

# %% get average
def get_avg(l):
    total = 0
    for ll in l:
        total = total + float(ll)
    return total / len(l)

# %% import data
all_files = glob.glob('input/*.bed')

# %% save avg
for ff in all_files:
    with open(ff,"r") as fin:
        for _ in range(1):
            next(fin)
        fout = open(Path("output") / (ff[6:]+'Graph'), "w")
        for line in fin:
            line = line.strip().split(sep = "\t") # csv with ";" ...
            avg = get_avg(line[3:][:]) 
            new_line = line[:3]
            new_line.append(str(avg))
            fout.write('\t'.join(new_line) + '\n')
        fout.close()
