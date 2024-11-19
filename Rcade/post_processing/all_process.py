#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 10:20:06 2021

@author: weiweijin
"""

# %% import lib & def functions
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import itertools

def chip_standardise(df):
    for i in range(1,len(df.columns)):
        df.iloc[:,i] =  (df.iloc[:,i]-df.iloc[:,i].min())/df.iloc[:,i].max()
    return(df)

def de_standardise(df):
    for i in range(len(df)):
        # if i != 0.0:
        df_std = (df.iloc[i,0]-df.iloc[:,0].min())/(df.iloc[:,0].max()-df.iloc[:,0].min())
        df.iloc[i,0] = df_std * 2 -1
    return(df)

def savefig(df,label,path):
    fig, ax = plt.subplots(figsize=(8,8))
    
    plt.rc('xtick',labelsize=16)
    plt.rc('ytick',labelsize=16) 
    
    sc = ax.scatter(df['logfc.DE'], df[label[8:]+'_M.ChIP'], c = df['B.ChIP.DE'], s =200)
    
    ax.set_xlabel('DE Log Ratio', fontsize=16)
    ax.set_ylabel(label[8:]+ ' Log Ratio', fontsize=16)
    
    plt.xlim(-18,18)
    plt.xticks(np.arange(-10, 12, 10)) 
    plt.ylim(0.5,6)
    plt.yticks(np.arange(1, 6, 2))
    
    sc.set_clim(vmin=1.5,vmax=4)
    cb = plt.colorbar(sc, ticks=[2,3,4])
    cb.ax.set_title('B', fontsize=16)
    cb.ax.set_yticklabels(['2', '3', '> 4'])  # vertically oriented colorbar
    
    f_path = os.path.join(path,'de_chip_'+label+'.png')
    fig.savefig(f_path)
    return(fig)

def heatmap_clustering(df,label,path):
    
    cmap = plt.get_cmap('summer')
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    
    lables = np.sign(df.iloc[:,0])
    # pal = sns.color_palette("seismic",len(lables.unique()))
    pal = sns.color_palette("seismic",2)
    # lut = dict(zip(lables.unique(), pal))
    lut = dict(zip((1,-1), pal))
    
    # lut = dict(zip(set(df.iloc[:,0]), sns.color_palette("bwr",len(set(df.iloc[:,0])))))
    sns.set(font_scale=2)
    fig = sns.clustermap(df.iloc[:,1:], metric = "euclidean", 
                         vmin = 0, vmax = 1, cmap=cmap,
                         #cbar = True,
                          # z_score=1, 
                         # cbar_kws={'shrink': 0.6},
                         row_colors=lables.map(lut),
                         row_cluster=True, col_cluster=False, yticklabels = False)
    plt.setp(fig.ax_heatmap.get_xticklabels(), rotation=75)
    # for lab in lables.unique():
    for lab in [1,-1]:
        if lab == 1:
            direction = 'Up'
        else:
            direction = 'Down'
        fig.ax_col_dendrogram.bar(0, 0, color=lut[lab],
                            label=direction, linewidth=0)
        fig.ax_col_dendrogram.legend(loc="center", ncol=2)
    
    f_path = os.path.join(path,'clustering_'+label+'.png')
    fig.savefig(f_path)
    return(fig)

def find_dup_val_in_dic(dict_org):
    flipped = {}
    
    for key, value in dict_org.items():
        if value not in flipped:
            flipped[value] = [key]
        else:
            flipped[value].append(key)
    
    return flipped

# %% input setting
input_gic_dir = 'output_GIC'
input_nsc_dir = 'output_NSC'
output_dir = 'all_process'
sub_id = ['ALL']

target_file_name = 'DEandChIP.csv'# file name of interest (the file include both de and chip)
th_B = 2.0 # threshhold for B (log-odds)
th_M = 1.0 # threshhold for M (log ratio for ChIP data)
antibody_names = ['H3K4me3', 'H3K36me3', 'H3K27me3', 'H3K27ac']

# %% import data for each patient and filter out data 
# GIC data
for sub in sub_id:
    sub_dir = os.path.join(input_gic_dir, sub)
    globals()[sub] = {}
    for i in range(4):
        f_path = os.path.join(sub_dir, str(i+1))
        temp = pd.read_csv(os.path.join(f_path, 'DEandChIP.csv'))
        temp = temp[temp['B.ChIP.DE']>th_B]
        temp = temp[temp['M.ChIP']>th_M]
        temp['names.DE'] = temp['names.DE'].str.split('.').str[0]
        temp = temp.drop_duplicates(subset=['names.DE'])
        temp = temp.rename(columns={"M.ChIP":antibody_names[i]+"_M.ChIP"})
        globals()[sub]['GIC_'+antibody_names[i]] = temp

# iNSC data
for sub in sub_id:
    sub_dir = os.path.join(input_nsc_dir, sub)
    for i in range(4):
        f_path = os.path.join(sub_dir, str(i+1))
        temp = pd.read_csv(os.path.join(f_path, 'DEandChIP.csv'))
        temp = temp[temp['B.ChIP.DE']>th_B]
        temp = temp[temp['M.ChIP']>th_M]
        temp['names.DE'] = temp['names.DE'].str.split('.').str[0]
        temp = temp.drop_duplicates(subset=['names.DE'])
        temp = temp.rename(columns={"M.ChIP":antibody_names[i]+"_M.ChIP"})
        globals()[sub]['NSC_'+antibody_names[i]] = temp
        
# # %% DE vs ChIP log ratio for individual patient and antibody
# # GIC data
# for sub in sub_id:
#     for i in range(4):
#         label = sub + '_GIC_' + antibody_names[i]
#         savefig(globals()[sub]['GIC_'+antibody_names[i]], label, output_dir)        

# # NSC data
# for sub in sub_id:
#     for i in range(4):
#         label = sub + '_NSC_' + antibody_names[i]
#         savefig(globals()[sub]['NSC_'+antibody_names[i]], label, output_dir) 

# %% Heatmap and clustering 
# prepare data
col_name = ['logfc.DE', 'GIC_H3K4me3', 'NSC_H3K4me3', 'GIC_H3K27ac', 'NSC_H3K27ac',
            'GIC_H3K27me3', 'NSC_H3K27me3','GIC_H3K36me3', 'NSC_H3K36me3'] # for re-arrange the col names
col_name2 = ['Gene Symbols', 'logfc.DE', 'GIC_H3K4me3', 'NSC_H3K4me3', 'GIC_H3K27ac', 'NSC_H3K27ac',
            'GIC_H3K27me3', 'NSC_H3K27me3','GIC_H3K36me3', 'NSC_H3K36me3'] # for re-arrange the col names
for sub in sub_id:
    for i in range(4):
        if i == 0:
            temp = globals()[sub]['GIC_'+antibody_names[i]][['logfc.DE', 'names.DE', antibody_names[i]+'_M.ChIP']]
            temp = temp.rename(columns={antibody_names[i]+'_M.ChIP':"GIC_"+antibody_names[i]})
        else:
            temp2 = globals()[sub]['GIC_'+antibody_names[i]][['logfc.DE', 'names.DE', antibody_names[i]+'_M.ChIP']]
            temp = pd.merge(temp,temp2,how='outer',on=['logfc.DE', 'names.DE'])
            temp = temp.rename(columns={antibody_names[i]+'_M.ChIP':"GIC_"+antibody_names[i]})
    for i in range(4):
        temp2 = globals()[sub]['NSC_'+antibody_names[i]][['logfc.DE', 'names.DE', antibody_names[i]+'_M.ChIP']]
        temp = pd.merge(temp,temp2,how='outer',on=['logfc.DE', 'names.DE'])
        temp = temp.rename(columns={antibody_names[i]+'_M.ChIP':"NSC_"+antibody_names[i]})
    temp.replace([np.inf, -np.inf], np.nan, inplace=True)
    temp = temp.fillna(0)
    temp = temp.rename(columns={'names.DE':'Gene Symbols'})
    # temp = temp.set_index('Gene Symbols')
    temp = temp[col_name2]
    temp.to_csv(os.path.join(output_dir,sub+'.csv'),index=False)
    temp = temp.drop(['Gene Symbols'], axis = 1)
    temp = de_standardise(temp)
    temp = chip_standardise(temp)
    temp = temp[col_name]
    temp = temp.fillna(0)
    temp = temp.rename(columns={'logfc.DE':'DE direction'})
    heatmap_clustering(temp,sub,output_dir)

# %% Histone Marks shared between GIC and iNSC
for sub in sub_id:
    for anti in antibody_names:
        temp = pd.merge(globals()[sub]['GIC_'+anti],globals()[sub]['NSC_'+anti],how='inner', on=['names.DE', 'logfc.DE'])
        temp2 = temp[['names.DE','logfc.DE', anti+'_M.ChIP_x', anti+'_M.ChIP_y']]
        globals()[sub][anti+'_shared'] = temp2.rename(columns={'names.DE':'Gene Symbols', anti+'_M.ChIP_x': anti+'_GIC_M.ChIP',
                                                               anti+'_M.ChIP_y': anti+'_NSC_M.ChIP'})
        
# %% heatmap for shared marks between GIC and iNSC
key_names = ['H3K4me3_shared', 'H3K27ac_shared', 'H3K27me3_shared', 'H3K36me3_shared']
for sub in sub_id:
    merged = globals()[sub][key_names[0]]
    for key in key_names[1:]:
        temp = globals()[sub][key]
        merged = pd.merge(merged, temp, how = 'outer', on = ['Gene Symbols','logfc.DE'])
    merged.to_csv(os.path.join(output_dir, 'shared',sub+'.csv'), index=False)
    merged = merged.drop('Gene Symbols', axis=1)
    merged.replace([np.inf, -np.inf], np.nan, inplace=True)
    merged = merged.fillna(0)
    merged = de_standardise(merged)
    merged = chip_standardise(merged)
    merged = merged.fillna(0)
    label = sub + '_shared'
    heatmap_clustering(merged,label,os.path.join(output_dir, 'shared'))
    
# %% shared marks among different patients
for anti in antibody_names:
    key = anti + '_shared'
    shared_stack = globals()[sub_id[0]][key].rename(columns={anti+'_GIC_M.ChIP':sub_id[0]+'_GIC', 
                                                             anti+'_NSC_M.ChIP':sub_id[0]+'_NSC'})
    for sub in sub_id[1:]:
        temp = globals()[sub][key].rename(columns={anti+'_GIC_M.ChIP':sub+'_GIC', anti+'_NSC_M.ChIP':sub+'_NSC'})
        shared_stack = pd.merge(shared_stack, temp, how = 'outer', on = ['Gene Symbols','logfc.DE'])
    shared_stack = shared_stack.drop('Gene Symbols', axis=1)
    shared_stack.replace([np.inf, -np.inf], np.nan, inplace=True)
    shared_stack = shared_stack.fillna(0)
    shared_stack = de_standardise(shared_stack)
    shared_stack = chip_standardise(shared_stack)
    shared_stack = shared_stack.fillna(0)
    shared_stack = shared_stack.rename(columns={'logfc.DE':'DE direction'})
    label_shared = anti+ '_shared_patient_shared'
    heatmap_clustering(shared_stack,label_shared,os.path.join(output_dir, 'shared'))

# %% find GIC and iNSC specific marks
for sub in sub_id:
    for anti in antibody_names:
        gic_marks = globals()[sub]['GIC_'+anti].geneID
        nsc_marks = globals()[sub]['NSC_'+anti].geneID
        gic_spc = list(set(gic_marks)-set(nsc_marks))
        nsc_spc = list(set(nsc_marks)-set(gic_marks))
        gic_data = globals()[sub]['GIC_'+anti].set_index('geneID')
        gic_spc_data = gic_data.loc[gic_spc]
        gic_spc_data = gic_spc_data[['names.DE','logfc.DE', anti+'_M.ChIP']]
        globals()[sub][anti+'_gic_speci'] = gic_spc_data.rename(columns={'names.DE':'Gene Symbols'})
        nsc_data = globals()[sub]['NSC_'+anti].set_index('geneID')
        nsc_spc_data = nsc_data.loc[nsc_spc]
        nsc_spc_data = nsc_spc_data[['names.DE','logfc.DE', anti+'_M.ChIP']]
        globals()[sub][anti+'_nsc_speci'] = nsc_spc_data.rename(columns={'names.DE':'Gene Symbols'})

# %% heatmap for GIC and iNSC specific marks
key_name_gic = ['H3K4me3_gic_speci', 'H3K27ac_gic_speci', 'H3K36me3_gic_speci', 'H3K27me3_gic_speci']
key_name_nsc = ['H3K4me3_nsc_speci', 'H3K27ac_nsc_speci', 'H3K36me3_nsc_speci', 'H3K27me3_nsc_speci']
gic_up_ratio = []
nsc_down_ratio = []
for sub in sub_id:
    gic_data = globals()[sub][key_name_gic[0]]
    for key in key_name_gic[1:]:
        temp = globals()[sub][key]
        gic_data = pd.merge(gic_data, temp, how = 'outer', on = ['Gene Symbols','logfc.DE'])
    label_gic = sub + '_gic_sp'
    gic_data.replace([np.inf, -np.inf], np.nan, inplace=True)
    gic_data = gic_data.fillna(0)
    gic_data.to_csv(os.path.join(output_dir,'gic_sp',label_gic+'.csv'), index=False)
    gic_data = gic_data.drop('Gene Symbols', axis=1)
    gic_data = de_standardise(gic_data)
    gic_data = chip_standardise(gic_data)
    gic_data = gic_data.fillna(0)
    gic_data = gic_data.rename(columns={'logfc.DE':'DE direction'})
    heatmap_clustering(gic_data,label_gic,os.path.join(output_dir,'gic_sp'))
    gic_up_ratio.append(round(len(gic_data.loc[gic_data.iloc[:,0]>0])/len(gic_data),3)*100)
    nsc_data = globals()[sub][key_name_nsc[0]]
    for key in key_name_nsc[1:]:
        temp = globals()[sub][key]
        nsc_data = pd.merge(nsc_data, temp, how = 'outer', on = ['Gene Symbols','logfc.DE'])
    label_nsc = sub + '_nsc_sp'
    nsc_data.replace([np.inf, -np.inf], np.nan, inplace=True)
    nsc_data = nsc_data.fillna(0)
    nsc_data.to_csv(os.path.join(output_dir,'nsc_sp',label_nsc+'.csv'), index=False)
    nsc_data = nsc_data.drop('Gene Symbols', axis=1)
    nsc_data = de_standardise(nsc_data)
    nsc_data = chip_standardise(nsc_data)
    nsc_data = nsc_data.fillna(0)
    nsc_data = nsc_data.rename(columns={'logfc.DE':'DE direction'})
    heatmap_clustering(nsc_data,label_nsc,os.path.join(output_dir,'nsc_sp'))
    nsc_down_ratio.append(round(len(nsc_data.loc[nsc_data.iloc[:,0]<0])/len(nsc_data),3)*100)
gic_up_df = pd.DataFrame({'Patients':sub_id,'Percentage of genes upregulated in GICs': gic_up_ratio})
gic_up_df.to_csv(os.path.join(output_dir,'gic_sp','GIC_up_ratio.csv'),index=False)
nsc_down_df = pd.DataFrame({'Patients':sub_id,'Percentage of genes downregulated in NSCs': nsc_down_ratio})
nsc_down_df.to_csv(os.path.join(output_dir,'nsc_sp','NSC_down_ratio.csv'),index=False)

# %% shared features in GIC/iNSC specifc data among patients
gic_shared_up = []
# subgroup = {"RTK I": ["P18", "P19", "P30", "P31"],
#             "RTK II": ["P17", "P50", "P54", "P61"],
#             "MES": ["P26", "P52"]}
for anti in antibody_names:
    key = anti + '_gic_speci'
    gic_stack = globals()[sub_id[0]][key].rename(columns={anti+'_M.ChIP':sub_id[0]})
    gic_stack['DE direction'] = np.sign(gic_stack['logfc.DE'])
    gic_stack = gic_stack.drop('logfc.DE', axis=1)
    gic_stack = gic_stack[["Gene Symbols","DE direction",sub_id[0]]]
    for sub in sub_id[1:]:
        temp = globals()[sub][key].rename(columns={anti+'_M.ChIP':sub})
        temp['DE direction'] = np.sign(temp['logfc.DE'])
        temp = temp.drop('logfc.DE', axis=1)
        gic_stack = pd.merge(gic_stack, temp, how = 'outer', on = ['Gene Symbols','DE direction'])
    label_gic = anti+ '_gic_patient_shared'
    gic_stack.replace([np.inf, -np.inf], np.nan, inplace=True)
    gic_stack = gic_stack.fillna(0)
    gic_stack1= gic_stack.copy()
    gic_stack.iloc[:,1]=['Up' if x > 0 else 'Down' for x in gic_stack.iloc[:,1]]
    gic_stack.to_csv(os.path.join(output_dir,'gic_sp', label_gic+'.csv'), index=False)
    # get shared genes between patients
    globals()[key+'_data'] = {}
    globals()[key+'_genes'] = {}
    for L in range(2,10):
        for subset in itertools.combinations(sub_id, L):
            col_names = list(subset)
            col_names.insert(0,'Gene Symbols')
            col_names.insert(1,'DE direction')
            sub_data = gic_stack[col_names]
            sub_data = sub_data[(sub_data.T != 0).all()].reset_index(drop = True)
            if len(sub_data) > 0:
                f_name = key 
                for ll in range(L):
                    f_name = f_name + '_' + subset[ll]
                globals()[key+'_data'][f_name] = sub_data
                globals()[key+'_genes'][f_name] = tuple(sub_data.iloc[:,0])
                # sub_data.to_csv(os.path.join(output_dir, f_name+'.csv'), index=False)
    # found subjects with same genes in the dic
    temp = find_dup_val_in_dic(globals()[key+'_genes'])
    f_name_list = []
    for k, v in temp.items():
        if len(v) == 1:
            f_name_list.append(v[0])
        else:
            logest = max(v, key=len)
            f_name_list.append(logest)
    globals()[key+'_final'] = {k2: globals()[key+'_data'][k2] for k2 in f_name_list}
    for k, v in globals()[key+'_final'].items():
        fn = k + '.csv'
        # v.iloc[:,1]=['Up' if x > 0 else 'Down' for x in v.iloc[:,1]]
        v.to_csv(os.path.join(output_dir,'gic_sp', key, fn), index=False)
    # # get shared genes in subgroups
    # for gn, pat in subgroup.items():
    #     c_name = pat.copy()
    #     c_name.insert(0,'Gene Symbols')
    #     c_name.insert(1,'DE direction')
    #     pat_data = gic_stack[c_name]
    #     pat_data = pat_data[(pat_data.T != 0).all()].reset_index(drop = True)
    #     # pat_data.iloc[:,1] = ['Up' if x > 0 else 'Down' for x in pat_data.iloc[:,1]]
    #     save_path = os.path.join(output_dir, 'gic_sp', key, gn, gn+'.csv')
    #     pat_data.to_csv(save_path,index=False)
    gic_stack1 = gic_stack1.drop('Gene Symbols', axis=1)
    # gic_stack = de_standardise(gic_stack)
    gic_stack1 = chip_standardise(gic_stack1)
    gic_stack1 = gic_stack1.fillna(0)
    # gic_stack = gic_stack.rename(columns={'logfc.DE':'DE direction'})
    heatmap_clustering(gic_stack1,label_gic,os.path.join(output_dir,'gic_sp'))
    gic_shared_up.append(round(len(gic_stack1.loc[gic_stack1.iloc[:,0]>0])/len(gic_stack1),3)*100)
gic_shared_up_df = pd.DataFrame({'Antibody':antibody_names,'Percentage of genes upregulated in GICs': gic_shared_up})
gic_shared_up_df.to_csv(os.path.join(output_dir,'gic_sp','GIC_up_ratio_shared.csv'),index=False)    

nsc_shared_down = []    
for anti in antibody_names:
    key = anti + '_nsc_speci'
    nsc_stack = globals()[sub_id[0]][key].rename(columns={anti+'_M.ChIP':sub_id[0]})
    nsc_stack['DE direction'] = np.sign(nsc_stack['logfc.DE'])
    nsc_stack = nsc_stack.drop('logfc.DE', axis=1)
    nsc_stack = nsc_stack[["Gene Symbols","DE direction",sub_id[0]]]
    for sub in sub_id[1:]:
        temp = globals()[sub][key].rename(columns={anti+'_M.ChIP':sub})
        temp['DE direction'] = np.sign(temp['logfc.DE'])
        temp = temp.drop('logfc.DE', axis=1)
        nsc_stack = pd.merge(nsc_stack, temp, how = 'outer', on = ['Gene Symbols','DE direction'])
    label_nsc = anti+ '_nsc_patient_shared'
    nsc_stack.replace([np.inf, -np.inf], np.nan, inplace=True)
    nsc_stack = nsc_stack.fillna(0)
    nsc_stack1 = nsc_stack.copy()
    nsc_stack.iloc[:,1]=['Up' if x > 0 else 'Down' for x in nsc_stack.iloc[:,1]]
    nsc_stack.to_csv(os.path.join(output_dir,'nsc_sp', label_nsc+'.csv'), index=False)
    nsc_stack1 = nsc_stack1.drop('Gene Symbols', axis=1)
    # nsc_stack = de_standardise(nsc_stack)
    nsc_stack1 = chip_standardise(nsc_stack1)
    nsc_stack1 = nsc_stack1.fillna(0)
    # nsc_stack = nsc_stack.rename(columns={'logfc.DE':'DE direction'})
    heatmap_clustering(nsc_stack1,label_nsc,os.path.join(output_dir,'nsc_sp'))
    nsc_shared_down.append(round(len(nsc_stack1.loc[nsc_stack1.iloc[:,0]<0])/len(nsc_stack1),3)*100)
nsc_shared_down_df = pd.DataFrame({'Antibody':antibody_names,'Percentage of genes downregulated in NSCs': nsc_shared_down})
nsc_shared_down_df.to_csv(os.path.join(output_dir,'nsc_sp','NSC_down_ratio_shared.csv'),index=False)

# # %% find genes shared between K4me3, K36me3 and K27ac but not in K27me3 or only in K27me3 
# for sub in sub_id:
#     gic_data = globals()[sub][key_name_gic[0]]
#     for key in key_name_gic[1:-1]:
#         temp = globals()[sub][key]
#         gic_data = pd.merge(gic_data, temp, how = 'inner', on = ['Gene Symbols','logfc.DE'])
#     # get the names of genes that shared between the 3 histone markers
#     shared_genes = set(gic_data.iloc[:,0])
#     # get the gene names for k27me3
#     k27me3_genes = set(globals()[sub][key_name_gic[3]].iloc[:,0])
#     # get the gene names only in the 3 histone markers 
#     shared_act = list(shared_genes-k27me3_genes)
#     gic_data = gic_data.set_index('Gene Symbols')
#     shared_data = gic_data.loc[shared_act]
#     shared_data.to_csv(os.path.join(output_dir, sub+'_shared_area.csv'),index=True)
#     # get the genes only in the K27me3
#     k27me3_only = list(k27me3_genes-shared_genes)
#     k27me3 = globals()[sub][key_name_gic[3]].set_index('Gene Symbols')
#     k27me3_data = k27me3.loc[k27me3_only]
#     k27me3_data.to_csv(os.path.join(output_dir, sub+'_k27me3.csv'),index=True)

# %% find genes in different histone mark combination:
'''
Potentical enhancer: K27ac - K4me3
Active promoter: K27ac + K4me3
Active gene bodies: K27ac + K36me3
Bivalent promoter: K27me3 + K4me3
Potential polycombs: K27me3 - K4me3
Potential polycombs: K27me3 - K36me3
Switch: K27me3 iNSC -> K27ac GIC
'''
key_add = [('H3K27ac_gic_speci','H3K4me3_gic_speci'), ('H3K27ac_gic_speci','H3K36me3_gic_speci'), ('H3K27me3_gic_speci','H3K4me3_gic_speci')]
key_subtract = [('H3K27ac_gic_speci','H3K4me3_gic_speci'), ('H3K27me3_gic_speci','H3K4me3_gic_speci'),('H3K27me3_gic_speci','H3K36me3_gic_speci')]

for sub in sub_id:
    gic_data = globals()[sub][key_name_gic[0]]
    for key in key_add:
        pd.merge(globals()[sub][key[0]],globals()[sub][key[1]], how = 'inner', on = ['Gene Symbols','logfc.DE']).to_csv(os.path.join(output_dir, sub+'_'+key[0][:-10]+'+'+key[1][:-10]+'_gic_speci.csv'),index=False)
    for key in key_subtract:
        gene_diff = list(set(globals()[sub][key[0]].iloc[:,0]) - set(globals()[sub][key[0]].iloc[:,1]))
        temp = globals()[sub][key[0]]
        temp = temp.set_index('Gene Symbols')
        temp.loc[gene_diff].to_csv(os.path.join(output_dir, sub+'_'+key[0][:-10]+'-'+key[1][:-10]+'_gic_speci.csv'),index=True)
    # for switch 
    insc_k27me3 = globals()[sub]['H3K27me3_nsc_speci']
    gic_k27ac = globals()[sub]['H3K27ac_gic_speci']
    pd.merge(insc_k27me3,gic_k27ac, how = 'inner', on = ['Gene Symbols','logfc.DE']).to_csv(os.path.join(output_dir, sub+'_K27me3_to_K27ac.csv'),index=False)

key_add = [('H3K27ac_nsc_speci','H3K4me3_nsc_speci'), ('H3K27ac_nsc_speci','H3K36me3_nsc_speci'), ('H3K27me3_nsc_speci','H3K4me3_nsc_speci')]
key_subtract = [('H3K27ac_nsc_speci','H3K4me3_nsc_speci'), ('H3K27me3_nsc_speci','H3K4me3_nsc_speci'),('H3K27me3_nsc_speci','H3K36me3_nsc_speci')]

for sub in sub_id:
    nsc_data = globals()[sub][key_name_nsc[0]]
    for key in key_add:
        pd.merge(globals()[sub][key[0]],globals()[sub][key[1]], how = 'inner', on = ['Gene Symbols','logfc.DE']).to_csv(os.path.join(output_dir, sub+'_'+key[0][:-10]+'+'+key[1][:-10]+'_nsc_speci.csv'),index=False)
    for key in key_subtract:
        gene_diff = list(set(globals()[sub][key[0]].iloc[:,0]) - set(globals()[sub][key[0]].iloc[:,1]))
        temp = globals()[sub][key[0]]
        temp = temp.set_index('Gene Symbols')
        temp.loc[gene_diff].to_csv(os.path.join(output_dir, sub+'_'+key[0][:-10]+'-'+key[1][:-10]+'_nsc_speci.csv'),index=True)
    