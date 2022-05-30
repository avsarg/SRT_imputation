# -*- coding: utf-8 -*-
"""
Created on Tue May 24 14:48:54 2022

@author: Gulben AVSAR
"""
import pandas as pd
import matplotlib.pyplot as plt
import NaiveDE
import SpatialDE
import numpy as np

def prep_data(data, dName):
    
    counts = data.copy()
    counts = counts.T[counts.sum(0) >= 0].T  # Filter practically unobserved genes
        
    sample_info = pd.read_csv('./datasets/'+dName+'/spatial/tissue_positions_list.csv',
                              index_col=0, header=None)
    
    inTissue = sample_info[sample_info.columns[[0]]] == 1
    sample_info = sample_info[inTissue[inTissue.columns[0]]]
    
    # sample_info[sample_info.columns[[0,1,2]]]
    sample_info = sample_info.drop(sample_info.columns[[0,1,2]], axis=1)
    sample_info.columns = ['x','y']
    
    tot_count = counts.sum(1)
    
    shrd = np.intersect1d(sample_info.index, tot_count.index)
    
    tot_count = tot_count[shrd]
    sample_info = sample_info.loc[shrd,:]
    
    sample_info['total_counts'] = tot_count[sample_info.index]
    
    # sample_info.head(5)
    counts = counts.loc[sample_info.index]  # Align count matrix with metadata table
    
    return(counts, sample_info)


def SpatialDE_run(data, dName):
    
    counts, sample_info = prep_data(data, dName)
    
    # the dataset were already normalized
    # norm_expr = NaiveDE.stabilize(counts.T).T
    resid_expr = NaiveDE.regress_out(sample_info, counts.T, 'np.log(total_counts)').T
    
    X = sample_info[['x', 'y']]
    results = SpatialDE.run(X, resid_expr)
    
    return(results)


def pull_signf_DEGs(df, metric):
    res_sorted = df.sort_values(metric)[['g', 'l','pval', 'qval']]
    res_sign = res_sorted[res_sorted[metric] < 0.05]
    
    return(res_sign)


def visualise_genes(data, dName, geneList,s_size, fig_size):
    # visualise 3 geens
    counts, sample_info = prep_data(data, dName)
    
    g = geneList
    fig = plt.subplots(1,3,figsize=fig_size)
    for i, g in enumerate(g):
        plt.subplot(1, 3, i + 1)
        plt.scatter(sample_info['x'], sample_info['y'], c=counts[g], s=s_size);
        plt.title(g)
        # plt.axis('equal')
        # plt.ylim([200, 650])
        # plt.gca().invert_yaxis()    
        plt.colorbar(ticks=[])
        plt.tight_layout()
    
    return(fig)




