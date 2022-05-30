#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 17:43:04 2022

Functions to load synthetic ST dataset as dataframe and anndata object
"""
import pandas as pd
import anndata
import scanpy as sc

# path = './benchmark/datasets/'
path = './datasets/'
# dname = 'SimulatedST'

def load_seq_df (dname):
    adata_file = '5705STDY8058280_filtered_feature_bc_matrix.h5'
    adata = sc.read_10x_h5(path + dname + '/' + adata_file)
    adata.var_names_make_unique()
    adata.var_names = [i.upper() for i in adata.var_names]
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    adata = adata[adata.obs.n_genes_by_counts < 6000, :] #(8804,31053)
    sc.pp.filter_genes(adata, min_cells=1) #(8804,24734)
    
    data = adata.to_df()
    return(data)


def load_seq_adata (dname):
    adata_file = '5705STDY8058280_filtered_feature_bc_matrix.h5'
    adata = sc.read_10x_h5(path + dname + '/' + adata_file)
    adata.var_names_make_unique()
    adata.var_names = [i.upper() for i in adata.var_names]
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    adata = adata[adata.obs.n_genes_by_counts < 6000, :] #(8804,31053)
    sc.pp.filter_genes(adata, min_cells=1) #(8804,24734)
    return(adata)


def load_st_df(dname):
    data = pd.read_csv(path + dname + '/' + 'synthetic_ST_seed467_1_counts.csv', sep=',', index_col=0)
    data.columns = [i.upper() for i in data.columns]
    data = data[data.sum(axis=1) > 0] #(717,31053)
    data = data.loc[:,data.sum(axis=0) > 0] #(717,21470)
    return(data)
    

def load_st_adata(dname):
    data = pd.read_csv(path + dname + '/' + 'synthetic_ST_seed467_1_counts.csv', sep=',', index_col=0)
    data.columns = [i.upper() for i in data.columns]
    data = data[data.sum(axis=1) > 0]
    data = data.loc[:,data.sum(axis=0) > 0]
    
    X = data #genes (vars) should be at the columns
    obs = data.index.to_frame()
    obs.columns = ['spots']
    var = data.columns.to_frame()
    adata = anndata.AnnData(X, obs, var)    
    return(adata)











