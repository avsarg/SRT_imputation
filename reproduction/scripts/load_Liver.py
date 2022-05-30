# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 17:04:53 2022

@author: Gulben AVSAR

Functions to load datasets (GSE185477) as anndata objects
"""
import scanpy as sc
import pandas as pd
#import stlearn as st

################## Load datasets as dfs (for SpaGE) ##################

def load_SEQ_adata():
    print('The dataset is loading')
    seqpath = './datasets/Liver/C58_SC/'
    adata = sc.read_10x_mtx(
        seqpath,  # the directory with the `.mtx` file
        var_names='gene_symbols')
    adata.var_names_make_unique()
    
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.n_genes_by_counts < 1500, :]
    adata = adata[adata.obs.n_genes_by_counts > 200, :]
    adata = adata[adata.obs.pct_counts_mt < 80, :]
    return(adata)


################## Load datasets as dfs (for gimVI) ##################
def load_ST_adata(dname):
    print('The dataset is loading')
    stpath = './datasets/Liver/'
    adata = sc.read_visium(stpath, library_id=dname)
    adata.var_names_make_unique()
    #adata #(737280,33538)
    obsNames = adata.obs_names
    positions = pd.read_csv(stpath + '/spatial/tissue_positions_list.csv', header=None)
    positions.columns = [ 'spot_id', 'in_tissue', 'array_col', 'array_row', 'x_pxl', 'y_pxl']
    positions.index = positions['spot_id']
    adata.obs = adata.obs.merge(positions, how="left")
    adata.obsm['spatial'] = adata.obs[['x_pxl', 'y_pxl']].to_numpy()
    
    adata.obs.drop(columns=['spot_id', 'x_pxl', 'y_pxl'], inplace=True)
    
    image_coor = adata.obsm["spatial"]
    adata.obs["imagecol"] = image_coor[:, 0]
    adata.obs["imagerow"] = image_coor[:, 1]
    
    adata.obs_names = obsNames
    
    adata.uns["spatial"][dname]["use_quality"] = 'hires'
    
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.n_genes_by_counts < 6000, :]
    adata = adata[adata.obs.pct_counts_mt < 4, :]
    return(adata)


