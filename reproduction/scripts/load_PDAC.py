# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 11:37:25 2021

@author: Gulben AVSAR

Functions to load datasets (GSE111672) as dataframes and anndata objects
"""
import pandas as pd
import numpy as np
import anndata
import warnings
import PIL.Image

path = './datasets/'
def load_seq_df (dname):
    if dname == 'PDAC-A':
        data = pd.read_csv(path + dname + '/' + 'GSE111672_PDAC-A-indrop-filtered-expMat.txt', sep='\t')
    elif dname == 'PDAC-B':
        data = pd.read_csv(path + dname + '/' + 'GSE111672_PDAC-B-indrop-filtered-expMat.txt', sep='\t')
    genes = data['Genes']
    data = data.drop(labels='Genes', axis=1)
    data.index = genes
    
    genes = [i for i in data.index]
    genes.sort()
    
    data = data.drop(genes[0:27], axis=0)
    return(data)

def load_ST_df (dname):
    if dname == 'PDAC-A':
        data = pd.read_csv(path + dname + '/' + 'GSM3036911_PDAC-A-ST1-filtered.txt', sep='\t')
    elif dname == 'PDAC-B':
        data = pd.read_csv(path + dname + '/' + 'GSM3405534_PDAC-B-ST1-filtered.txt', sep='\t')
    data.index=data['Genes']
    data = data.drop(['Genes'], axis=1)
    
    genes = [i for i in data.index]
    genes.sort()
    
    data = data.drop(genes[0:27], axis=0)
    return(data)



def load_ST_adata(dname):
    #file = [file for file in os.listdir(path) if file.endswith(".txt")]
    print('The dataset is loading')
    if dname == 'PDAC-A':
        data = pd.read_csv(path + dname + '/' + 'GSM3036911_PDAC-A-ST1-filtered.txt', sep='\t')
    elif dname == 'PDAC-B':
        data = pd.read_csv(path + dname + '/' + 'GSM3405534_PDAC-B-ST1-filtered.txt', sep='\t')
    
    data.index=data['Genes']
    data = data.drop(['Genes'], axis=1)
    
    genes = [i for i in data.index]
    genes.sort()
    
    data = data.drop(genes[0:27], axis=0)
    del genes
    print('There are {0} variables and {1} observations'.format(data.shape[0], data.shape[1]))
    
    X = data.T
    obs = data.columns.to_frame()
    obs.columns = ['spots']
    var = data.index.to_frame()
    adata = anndata.AnnData(X, obs, var)
    adata.obs['batch'] = np.zeros((len(obs),), dtype=int)
    
    return(adata)

def load_SEQ_adata(dname):
    warnings.filterwarnings("ignore")
    print('The dataset is loading')
    if dname == 'PDAC-A':
        data = pd.read_csv(path + dname + '/' + 'GSE111672_PDAC-A-indrop-filtered-expMat.txt', sep='\t')
    elif dname == 'PDAC-B':
        data = pd.read_csv(path + dname + '/' + 'GSE111672_PDAC-B-indrop-filtered-expMat.txt', sep='\t')
    
    genes = data['Genes']
    data = data.drop(labels='Genes', axis=1)
    data.index = genes
    del genes
    
    genes = [i for i in data.index]
    genes.sort()
    
    data = data.drop(genes[0:27], axis=0)
    del genes
    print('There are {0} variables and {1} observations'.format(data.shape[0], data.shape[1]))
    
    X = data.T
    X = X.reset_index(drop=True)
    obs = data.columns.to_frame()
    obs.columns = ['cell_type']
    obs = obs.reset_index(drop=True)
    var = data.index.to_frame()
    adata = anndata.AnnData(X, obs, var)
    return(adata)    


################## Load ST dataset and image as AnnData object (for stLearn) ##################
def STadata_imaged(dname, quality):
    
    print('The dataset is loading')
    if dname == 'PDAC-A':
        data = pd.read_csv(path + dname + '/' + 'GSM3036911_PDAC-A-ST1-filtered.txt', sep='\t')
    elif dname == 'PDAC-B':
        data = pd.read_csv(path + dname + '/' + 'GSM3405534_PDAC-B-ST1-filtered.txt', sep='\t')
    
    data.index=data['Genes']
    data = data.drop(['Genes'], axis=1)
    
    genes = [i for i in data.index]
    genes.sort()
    
    data = data.drop(genes[0:27], axis=0) #(19710,428)
    del genes
    print('There are {0} variables and {1} observations'.format(data.shape[0], data.shape[1]))
    
    X = data.T
    obs = data.columns.to_frame()
    obs.columns = ['spots']
    var = data.index.to_frame()
    adata = anndata.AnnData(X, obs, var)
    adata.obs['batch'] = np.zeros((len(obs),), dtype=int)
    #del adata.obs['spots']
    adata.obs['sum_counts'] = np.array(adata.X.sum(axis=1))
    
    adata.uns["spatial"] = dict()
    library_id = dname
    adata.uns["spatial"][library_id] = dict()
    

    adata.uns["spatial"][library_id]['images'] = dict()
    
    rgba_img = PIL.Image.open(path + dname + '/' + '/spatial/tissue_fullres_image.png')
    img = rgba_img.convert('RGB')
    adata.uns["spatial"][library_id]['images'][quality] = np.array(img)

    adata.uns["spatial"][library_id]["use_quality"] = quality    

    positions = pd.read_csv(path + dname + '/' + '/spatial/tissue_positions_list.csv', header=None)
    positions.columns = [ 'spot_id', 'in_tissue', 'array_col', 'array_row', 'x_pxl', 'y_pxl']
    positions.index = positions['spot_id']
    adata.obs = adata.obs.join(positions, how="left")
    adata.obsm['spatial'] = adata.obs[['x_pxl', 'y_pxl']].to_numpy()
    
    adata.obs.drop(columns=['spot_id', 'x_pxl', 'y_pxl'], inplace=True)
    
    image_coor = adata.obsm["spatial"]
    adata.obs["imagecol"] = image_coor[:, 0]
    adata.obs["imagerow"] = image_coor[:, 1]
    adata.uns["spatial"]["use_quality"] = quality
    
    return(adata)


################## Load imputed ST dataset and image as AnnData object (for evaluation) ##################
def impSTadata_imaged(dname, imp_data, list_obs_names, quality):    
    
    imp_data.index = list_obs_names
    X = imp_data
    obs = imp_data.index.to_frame()
    obs.columns = ['spots']
    var = imp_data.columns.to_frame()
    adata = anndata.AnnData(X, obs, var)
    adata.obs['batch'] = np.zeros((len(obs),), dtype=int)
    adata.var['Genes'] = imp_data.columns.to_frame()
    del adata.var[0]
    
    adata.obs['batch'] = np.zeros((len(obs),), dtype=int)
    adata.obs['sum_counts'] = np.array(adata.X.sum(axis=1))
       
    adata.uns["spatial"] = dict()
    library_id = dname
    adata.uns["spatial"][library_id] = dict()
    

    adata.uns["spatial"][library_id]['images'] = dict()
    
    rgba_img = PIL.Image.open(path + dname + '/' + '/spatial/tissue_fullres_image.png')
    img = rgba_img.convert('RGB')
    adata.uns["spatial"][library_id]['images'][quality] = np.array(img)

    adata.uns["spatial"][library_id]["use_quality"] = quality    

    positions = pd.read_csv(path + dname + '/' + '/spatial/tissue_positions_list.csv', header=None)
    positions.columns = [ 'spot_id', 'in_tissue', 'array_col', 'array_row', 'x_pxl', 'y_pxl']
    positions.index = positions['spot_id']
    adata.obs = adata.obs.join(positions, how="left")
    adata.obsm['spatial'] = adata.obs[['x_pxl', 'y_pxl']].to_numpy()
    
    adata.obs.drop(columns=['spot_id', 'x_pxl', 'y_pxl'], inplace=True)
    
    image_coor = adata.obsm["spatial"]
    adata.obs["imagecol"] = image_coor[:, 0]
    adata.obs["imagerow"] = image_coor[:, 1]
    
    adata.uns["spatial"][library_id]['scalefactors'] = dict()
    adata.uns["spatial"][library_id]['scalefactors']['spot_diameter_fullres'] = 7.0/2
    adata.uns["spatial"][library_id]['scalefactors']['tissue_{}_scalef'.format(quality)] = 1.0

    return(adata)











