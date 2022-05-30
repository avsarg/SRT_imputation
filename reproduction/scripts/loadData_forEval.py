# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 18:01:49 2022

@author: Gulben AVSAR

Function to upload preprocessed STR datasets for the evaluation

"""
from scripts import load_PDAC as ld
from scripts import load_Liver as ldLiv
from scripts import load_SimulatedST as ldSst
import numpy as np
import scanpy as sc
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

def getGenes():
    df = pd.read_csv("genes_for_imputation.txt", sep='\t')
    df = df.set_index('Genes')
    df = df.sort_index(axis = 0)
    return(df)

def load_datasets (dName, mName):
    # dName: name of dataset: 'PDAC-A', 'PDAC-B', or 'Liver'
    # mName: name of method: 'spage', 'stplus', 'gimvi', 'tangram','5stlearn' or '30stlearn'
    imp_paths = './predictions/'
    geneSet = getGenes()
    if mName == 'spage' or mName == 'stplus':
        if dName == 'PDAC-A' or dName == 'PDAC-B':
            seq = ld.load_seq_df(dName)
            seq = seq[seq.sum(axis=1) > 0]
            
            # filtering
            Genes_count = np.sum(seq > 0, axis=1)
            seq = seq.loc[Genes_count >=10,:]
            
            # Normalization of seq data
            def Log_Norm_cpm(x):
                return np.log(((x/np.sum(x))*1000000) + 1)
            seq = seq.apply(Log_Norm_cpm,axis=0)
            
            # upload ST dataset
            st = ld.load_ST_df(dName)
            
            #subsetting
            genes_to_keep = [st.index[i] for i in range(0,len(st)) if np.any(seq.index == st.index[i])]
            st = st.loc[genes_to_keep, :]
            
            # Normalization for ST datasets
            cell_count = np.sum(st, axis=0)
            def Log_Norm_spatial(x):
                return np.log(((x/np.sum(x))*np.median(cell_count)) + 1)
            st = st.apply(Log_Norm_spatial,axis=0)
        
            st = st.T
            org_st = st.copy()
            
            dName = dName.split('-')[0].lower() + dName.split('-')[1]
            imp_st = pd.read_csv(imp_paths + '{}/'.format(mName) +
                                 'Imp_{}_{}.txt.gz'.format(mName,dName)
                                 ,index_col=0)
            midx ='nan'
        elif dName == 'Liver':
            SEQadata = ldLiv.load_SEQ_adata()
            #filter genes with no expression
            sc.pp.filter_genes(SEQadata,min_cells=10)
        
            # convert to dataframe
            seq = SEQadata.to_df()
            
            ####### Load the Liver ST dataset as dataframes
            STadata = ldLiv.load_ST_adata(dName)
        
            #filter genes with no expression
            sc.pp.filter_genes(STadata,min_cells=1)
        
            # convert to dataframe
            st = STadata.to_df()
        
            # the genes of ST data should be a subset of scRNA-seq data in SpaGE
            # so, filter the genes that are not found in scRNA-seq data
            shrd = np.intersect1d(seq.columns,st.columns)
            st = st.loc[:,shrd]
            org_st = st.copy()
            
            imp_st = pd.read_csv(imp_paths + '{}/'.format(mName) +
                                 'Imp_{}_{}.txt.gz'.format(mName,dName)
                                 ,index_col=0)
        
            # Normalization for ST datasets
            cell_count = np.sum(st, axis=1)
            def Log_Norm_spatial(x):
                return np.log(((x/np.sum(x))*np.median(cell_count)) + 1)
            st = st.apply(Log_Norm_spatial,axis=1)
            org_st = st.copy()
            midx ='nan'
        elif dName == 'SimulatedST':
            seq = ldSst.load_seq_df(dName) #(8804,24773)
            st = ldSst.load_st_df(dName) #(717,21470)
            
            seq = seq.T
            st = st.T
            
            seq = seq[seq.sum(axis=1) > 0]
            st = st[st.sum(axis=1) > 0]
            
            Genes_count = np.sum(seq > 0, axis=1)
            seq = seq.loc[Genes_count >=10,:]
            
            genes_to_keep = np.intersect1d(st.index,seq.index)
            st = st.loc[genes_to_keep, :] #(19669,717)
            seq = seq.loc[genes_to_keep, :] #(19669,8804)
            
            def Log_Norm_cpm(x):
                return np.log(((x/np.sum(x))*1000000) + 1)
            seq = seq.apply(Log_Norm_cpm,axis=0)
            
            
            cell_count = np.sum(st, axis=0) #herbir hücredeki toplam count değerini verir
            def Log_Norm_spatial(x):
                return np.log(((x/np.sum(x))*np.median(cell_count)) + 1)
            st = st.apply(Log_Norm_spatial,axis=0) #(12576,428)
            
            st = st.T
            org_st = st.copy()
            
            imp_st = pd.read_csv(imp_paths + '{}/'.format(mName) +
                                 'Imp_{}_{}.txt.gz'.format(mName,dName)
                                 ,index_col=0)
            midx ='nan'

    if mName == 'gimvi':
        if dName == 'PDAC-A' or dName == 'PDAC-B':
            dName = dName.split('-')[0].lower() + dName.split('-')[1]
        org_st = pd.read_csv(imp_paths + '{}/'.format(mName) + 'Orig_gimvi_{}.txt.gz'.format(dName))
        imp_st = pd.read_csv(imp_paths + '{}/'.format(mName) + 'Imp_gimvi_{}.txt.gz'.format(dName))
        midx ='nan'
    
    if mName == 'tangram':
        if dName == 'PDAC-A' or dName == 'PDAC-B':
            # Load the datasets as dataframes
            st = ld.load_ST_df(dName)
            
            #filter genes with no expression
            st = st[st.sum(axis=1) > 0]
            
            st = st.T
            
            dName = dName.split('-')[0].lower() + dName.split('-')[1]
            imp_st = pd.read_csv(imp_paths + '{}/'.format(mName) +
                                 'Imp_{}_{}.txt.gz'.format(mName,dName),
                                 index_col=None)
        elif dName == 'Liver':
            STadata = ldLiv.load_ST_adata(dName)

            #filter genes with no expression
            sc.pp.filter_genes(STadata,min_cells=1)
            
            # convert to dataframe
            st = STadata.to_df()
            
            
            imp_st = pd.read_csv(imp_paths + '{}/'.format(mName) +
                                 'Imp_{}_{}.txt.gz'.format(mName,dName),
                                 index_col=None)
        elif dName == 'SimulatedST':
            st = ldSst.load_st_adata(dName) #(717,21470)
            
            sc.pp.filter_genes(st,min_cells=1)
    
            # convert to dataframe
            st = st.to_df()
            imp_st = pd.read_csv(imp_paths + '{}/'.format(mName) +
                                 'Imp_{}_{}.txt.gz'.format(mName,dName),
                                 index_col=None)
        ###### Normalization
        def Log_Norm_cpm(x):
            return np.log((x/np.sum(x)) + 1)
        st = st.apply(Log_Norm_cpm,axis=1)
        
        org_st = st.copy()
        midx ='nan'
        
    
    if mName == '5stlearn' or mName == '30stlearn':
        if dName == 'Liver':
            st = ldLiv.load_ST_adata(dName)
        elif dName == 'SimulatedST':
            st = ldSst.load_st_adata(dName) #(717,21470)
        else:
            st = ld.STadata_imaged(dName, quality='fullres')
            dName = dName.split('-')[0].lower() + dName.split('-')[1]
        
        #filter genes with no expression
        sc.pp.filter_genes(st,min_cells=1)
        
        # normalize the count values
        sc.pp.normalize_total(st)
        
        # tranfrom the normalized count values
        sc.pp.log1p(st)
        
        # Converts anndata objects to dataframes
        st = st.to_df()
        
        org_st = st.copy()
        
        imp_st = pd.read_csv(imp_paths + '{}/'.format(mName) +
                             'Imp_{}_{}.txt.gz'.format(mName,dName)
                             ,index_col=None)
        
        midx = pd.read_csv(imp_paths+ '{}/'.format(mName) +
                           'MaskedIdx_{}_{}.txt.gz'.format(mName,dName),
                           header=None, index_col=0)

    org_st = org_st[geneSet.index]
    imp_st = imp_st[geneSet.index]
    
    return(org_st, imp_st, midx)


