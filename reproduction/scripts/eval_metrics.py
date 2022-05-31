# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 12:05:44 2022

@author: Gulben AVSAR
"""
import pandas as pd
from scripts import evals as evl

def getGenes():
    df = pd.read_csv("genes_for_imputation.txt", sep='\t')
    df = df.set_index('Genes')
    df = df.sort_index(axis = 0)
    return(df)

def eval_GW(org_data, imp_data, dName, mName, masked_idxs, eval_path):
    #gene-wise evaluation
    imp_data.index = org_data.index
    genes = getGenes()
    geneSet = list(genes.index)
    metrics = {}
    if mName in ['spage','stplus','gimvi','tangram']:
        pcc = evl.eval_PCC(org_data, imp_data, geneSet)
        cs = evl.eval_CosSim(org_data, imp_data, geneSet)
        rmsle = evl.eval_RMSLE(org_data, imp_data, geneSet)
    elif mName in ['5stlearn', '30stlearn']:
        pcc = evl.eval_maskedPCC(org_data, imp_data, geneSet, masked_idxs)
        cs = evl.eval_maskedCosSim(org_data, imp_data, geneSet, masked_idxs)
        rmsle = evl.eval_maskedRMSLE(org_data, imp_data, geneSet, masked_idxs)
    else:
        pcc = evl.eval_PCC(org_data, imp_data, geneSet)
        cs = evl.eval_CosSim(org_data, imp_data, geneSet)
        rmsle = evl.eval_RMSLE(org_data, imp_data, geneSet)
    metrics = {'pcc':pcc, 'cs':cs, 'rmsle':rmsle}
    
    pcc.to_csv(eval_path+'genewise_PCC_{}_{}.txt'.format(dName,mName))
    cs.to_csv(eval_path+'genewise_CS_{}_{}.txt'.format(dName,mName))
    rmsle.to_csv(eval_path+'genewise_RMSLE_{}_{}.txt'.format(dName,mName))
    return(metrics)

def eval_SW(org_data, imp_data, dName, mName, masked_idxs, eval_path):
    #spot-wise evaluation
    imp_data.index = org_data.index
    genes = getGenes()
    LSgenes = genes.Groups.copy()
    LSgenes = LSgenes[LSgenes>=3]
    geneSet = list(LSgenes.index)
    metrics = {}
    if mName in ['spage','stplus','gimvi','tangram']:
        pcc = evl.eval_spotPCC(org_data, imp_data, geneSet)
        cs = evl.eval_spotCS(org_data, imp_data, geneSet)
        rmsle = evl.eval_spotRMSLE(org_data, imp_data, geneSet)
    elif mName in ['5stlearn', '30stlearn']:
        pcc = evl.eval_maskedspotPCC(org_data, imp_data, geneSet, masked_idxs)
        cs = evl.eval_maskedspotCS(org_data, imp_data, geneSet, masked_idxs)
        rmsle = evl.eval_maskedspotRMSLE(org_data, imp_data, geneSet, masked_idxs)
    else:
        pcc = evl.eval_PCC(org_data, imp_data, geneSet)
        cs = evl.eval_CosSim(org_data, imp_data, geneSet)
        rmsle = evl.eval_RMSLE(org_data, imp_data, geneSet)
    metrics = {'pcc':pcc, 'cs':cs, 'rmsle':rmsle}
    
    pcc.to_csv(eval_path+'spotwise_PCC_{}_{}.txt'.format(dName,mName))
    cs.to_csv(eval_path+'spotwise_CS_{}_{}.txt'.format(dName,mName))
    rmsle.to_csv(eval_path+'spotwise_RMSLE_{}_{}.txt'.format(dName,mName))
    return(metrics)


def CS_pvals(org_data, imp_data, dName, mName, masked_idxs, eval_path):
    #gene-wise evaluation
    imp_data.index = org_data.index
    genes = getGenes()
    geneSet = list(genes.index)
    metrics = {}
    if mName in ['spage','stplus','gimvi','tangram']:
        CS = evl.CS_pval(org_data, imp_data, geneSet)
    elif mName in ['5stlearn', '30stlearn']:
        CS = evl.maskedCS_pval(org_data, imp_data, geneSet, masked_idxs)
    metrics = {'cs':CS}
    
    CS.to_csv(eval_path+'genewise_CSpvals_{}_{}.txt'.format(dName,mName))
    return(metrics)

def PCC_pvals(org_data, imp_data, dName, mName, masked_idxs, eval_path):
    #gene-wise evaluation
    imp_data.index = org_data.index
    genes = getGenes()
    geneSet = list(genes.index)
    metrics = {}
    if mName in ['spage','stplus','gimvi','tangram']:
        CS = evl.PCC_pval(org_data, imp_data, geneSet)
    elif mName in ['5stlearn', '30stlearn']:
        CS = evl.maskedPCC_pval(org_data, imp_data, geneSet, masked_idxs)
    metrics = {'cs':CS}
    
    CS.to_csv(eval_path+'genewise_PCCpvals_{}_{}.txt'.format(dName,mName))
    return(metrics)


def RHO_pvals(org_data, imp_data, dName, mName, masked_idxs, eval_path):
    #gene-wise evaluation
    imp_data.index = org_data.index
    genes = getGenes()
    geneSet = list(genes.index)
    metrics = {}
    if mName in ['spage','stplus','gimvi','tangram']:
        CS = evl.RHO_pval(org_data, imp_data, geneSet)
    elif mName in ['5stlearn', '30stlearn']:
        CS = evl.maskedRHO_pval(org_data, imp_data, geneSet, masked_idxs)
    metrics = {'cs':CS}
    
    CS.to_csv(eval_path+'genewise_RHOpvals_{}_{}.txt'.format(dName,mName))
    return(metrics)



