# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 12:05:03 2022

@author: Gulben AVSAR

Plotting and PCC
"""
import numpy as np
import pandas as pd

################################################################################################
################################ EVALUATION METRICS ############################################

def nz4genes(org_data, imp_data, geneSet):
    idx = {} #indices of non-zero values for given genes
    for i in geneSet:
        if i not in org_data.columns:
            print(i, '    not found in ST dataset')
            continue
        else:
            indices = [idx for idx, element in enumerate(org_data.loc[:,i]) if element != 0]
            idx[i] = indices
    return(idx)

def nz4genesMasked(org_data, imp_data, geneSet, masked_idxs):
    idx = {} #indices of non-zero values for given genes
    for i in geneSet:
        indices = masked_idxs.loc[i,:].dropna().tolist()
        idx[i] = [int(x) for x in indices]
    return(idx)

def nz4spots(org_data, imp_data, geneSet, spots):
    idx = {} #indices of non-zero values for given spots
    for i in spots:
        if i not in org_data.index:
            print(i, '    not found in ST dataset')
            continue
        else:
            indices = [idx for idx, element in enumerate(org_data.loc[i,geneSet]) if element != 0]
            idx[i] = indices
    return(idx)

################################################################################################
######################## Gene-based evals ################################################
######################## PCC (gw)
#Pearson's correlation coefficients (PCC) for each gene
#between the NON-ZERO original values and the predicted values
def eval_PCC(org_data, imp_data, geneSet):
    spots = np.array(org_data.index)
    corrs = pd.DataFrame(index = geneSet, columns=['corr_val'])
    idx = nz4genes(org_data, imp_data, geneSet)
    
    for i in range(0,len(geneSet)):
        c = org_data.loc[spots[idx[geneSet[i]]],geneSet[i]].corr(imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]])
        corrs.loc[geneSet[i],'corr_val'] = c    
    return(corrs)


#Pearson's correlation coefficients (PCC) for each gene for masked experiment
#between the original values and the predicted values
def eval_maskedPCC(org_data, imp_data, geneSet, masked_idxs):
    spots = np.array(org_data.index)
    corrs = pd.DataFrame(index = geneSet, columns=['corr_val'])
    idx = nz4genesMasked(org_data, imp_data, geneSet, masked_idxs)
    
    for i in range(0,len(geneSet)):
        c = org_data.loc[spots[idx[geneSet[i]]],geneSet[i]].corr(imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]])
        corrs.loc[geneSet[i],'corr_val'] = c
    return(corrs)


######################## RMSLE (gw)
#Root Mean Squared Log Error (RMSLE) for each gene
#between the NON-ZERO original values and the predicted values
def eval_RMSLE(org_data, imp_data, geneSet):
    spots = np.array(org_data.index)
    rmsles = pd.DataFrame(index = geneSet, columns=['rmsle_val'])
    idx = nz4genes(org_data, imp_data, geneSet)
    
    for i in range(0,len(geneSet)):
        org = org_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        imp = imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        rmsle = np.sqrt(np.square(np.log(imp + 1e-12) - np.log(org + 1e-12)).mean())
        rmsles.loc[geneSet[i],'rmsle_val'] = rmsle
    return(rmsles)


#Root Mean Squared Log Error (RMSLE) for each gene for masked experiment
#between the NON-ZERO original values and the predicted values
def eval_maskedRMSLE(org_data, imp_data, geneSet, masked_idxs):
    spots = np.array(org_data.index)
    rmsles = pd.DataFrame(index = geneSet, columns=['rmsle_val'])
    idx = nz4genesMasked(org_data, imp_data, geneSet, masked_idxs)
    
    for i in range(0,len(geneSet)):
        org = org_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        imp = imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        rmsle = np.sqrt(np.square(np.log(imp + 1e-12) - np.log(org + 1e-12)).mean())
        rmsles.loc[geneSet[i],'rmsle_val'] = rmsle
    return(rmsles)

######################## CosSim (gw)
#Cosine similarity for each gene
#between the NON-ZERO original values and the predicted values
def eval_CosSim(org_data, imp_data, geneSet):
    spots = np.array(org_data.index)
    cos_sims = pd.DataFrame(index = geneSet, columns=['cos_similarity'])
    idx = nz4genes(org_data, imp_data, geneSet)
        
    for i in range(len(geneSet)):
        org = org_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        imp = imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        norm_sq = np.linalg.norm(org) * np.linalg.norm(imp)
        cos_sims.loc[geneSet[i],'cos_similarity'] = (org @ imp) / norm_sq
    return(cos_sims)


#Cosine similarity for each gene for masked experiment
#between the NON-ZERO original values and the predicted values
def eval_maskedCosSim(org_data, imp_data, geneSet, masked_idxs):
    spots = np.array(org_data.index)
    cos_sims = pd.DataFrame(index = geneSet, columns=['cos_similarity'])
    idx = nz4genesMasked(org_data, imp_data, geneSet, masked_idxs)
        
    for i in range(len(geneSet)):
        org = org_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        imp = imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        norm_sq = np.linalg.norm(org) * np.linalg.norm(imp)
        cos_sims.loc[geneSet[i],'cos_similarity'] = (org @ imp) / norm_sq
    return(cos_sims)



######################## CosSim with p-values (gw)
import math, random
from scipy import stats
def CS_pval(org_data, imp_data, geneSet):
    spots = np.array(org_data.index)
    cos_sims = pd.DataFrame(index = geneSet, columns=['cos_similarity','pvals'])
    idx = nz4genes(org_data, imp_data, geneSet)
        
    for i in range(len(geneSet)):
        org = org_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        imp = imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        # print('{}  gene {}'.format(i,geneSet[i]))
        
        similarity = lambda x1, x2: sum(xj*xk for xj,xk in zip(x1, x2))/math.sqrt(sum(xj**2 for xj in x1)*sum(xk**2 for xk in x2))
        x1 = list(org.values)
        x2 = list(imp.values)
        s = similarity(x1, x2)
        
        if np.isnan(s):
            cos_sims.loc[geneSet[i],'cos_similarity'] = s
            cos_sims.loc[geneSet[i],'pvals'] = np.nan
        else:
            ## permutation test
            lx, sr = len(x1), []
            for j in range(10000):
                mj = random.sample(x1, lx)
                sr.append(similarity(mj, x1))
            shape, loc, scale = stats.weibull_min.fit(sr)
            ## -log10(p)
            ej = ((s-loc)/scale)**shape*math.log10(math.exp(1.))
            p = 10**(-ej)
    
            cos_sims.loc[geneSet[i],'cos_similarity'] = s
            cos_sims.loc[geneSet[i],'pvals'] = p
        
    return(cos_sims)


def maskedCS_pval(org_data, imp_data, geneSet, masked_idxs):
    spots = np.array(org_data.index)
    cos_sims = pd.DataFrame(index = geneSet, columns=['cos_similarity','pvals'])
    idx = nz4genesMasked(org_data, imp_data, geneSet, masked_idxs)
        
    for i in range(len(geneSet)):
        org = org_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        imp = imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        # print('{}  gene {}'.format(i,geneSet[i]))
        
        similarity = lambda x1, x2: sum(xj*xk for xj,xk in zip(x1, x2))/math.sqrt(sum(xj**2 for xj in x1)*sum(xk**2 for xk in x2))
        x1 = list(org.values)
        x2 = list(imp.values)
        s = similarity(x1, x2)
        if np.isnan(s):
            cos_sims.loc[geneSet[i],'cos_similarity'] = s
            cos_sims.loc[geneSet[i],'pvals'] = np.nan
        else:
            ## permutation test
            lx, sr = len(x1), []
            for j in range(10000):
                mj = random.sample(x1, lx)
                sr.append(similarity(mj, x2))
            shape, loc, scale = stats.weibull_min.fit(sr)
            ## -log10(p)
            ej = ((s-loc)/scale)**shape*math.log10(math.exp(1.))
            p = 10**(-ej)
    
            cos_sims.loc[geneSet[i],'cos_similarity'] = s
            cos_sims.loc[geneSet[i],'pvals'] = p
        
    return(cos_sims)

######################## PCC with p-values (gw)
#Pearson's correlation coefficients (PCC) for each gene
#between the NON-ZERO original values and the predicted values
def PCC_pval(org_data, imp_data, geneSet):
    spots = np.array(org_data.index)
    corrs = pd.DataFrame(index = geneSet, columns=['corr_val','pvals'])
    idx = nz4genes(org_data, imp_data, geneSet)
    
    for i in range(0,len(geneSet)):
        org = org_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        imp = imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        if len(org) < 2:
            c = [np.nan, np.nan]
        else:
            c = stats.pearsonr(org,imp)
        corrs.loc[geneSet[i],'corr_val'] = c[0]
        corrs.loc[geneSet[i],'pvals'] = c[1]
    return(corrs)


#Pearson's correlation coefficients (PCC) for each gene for masked experiment
#between the original values and the predicted values
def maskedPCC_pval(org_data, imp_data, geneSet, masked_idxs):
    spots = np.array(org_data.index)
    corrs = pd.DataFrame(index = geneSet, columns=['corr_val','pvals'])
    idx = nz4genesMasked(org_data, imp_data, geneSet, masked_idxs)
    
    for i in range(0,len(geneSet)):
        org = org_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        imp = imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        if len(org) < 2:
            c = [np.nan, np.nan]
        else:
            c = stats.pearsonr(org,imp)
        corrs.loc[geneSet[i],'corr_val'] = c[0]
        corrs.loc[geneSet[i],'pvals'] = c[1]
    return(corrs)

######################## Spearman`s rho with p-values (gw)
#Pearson's correlation coefficients (PCC) for each gene
#between the NON-ZERO original values and the predicted values
def RHO_pval(org_data, imp_data, geneSet):
    spots = np.array(org_data.index)
    corrs = pd.DataFrame(index = geneSet, columns=['rho','pvals'])
    idx = nz4genes(org_data, imp_data, geneSet)
    
    for i in range(0,len(geneSet)):
        org = org_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        imp = imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        if len(org) < 2:
            c = [np.nan, np.nan]
        else:
            c = stats.spearmanr(org,imp)
        corrs.loc[geneSet[i],'rho'] = c[0]
        corrs.loc[geneSet[i],'pvals'] = c[1]
    return(corrs)


#Pearson's correlation coefficients (PCC) for each gene for masked experiment
#between the original values and the predicted values
def maskedRHO_pval(org_data, imp_data, geneSet, masked_idxs):
    spots = np.array(org_data.index)
    corrs = pd.DataFrame(index = geneSet, columns=['rho','pvals'])
    idx = nz4genesMasked(org_data, imp_data, geneSet, masked_idxs)
    
    for i in range(0,len(geneSet)):
        org = org_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        imp = imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        if len(org) < 2:
            c = [np.nan, np.nan]
        else:
            c = stats.spearmanr(org,imp)
        corrs.loc[geneSet[i],'rho'] = c[0]
        corrs.loc[geneSet[i],'pvals'] = c[1]
    return(corrs)


################################################################################################
######################## Spot-based evals ################################################
######################## PCC (sw)
#Pearson's correlation coefficients (PCC) for each gene
#between the NON-ZERO original values and the predicted values
def eval_spotPCC(org_data, imp_data, geneSet):
    spots = np.array(org_data.index)
    idx = nz4spots(org_data, imp_data, geneSet, spots)
    
    corrs = pd.DataFrame(index = spots, columns=['corr_val'])
    for k,v in idx.items():
        #print(k)
        org_data.loc[k,[geneSet[a] for a in v]]
        c = org_data.loc[k,[geneSet[a] for a in v]].corr(imp_data.loc[k,[geneSet[a] for a in v]])
        corrs.loc[k,'corr_val'] = c
    return(corrs)


#Pearson's correlation coefficients (PCC) for each gene for masked experiment
#between the NON-ZERO original values and the predicted values
def eval_maskedspotPCC(org_data, imp_data, geneSet, masked_idxs):
    spots = np.array(org_data.index)
    idx = {}
    for i in geneSet:
        indices = masked_idxs.loc[i,:].dropna().tolist()
        idx[i] = [int(x) for x in indices]
    
    s = list(idx.values())[0]
    for i in range(1,len(list(idx.values()))):
        s = s + list(idx.values())[i]
    s = list(set(s))
    
    corrs = pd.DataFrame(index = s, columns=['corr_val'])
    
    for j in range(0,len(s)):
        genes = [k for k,v in idx.items() if s[j] in v]
        c = org_data.loc[spots[j],genes].corr(imp_data.loc[spots[j],genes])
        corrs.loc[s[j],'corr_val'] = c
    return(corrs)



######################## CosSim (sw)
#Cosine similarity for each gene
#between the NON-ZERO original values and the predicted values
def eval_spotCS(org_data, imp_data, geneSet):
    spots = np.array(org_data.index)
    idx = nz4spots(org_data, imp_data, geneSet, spots)
    
    cos_sims = pd.DataFrame(index = spots, columns=['cos_similarity'])
    for k,v in idx.items():
        #print(k)
        org = org_data.loc[k,[geneSet[a] for a in v]]
        imp = imp_data.loc[k,[geneSet[a] for a in v]]
        norm_sq = np.linalg.norm(org) * np.linalg.norm(imp)
        cos_sims.loc[k,'cos_similarity'] = (org @ imp) / norm_sq
    return(cos_sims)



#Cosine similarity for each gene for masked experiment
#between the NON-ZERO original values and the predicted values
def eval_maskedspotCS(org_data, imp_data, geneSet, masked_idxs):
    spots = np.array(org_data.index)
    idx = {}
    for i in geneSet:
        indices = masked_idxs.loc[i,:].dropna().tolist()
        idx[i] = [int(x) for x in indices]
    
    s = list(idx.values())[0]
    for i in range(1,len(list(idx.values()))):
        s = s + list(idx.values())[i]
    s = list(set(s))
    
    cos_sims = pd.DataFrame(index = s, columns=['cos_similarity'])
    for j in range(0,len(s)):
        genes = [k for k,v in idx.items() if s[j] in v]
        org = org_data.loc[spots[j],genes]
        imp = imp_data.loc[spots[j],genes]
        norm_sq = np.linalg.norm(org) * np.linalg.norm(imp)
        cos_sims.loc[s[j],'cos_similarity'] = (org @ imp) / norm_sq
    cos_sims.index = spots[cos_sims.index]
    return(cos_sims)


######################## RMSLE (sw)
#Root Mean Squared Log Error (RMSLE) for each gene
#between the NON-ZERO original values and the predicted values
def eval_spotRMSLE(org_data, imp_data, geneSet):
    spots = np.array(org_data.index)
    idx = nz4spots(org_data, imp_data, geneSet, spots)
    
    rmsles = pd.DataFrame(index = spots, columns=['rmsle_val'])
    for k,v in idx.items():
        org = org_data.loc[k,[geneSet[a] for a in v]]
        imp = imp_data.loc[k,[geneSet[a] for a in v]]
        rmsle = np.sqrt(np.square(np.log(imp + 1e-12) - np.log(org + 1e-12)).mean())
        rmsles.loc[k,'rmsle_val'] = rmsle

    return(rmsles)


#Root Mean Squared Log Error (RMSLE) for each gene for masked experiment
#between the NON-ZERO original values and the predicted values
def eval_maskedspotRMSLE(org_data, imp_data, geneSet, masked_idxs):
    spots = np.array(org_data.index)
    idx = {}
    for i in geneSet:
        indices = masked_idxs.loc[i,:].dropna().tolist()
        idx[i] = [int(x) for x in indices]
    
    s = list(idx.values())[0]
    for i in range(1,len(list(idx.values()))):
        s = s + list(idx.values())[i]
    s = list(set(s))
    
    rmsles = pd.DataFrame(index = s, columns=['rmsle_val'])
    for j in range(0,len(s)):
        genes = [k for k,v in idx.items() if s[j] in v]
        org = org_data.loc[spots[j],genes]
        imp = imp_data.loc[spots[j],genes]
        rmsle = np.sqrt(np.square(np.log(imp + 1e-12) - np.log(org + 1e-12)).mean())
        rmsles.loc[s[j],'rmsle_val'] = rmsle
    rmsles.index = spots[rmsles.index]
    return(rmsles)


