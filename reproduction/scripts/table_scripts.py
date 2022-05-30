# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 11:28:38 2022

@author: Gulben AVSAR

Codes for the tables
"""
import os
import numpy as np
import pandas as pd


def getGenes():
    df = pd.read_csv("genes_for_imputation.txt", sep='\t')
    df = df.set_index('Genes')
    df = df.sort_index(axis = 0)
    return(df)

# define the path for figures
# fig_path = './figures/'

#### Functions
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


def all_for_metrics(file_path, mtrcName, geneORspot):
    a = [i for i in os.listdir(file_path) if '{}_{}_pdacA'.format(geneORspot,mtrcName) in i]
    b = [i for i in os.listdir(file_path) if '{}_{}_pdacB'.format(geneORspot,mtrcName) in i]
    l = [i for i in os.listdir(file_path) if '{}_{}_Liver'.format(geneORspot,mtrcName) in i]
    all_mtrc = a + b + l
    mtrc = {}
    for i in all_mtrc:
        mtrc[i]=pd.read_csv(file_path+i, index_col=0)
    return(mtrc)


def meanMetricIndv_df(dataDict, sprsity_df, mthdsName, mtrcName):
    df = sprsity_df.copy()
    for i in mthdsName:
        df = pd.concat([dataDict['genewise_{}_pdacA_{}.txt'.format(mtrcName,i)],
                    dataDict['genewise_{}_pdacB_{}.txt'.format(mtrcName,i)],
                    dataDict['genewise_{}_Liver_{}.txt'.format(mtrcName,i)]], axis=1)
    df = df.dropna()
    return(df)


def perc_zeroPred(org_data, imp_data, geneSet):
    spots = np.array(org_data.index)
    idx = nz4genes(org_data, imp_data, geneSet)
    
    count = 0
    tot = 0
    for i in range(len(geneSet)):
        imp = imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        numb = imp[imp==0].count()
        count += numb
        tot += len(imp)
    
    return(count/tot)

def perc_maskedzeroPred(org_data, imp_data, geneSet, masked_idxs):
    spots = np.array(org_data.index)
    idx = nz4genesMasked(org_data, imp_data, geneSet)
    
    count = 0
    tot = 0
    for i in range(len(geneSet)):
        imp = imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        numb = imp[imp==0].count()
        count += numb
        tot += len(imp)
    
    return(count/tot)


def vol_dataPoints(org_data, imp_data, geneSet, mIdx):
    imp_data.index = org_data.index
    spots = np.array(org_data.index)
    if len(mIdx) == 3:
        idx = nz4genes(org_data, imp_data, geneSet)
    else:
        idx = nz4genesMasked(org_data, imp_data, geneSet, mIdx)
    tot = 0
    for i in range(len(geneSet)):
        imp = imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
        tot += len(imp)
    return(tot)


######################## Figure 1.b
def Fig1b (org_data, imp_data, mIdx, fig_path):
    geneSet = getGenes()
    geneSet = list(geneSet.index)
    
    mthds = ['spage', 'stplus', 'gimvi', 'tangram','5stlearn','30stlearn']
    labls = ['PDAC-A','PDAC-B','Liver']
    df = pd.DataFrame(index=mthds, columns=labls)
    
    for i in range(len(labls)):
        for j in range(len(mthds)):
            vol = vol_dataPoints(org_data[i][j], imp_data[i][j], geneSet, mIdx[i][j])
            df.loc[mthds[j],labls[i]] = vol
    df.to_csv(fig_path+'Figure1b_DataPointsVolume.txt')
    return(df)



######################## Table 1
def TabA(fig_path):
    Gene_set = getGenes()
    df = pd.DataFrame(data=Gene_set[['s_stA', 's_stB', 's_stL']].mean(axis=1),
                  index=Gene_set.index, columns=['sparsity'])
    file_path = './evaluations/'
    
    pcc = all_for_metrics(file_path, 'PCC', 'genewise')
    cs = all_for_metrics(file_path, 'CS', 'genewise')
    rmsle = all_for_metrics(file_path, 'RMSLE', 'genewise')
    
    mthds = ['PCC spage', 'PCC stplus', 'PCC gimvi', 'PCC tangram',
             'PCC 5stlearn','PCC 30stlearn', 'CS spage', 'CS stplus',
             'CS gimvi', 'CS tangram','CS 5stlearn','CS 30stlearn',
             'RMSLE spage', 'RMSLE stplus', 'RMSLE gimvi', 'RMSLE tangram',
             'RMSLE 5stlearn','RMSLE 30stlearn']
    
    labls = ['PDAC-A mean', 'PDAC-A std','PDAC-B mean', 'PDAC-B std',
             'Liver mean', 'Liver std']
    tabAll = pd.DataFrame(index=mthds, columns=labls)
    
    for i in mthds[0:6]:
        m = i.split(' ')[1]
        # meanMetricIndv_df(dataDict, sprsity_df, mthdsName, mtrcName)
        a = meanMetricIndv_df(pcc, df, [m], 'PCC')
        a.columns = ['A', 'B', 'L']
        tabAll.loc[i,'PDAC-A mean'] = a.mean()[0]
        tabAll.loc[i,'PDAC-B mean'] = a.mean()[1]
        tabAll.loc[i,'Liver mean'] = a.mean()[2]
        tabAll.loc[i,'PDAC-A std'] = a.std()[0]
        tabAll.loc[i,'PDAC-B std'] = a.std()[1]
        tabAll.loc[i,'Liver std'] = a.std()[2]
    
    for i in mthds[6:12]:
        m = i.split(' ')[1]
        a = meanMetricIndv_df(cs, df, [m], 'CS')
        a.columns = ['A', 'B', 'L']
        tabAll.loc[i,'PDAC-A mean'] = a.mean()[0]
        tabAll.loc[i,'PDAC-B mean'] = a.mean()[1]
        tabAll.loc[i,'Liver mean'] = a.mean()[2]
        tabAll.loc[i,'PDAC-A std'] = a.std()[0]
        tabAll.loc[i,'PDAC-B std'] = a.std()[1]
        tabAll.loc[i,'Liver std'] = a.std()[2]
    
    for i in mthds[12:18]:
        m = i.split(' ')[1]
        a = meanMetricIndv_df(rmsle, df, [m], 'RMSLE')
        a.columns = ['A', 'B', 'L']
        tabAll.loc[i,'PDAC-A mean'] = a.mean()[0]
        tabAll.loc[i,'PDAC-B mean'] = a.mean()[1]
        tabAll.loc[i,'Liver mean'] = a.mean()[2]
        tabAll.loc[i,'PDAC-A std'] = a.std()[0]
        tabAll.loc[i,'PDAC-B std'] = a.std()[1]
        tabAll.loc[i,'Liver std'] = a.std()[2]

    tabAll.to_csv(fig_path+'Table1_EvaluationMetrics.txt')    
    
    return(tabAll)

######################## Table 2
def TabB (org_data, imp_data, mIdx, fig_path):
    
    geneSet = getGenes()
    geneSet = list(geneSet.index)
    
    mthds = ['spage', 'stplus', 'gimvi', 'tangram','5stlearn','30stlearn']
    labls = ['PDAC-A','PDAC-B','Liver']
    df = pd.DataFrame(index=mthds, columns=labls)
    
    for n in range(3):
        org_data[n][2].index = org_data[n][0].index
        imp_data[n][2].index = imp_data[n][0].index
    
    for i in range(len(labls)):
        spots = np.array(org_data[i][0].index)
        for j in range(len(mthds)):
            
            if len(mIdx[i][j]) == 3:
                idx = nz4genes(org_data[i][j], imp_data[i][j], geneSet)
            else:
                idx = nz4genesMasked(org_data[i][j], imp_data[i][j], geneSet, mIdx[i][j])
            
            count = 0
            tot = 0
            for k in range(len(geneSet)):
                imp = imp_data[i][j].loc[spots[idx[geneSet[k]]],geneSet[k]]
                numb = imp[imp==0].count()
                count += numb
                tot += len(imp)
            df.loc[mthds[j],labls[i]] = count/tot
    
    df.to_csv(fig_path+'Table2_ZeroPredPerc.txt')
    return(df)


########################################################################
######################## Extra in text
# predictions with only at most 3% error rate
# find the number of spots that is the original and imp data are too close
# with a defined error ratio
def find_numbSpots (org_data, imp_data, geneName, Eratio):
    g = geneName
    x = org_data.loc[:,g].copy()
    x1 = x - x*Eratio
    x2 = x + x*Eratio
    y = imp_data.loc[:,g].copy()
    y = y[x1 < y]
    x2 = x2[y.index]
    y = y[y < x2]
    
    return(len(y)/len(org_data))


def preds_w_Eratio(org_data, imp_data, Eratio, geneSet):
    if len(geneSet) == 0:
        geneSet = getGenes()
        geneSet = list(geneSet.index)
    
    imp_data[0].index = org_data[0].index
    org_data[1].index= org_data[0].index
    imp_data[1].index = org_data[0].index
    
    df = pd.DataFrame(index= geneSet, columns=['in_stplus', 'in_gimvi'])
    for g in df.index:
        df.loc[g,'in_stplus'] = find_numbSpots(org_data[0], imp_data[0], g, 0.03)
        df.loc[g,'in_gimvi'] = find_numbSpots(org_data[1], imp_data[1], g, 0.03)
    
    return(print(str('\n#######################☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻##################')+
                 '\nRatio of predictions with at most '+str(Eratio*100)+'% error rate:\n'+
                  '      stPlus --->  '+str(df.in_stplus.mean())+'\n'+
                  '      gimVI --->  '+str(df.in_gimvi.mean())+'\n'+
                  str('#######################☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻☻##################')))

########################################################################
















