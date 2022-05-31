# -*- coding: utf-8 -*-
"""
@author: Gulben AVSAR
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def getGenes():
    df = pd.read_csv("genes_for_imputation.txt", sep='\t')
    df = df.set_index('Genes')
    df = df.sort_index(axis = 0)
    return(df)


##### Figure 1.a
def Fig1a(fig_path):
    Gene_set = getGenes()
    
    sparsity = [Gene_set['s_stA'].to_list(),Gene_set['s_stB'].to_list(),Gene_set['s_stL'].to_list()]

    np.random.seed(20)
    X = list(np.random.uniform(0.5, 3.5, 8))
    X.extend([0.5, 3.5])
    pos = [1, 2, 3]
    
    fig, ax = plt.subplots(tight_layout=True)
    parts = ax.violinplot(sparsity, pos, showmedians=False)
    for pc in parts['bodies']:
        pc.set_facecolor('cornflowerblue')
        pc.set_edgecolor('black')
        pc.set_linewidth(0.5)
        pc.set_alpha(0.7)
    parts['cbars'].set_color('black')
    parts['cmaxes'].set_color('black')
    parts['cmins'].set_color('black')
    #parts['cmedians'].set_color('white')
    parts['cbars'].set_linewidth(0.5)
    #parts['cmedians'].set_linewidth(0.8)
    medians = np.percentile(sparsity, [50], axis=1)
    ax.scatter(pos, medians, marker='o', color='white', s=30, zorder=3)
    ax.plot(sorted(X), [0.97]*10, '--', color='black', markersize=0.005)
    ax.plot(sorted(X), [0.90]*10, '--', color='black', markersize=0.005)
    ax.plot(sorted(X), [0.75]*10, '--', color='black', markersize=0.005)
    ax.set_xticks(pos)
    ax.set_xticklabels(['PDAC-A', 'PDAC-B', 'Liver'])
    ax.set_ylabel('sparsity')
    plt.savefig(fig_path +'/Figure1a_Sparsity_'+str(600)+'dpi.tiff', dpi=600, format="tiff")
    return(fig)


##### Figure 2
# dict for the metric files
def all_for_metrics(file_path, mtrcName, geneORspot):
    a = [i for i in os.listdir(file_path) if '{}_{}_pdacA'.format(geneORspot,mtrcName) in i]
    b = [i for i in os.listdir(file_path) if '{}_{}_pdacB'.format(geneORspot,mtrcName) in i]
    l = [i for i in os.listdir(file_path) if '{}_{}_Liver'.format(geneORspot,mtrcName) in i]
    all_mtrc = a + b + l
    mtrc = {}
    for i in all_mtrc:
        mtrc[i]=pd.read_csv(file_path+i, index_col=0)
    return(mtrc)


# mean of metric from all datasets
def meanMetric_df(dataDict, sprsity_df, mthdsName, mtrcName):
    df = sprsity_df.copy()
    for i in mthdsName:
        mtrc_df = pd.concat([dataDict['genewise_{}_pdacA_{}.txt'.format(mtrcName,i)],
                    dataDict['genewise_{}_pdacB_{}.txt'.format(mtrcName,i)],
                    dataDict['genewise_{}_Liver_{}.txt'.format(mtrcName,i)]], axis=1)
        df['{}'.format(i)] = mtrc_df.mean(axis = 1, skipna = True)
    df = df.dropna()
    df = df.sort_values('sparsity')
    return(df)



# plot sparsity vs metric for each tool
def pl_sprsVSmtrc(df, mthdsName, sprsLim, mtrcLim, xlabel, ylabel, figSize):
    fig = plt.figure(tight_layout=True, figsize=figSize)
    plt.subplots_adjust(wspace= 0.55, hspace= 0.55)
    for i in range(0,len(mthdsName)):
        sub=fig.add_subplot(2,3,i+1)
        sub.plot(df['sparsity'], df[mthdsName[i]], 'ro', markersize=2.5)
        sub.plot(df['sparsity'], [sprsLim]*len(df), '--', color='dimgray', markersize=0.0005)
        xline = [min(df[mthdsName[i]])-0.001, max(df[mthdsName[i]])+0.001]
        xline.extend(sorted(np.random.uniform(min(df[mthdsName[i]]),max(df[mthdsName[i]]), len(df))))
        sub.plot([mtrcLim]*len(xline), sorted(xline), '--', color='dimgray', markersize=0.0005)
        # sub.set_ylim(ylims[i])
        sub.set_title(mthdsName[i], fontweight='bold', fontsize=20)
        sub.set_xlabel(xlabel, fontsize=15)
        sub.set_ylabel(ylabel, fontsize=15)
        sub.tick_params(axis='y', labelsize=15)
        sub.tick_params(axis='x', labelsize=15)
    return(fig)


def Fig2(fig_path):
    Gene_set = getGenes()
    df = pd.DataFrame(data=Gene_set[['s_stA', 's_stB', 's_stL']].mean(axis=1),
                  index=Gene_set.index, columns=['sparsity'])
    file_path = './evaluations/'
    mthds = ['spage', 'stplus', 'gimvi', 'tangram', '5stlearn', '30stlearn']
    
    pcc = all_for_metrics(file_path, 'PCC', 'genewise')
    cs = all_for_metrics(file_path, 'CS', 'genewise')
    rmsle = all_for_metrics(file_path, 'RMSLE', 'genewise')
    
    ##### mean of metric from all datasets for each gene
    df_pcc =  meanMetric_df(pcc, df, mthds, 'PCC')
    df_cs =  meanMetric_df(cs, df, mthds, 'CS')
    df_rmsle =  meanMetric_df(rmsle, df, mthds, 'RMSLE')
    
    fig_genepcc = pl_sprsVSmtrc(df_pcc, mthds, 0, 0.90, 'sparsity', 'PCC', (16,8))
    fig_geneCosSim = pl_sprsVSmtrc(df_cs, mthds, 0.70, 0.90, 'sparsity', 'CS', (16,8))
    fig_genermsle = pl_sprsVSmtrc(df_rmsle, mthds, 0, 0.90, 'sparsity', 'RMSLE', (16,8))
    
    fig_genepcc.savefig(fig_path+'Figure2a_genewise_PCC.tiff', dpi=600, format="tiff")
    fig_geneCosSim.savefig(fig_path+'Figure2b_genewise_CS.tiff', dpi=600, format="tiff")
    fig_genermsle.savefig(fig_path+'Figure2c_genewise_RMSLE.tiff', dpi=600, format="tiff")    
    
    return()


##### Figure 3
def FigC(org_data, imp_data, mIdx, dName, mName, fig_size, a_org, a_imp, fig_path):
    # if there is only one value, it plots barplot otherwise filled scatter plot is done
    Gene_set = getGenes()
    np.random.seed(10)
    geneSet = pd.DataFrame()
    for i in [1,2,3,4]:
        s = Gene_set[Gene_set['Groups'] == i].iloc[np.random.choice(np.arange(len(Gene_set[Gene_set['Groups'] == i])), size=1)]
        geneSet = geneSet.append(s, ignore_index=False)
    del i, s
    geneSet = list(geneSet.index)
        
    if mName in ['spage', 'stplus', 'gimvi', 'tangram']:
        gene_idx_dict = {} #indices of non-zero values for given genes
        for i in geneSet:
            if i not in org_data.columns:
                print(i, '    not found in ST dataset')
                continue
            else:
                indices = [idx for idx, element in enumerate(org_data.loc[:,i]) if element != 0]
                gene_idx_dict[i] = indices
    
        fig = plt.figure(tight_layout=True, figsize=fig_size)
        plt.subplots_adjust(wspace= 0.25, hspace= 0.25)
        for i in range(0,len(geneSet)):
            sub=fig.add_subplot(1,4,i+1)
            cmap=['brown', 'teal']
            ld = gene_idx_dict[geneSet[i]]
            if len(ld) <= 5:
                 X = np.arange(len(ld))
                 org = org_data.iloc[ld,org_data.columns == geneSet[i]][geneSet[i]].values
                 imp = imp_data.iloc[ld,imp_data.columns == geneSet[i]][geneSet[i]].values
                 sub.bar(X, org, 0.10, color='lightcoral', alpha=a_org)
                 sub.bar(X, imp, 0.10, color='turquoise', alpha=a_imp)
                 sub.set_title(geneSet[i], fontsize=20)
                 sub.set_xlim([-0.5, len(X)+0.05])
                 sub.set_xticks(np.arange(len(X)))
                 sub.tick_params(axis='y', labelsize=18)
                 sub.tick_params(axis='x', labelsize=18)
                 sub.grid(False)
            else:
                 org = org_data.iloc[ld,org_data.columns == geneSet[i]]
                 imp = imp_data.iloc[ld,imp_data.columns == geneSet[i]]
                 sub.plot(range(len(ld)), org,
                          'o-', c=cmap[0], markersize=0.5)
                 sub.plot(range(len(ld)), imp,
                          'o-', c=cmap[1], markersize=0.5)
                 sub.fill_between(range(len(ld)),
                                  org.to_numpy().reshape(1,-1)[0],
                                  0, color='lightcoral',       # The outline color
                                  alpha=0.85)          # Transparency of the fill
                 sub.fill_between(range(len(ld)),
                                  imp.to_numpy().reshape(1,-1)[0],
                                  0, color='turquoise',       # The outline color
                                  alpha=0.2)          # Transparency of the fill
                 sub.set_title(geneSet[i], fontsize=20)
                 sub.tick_params(axis='y', labelsize=18)
                 sub.tick_params(axis='x', labelsize=18)
                 sub.set_facecolor('white')
            if (mName == 'spage') and (i == len(geneSet)-1):
                sub.legend(['original','imputed'],bbox_to_anchor=(1.05, 1),
                           loc='upper left', fontsize=20)
    elif mName in ['5stlearn','30stlearn']:
        
        mIdxs={}
        for i in geneSet:
            indices = mIdx.loc[i,:].dropna().tolist()
            mIdxs[i] = indices
        
        fig = plt.figure(tight_layout=True, figsize=fig_size)
        plt.subplots_adjust(wspace= 0.25, hspace= 0.25)
        for i in range(0,len(geneSet)):
            sub=fig.add_subplot(1,4,i+1)
            cmap=['brown', 'teal']
            ld = mIdxs[geneSet[i]]
            if len(ld) <= 5:
                 X = np.arange(len(ld))
                 org = org_data.iloc[ld,org_data.columns == geneSet[i]][geneSet[i]].values
                 imp = imp_data.iloc[ld,imp_data.columns == geneSet[i]][geneSet[i]].values
                 sub.bar(X, org, 0.10, color='lightcoral', alpha=a_org)
                 sub.bar(X, imp, 0.10, color='turquoise', alpha=a_imp)
                 sub.set_title(geneSet[i], fontsize=20)
                 sub.set_xlim([-0.5, len(X)+0.05])
                 sub.set_xticks(np.arange(len(X)))
                 sub.tick_params(axis='y', labelsize=18)
                 sub.tick_params(axis='x', labelsize=18)
                 sub.grid(False)
            else:
                 org = org_data.iloc[ld,org_data.columns == geneSet[i]]
                 imp = imp_data.iloc[ld,imp_data.columns == geneSet[i]]
                 sub.plot(range(len(ld)), org,
                          'o-', c=cmap[0], markersize=0.5)
                 sub.plot(range(len(ld)), imp,
                          'o-', c=cmap[1], markersize=0.5)
                 sub.fill_between(range(len(ld)),
                                  org.to_numpy().reshape(1,-1)[0],
                                  0, color='lightcoral',       # The outline color
                                  alpha=0.85)          # Transparency of the fill
                 sub.fill_between(range(len(ld)),
                                  imp.to_numpy().reshape(1,-1)[0],
                                  0, color='turquoise',       # The outline color
                                  alpha=0.2)          # Transparency of the fill
                 sub.set_title(geneSet[i], fontsize=20)
                 sub.tick_params(axis='y', labelsize=18)
                 sub.tick_params(axis='x', labelsize=18)
                 sub.set_facecolor('white')
            # if i == len(geneSet)-1:
            #     sub.legend(['original','imputed'],bbox_to_anchor=(1.05, 1),
            #                loc='upper left', fontsize=20)
    
    if dName == 'PDAC-A':
        dName = dName.split('-')[0].lower()+dName.split('-')[1]
        fig.savefig(fig_path+'Figure3parts_{}_{}.tiff'.format(dName,mName), dpi=600, format="tiff")
    elif dName == 'PDAC-B':
        dName = dName.split('-')[0].lower()+dName.split('-')[1]
        fig.savefig(fig_path+'SupplFigure3parts_{}_{}.tiff'.format(dName,mName), dpi=600, format="tiff")
    elif dName == 'Liver':
        fig.savefig(fig_path+'SupplFigure4parts_{}_{}.tiff'.format(dName,mName), dpi=600, format="tiff")
    
    return(fig)


##### Figure 4
def FigD(org_data, imp_data, mIdx, dName, mName, fig_size, PointColor, LineColor, fig_path):
    Gene_set = getGenes()
    np.random.seed(10)
    geneSet = pd.DataFrame()
    for i in [3,4]:
        s = Gene_set[Gene_set['Groups'] == i].iloc[np.random.choice(np.arange(len(Gene_set[Gene_set['Groups'] == i])), size=2)]
        geneSet = geneSet.append(s, ignore_index=False)
    del i, s
    geneSet = list(geneSet.index)
        
    if mName in ['spage', 'stplus', 'gimvi', 'tangram']:
        idx = {} #indices of non-zero values for given genes
        for i in geneSet:
            if i not in org_data.columns:
                print(i, '    not found in ST dataset')
                continue
            else:
                indices = [idx for idx, element in enumerate(org_data.loc[:,i]) if element != 0]
                idx[i] = indices
        
        imp_data.index = org_data.index
        spots = np.array(org_data.index)
        
        fig = plt.figure(tight_layout=True, figsize=fig_size)
        plt.subplots_adjust(wspace= 0.25, hspace= 0.25)
        for i in range(len(geneSet)):
            sub=fig.add_subplot(1,4,i+1)
            org = org_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
            imp = imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
            sub.scatter(org, imp, c=PointColor)
            sub.plot(org, org, c=LineColor)
            sub.tick_params(axis='y', labelsize=15)
            sub.tick_params(axis='x', labelsize=15)
            sub.set_xlabel('original', fontsize=15)
            sub.set_ylabel('imputed', fontsize=15)
            sub.set_title(geneSet[i], fontsize=20)
    
    elif mName in ['5stlearn','30stlearn']:
        
        mIdxs={}
        for i in geneSet:
            indices = mIdx.loc[i,:].dropna().tolist()
            mIdxs[i] = indices
        
        imp_data.index = org_data.index
        
        fig = plt.figure(tight_layout=True, figsize=fig_size)
        plt.subplots_adjust(wspace= 0.25, hspace= 0.25)
        for i in range(len(geneSet)):
            sub=fig.add_subplot(1,4,i+1)
            ld = idx[geneSet[i]]
            org = org_data.iloc[ld,org_data.columns == geneSet[i]]
            imp = imp_data.iloc[ld,org_data.columns == geneSet[i]]
            sub.scatter(org, imp, c=PointColor)
            sub.plot(org, org, c=LineColor)
            sub.tick_params(axis='y', labelsize=15)
            sub.tick_params(axis='x', labelsize=15)
            sub.set_xlabel('original', fontsize=15)
            sub.set_ylabel('imputed', fontsize=15)
            sub.set_title(geneSet[i], fontsize=20)
    
    if dName == 'PDAC-A':
        dName = dName.split('-')[0].lower()+dName.split('-')[1]
        fig.savefig(fig_path+'Figure4_{}_{}.tiff'.format(dName,mName), dpi=600, format="tiff")
    elif dName == 'PDAC-B':
        dName = dName.split('-')[0].lower()+dName.split('-')[1]
        fig.savefig(fig_path+'Figure4_{}_{}.tiff'.format(dName,mName), dpi=600, format="tiff")
    elif dName == 'Liver':
        fig.savefig(fig_path+'Figure4_{}_{}.tiff'.format(dName,mName), dpi=600, format="tiff")

    return(fig)


##### Suppl.Figure 5
def SupplFigC(org_data, imp_data, mIdx, geneName, dName, mName, dataLims,
              fig_size, PointColor, LineColor, PointAlpha, fig_path):
    g = geneName
    idx = {} #indices of non-zero values for given genes
    indices = [idx for idx, element in enumerate(org_data.loc[:,g]) if element != 0]
    idx[g] = indices
    imp_data.index = org_data.index
    spots = np.array(org_data.index)
    
    org = org_data.loc[spots[idx[g]],g]
    imp = imp_data.loc[spots[idx[g]],g]
        
    org = org[dataLims[0]<org]
    org = org[org<dataLims[1]]
    imp = imp[org.index]
    
    plt.figure(tight_layout=True, figsize=fig_size)
    plt.scatter(org, imp, c=PointColor, s=100, alpha=PointAlpha)
    plt.plot(org, org, c=LineColor)
    plt.tick_params(axis='y', labelsize=15)
    plt.tick_params(axis='x', labelsize=15)
    plt.xlabel('original', fontsize=15)
    plt.ylabel('imputed', fontsize=15)
    plt.title(g, fontsize=20)
    plt.savefig(fig_path+'SupplFigure5_{}_{}_{}.tiff'.format(dName,mName,g),
                dpi=600, format="tiff")
    
    return(plt.show())



##### Figure 5
def plot4spotMetrics(dataDict, mName, mtrcLim, figsize, yLimits):
    dataA = [v.values.reshape(-1,) for k,v in dataDict.items() if 'A' in k]
    dataA = [dataA[3][~np.isnan(dataA[3])], dataA[4][~np.isnan(dataA[4])],
             dataA[2][~np.isnan(dataA[2])], dataA[5][~np.isnan(dataA[5])],
             dataA[1][~np.isnan(dataA[1])], dataA[0][~np.isnan(dataA[0])]]
    dataB = [v.values.reshape(-1,) for k,v in dataDict.items() if 'B' in k]
    dataB = [dataB[3][~np.isnan(dataB[3])], dataB[4][~np.isnan(dataB[4])],
             dataB[2][~np.isnan(dataB[2])], dataB[5][~np.isnan(dataB[5])],
             dataB[1][~np.isnan(dataB[1])], dataB[0][~np.isnan(dataB[0])]]
    dataL = [v.values.reshape(-1,) for k,v in dataDict.items() if 'L' in k]
    dataL = [dataL[3][~np.isnan(dataL[3])], dataL[4][~np.isnan(dataL[4])],
             dataL[2][~np.isnan(dataL[2])], dataL[5][~np.isnan(dataL[5])],
             dataL[1][~np.isnan(dataL[1])], dataL[0][~np.isnan(dataL[0])]]
    data = dataA + dataB + dataL
    labls = [k for k,v in dataDict.items() if 'A' in k]
    labls = [labls[3], labls[4], labls[2], labls[5], labls[1], labls[0]]
    labls = [i.split('_')[3] for i in labls]
    labls = [i.split('.')[0] for i in labls]*3
    colors = ['#0b9deb','#427def','#7773f3','#9268f0','#a553d5','purple']*3
    locs = [1,2,3,4,5,6]
    xlocs = locs +  list(map(lambda x:x+8, locs)) + list(map(lambda x:x+16, locs))

    fig, ax = plt.subplots(tight_layout=True, figsize=figsize)
    parts= ax.violinplot(data, xlocs, showmedians=True)
    count=0
    for pc in parts['bodies']:
        pc.set_facecolor(colors[count])
        pc.set_edgecolor('black')
        pc.set_linewidth(0.5)
        pc.set_alpha(0.7)
        count += 1
    parts['cbars'].set_color('black')
    parts['cmaxes'].set_color('black')
    parts['cmins'].set_color('black')
    parts['cmedians'].set_color('black')
    parts['cbars'].set_linewidth(0.5)
    parts['cmedians'].set_linewidth(0.8)
    xline = yLimits
    xline.extend(np.random.uniform(xline[0],xline[1], 28))
    ax.plot([7.5]*30, sorted(xline), '--', color='dimgray', markersize=0.0005)
    ax.plot([15.5]*30, sorted(xline), '--', color='dimgray', markersize=0.0005)
    xline = [0, max(xlocs)+1]
    xline.extend(np.random.uniform(0, max(xlocs), 28))
    ax.plot(sorted(xline),[mtrcLim]*len(xline), '--', color='dimgray', markersize=0.0005)
    ax.set_xlim(0, max(xlocs) + 1)
    ax.set_xticks(xlocs)
    ax.set_xticklabels(labls, rotation=45 , fontsize=20)
    #ax.set_ylabel(mName, fontsize=20)
    ax.tick_params(axis='y', labelsize=18)
    
    return(fig)


def FigE(fig_path):
    file_path = './evaluations/'
    
    pcc = all_for_metrics(file_path, 'PCC', 'spotwise')
    cs = all_for_metrics(file_path, 'CS', 'spotwise')
    rmsle = all_for_metrics(file_path, 'RMSLE', 'spotwise')
    
    allData = [pcc, cs, rmsle]
    labl = ['PCC', 'CS', 'RMSLE']
    mLims = [0.0, 0.7, 0.0]
    yLims = [[-1.1, 1.05], [0, 1.05], [0, 1.62]]
    for i in range(len(allData)):
        fig = plot4spotMetrics(allData[i], labl[i], mLims[i], (18,7), yLims[i])
        fig.savefig(fig_path+'Figure5_{}.tiff'.format(labl[i]), dpi=600, format='tiff')
    

##### Figure 6
def plot4synST(data,dataerr, mthds):
    
    data['color'] = ['#0b9deb','#427def','#7773f3','#9268f0']
    x = range(4)
    fig = plt.figure(tight_layout=True, figsize=(18,7))
    ax1 = plt.subplot(131)
    ax1.bar(x, data['PCC'], yerr=dataerr['PCC'], width=0.8, color=data['color'],
            align='center', ecolor='black', capsize=10)
    plt.xticks(x, mthds, fontsize=20)
    plt.ylabel('PCC', fontsize=20)
    ax2 = plt.subplot(132)
    ax2.bar(x, data['CS'], yerr=dataerr['CS'], width=0.8, color=data['color'],
            align='center', ecolor='black', capsize=10)
    plt.xticks(x, mthds, fontsize=20)
    plt.ylabel('CS', fontsize=20)
    ax3 = plt.subplot(133)
    ax3.bar(x, data['RMSLE'], yerr=dataerr['RMSLE'], width=0.8, color=data['color'],
            align='center', ecolor='black', capsize=10)
    plt.xticks(x, mthds, fontsize=20)
    plt.ylabel('RMSLE', fontsize=20)
    
    return(fig)


def metrics4Syn(file_path, mtrcName, geneORspot):
    s = [i for i in os.listdir(file_path) if '{}_{}_SimulatedST'.format(geneORspot,mtrcName) in i]
    mtrc = {}
    for i in s:
        mtrc[i]=pd.read_csv(file_path+i, index_col=0)
    return(mtrc)


def FigF(fig_path):
    file_path = './evaluations/'
    mthds = ['spage', 'stplus', 'gimvi', 'tangram']
    
    genepcc = metrics4Syn(file_path, 'PCC', 'genewise')
    genecs = metrics4Syn(file_path, 'CS', 'genewise')
    genermsle = metrics4Syn(file_path, 'RMSLE', 'genewise')
    
    spotpcc = metrics4Syn(file_path, 'PCC', 'spotwise')
    spotcs = metrics4Syn(file_path, 'CS', 'spotwise')
    spotrmsle = metrics4Syn(file_path, 'RMSLE', 'spotwise')
    
    # for gene-wise
    labels = ['PCC', 'CS', 'RMSLE']
    data = pd.DataFrame(index=mthds, columns=labels)
    m = 0
    for i in [1,2,0,3]:
        data.loc[mthds[m],'PCC'] = genepcc[list(genepcc.keys())[i]].mean(skipna=True).values[0]
        data.loc[mthds[m],'CS'] = genecs[list(genecs.keys())[i]].mean(skipna=True).values[0]
        data.loc[mthds[m],'RMSLE'] = genermsle[list(genermsle.keys())[i]].mean(skipna=True).values[0]
        m += 1
    
    err = pd.DataFrame(index=mthds, columns=labels)
    m = 0
    for i in [1,2,0,3]:
        err.loc[mthds[m],'PCC'] = genepcc[list(genepcc.keys())[i]].std(skipna=True).values[0]
        err.loc[mthds[m],'CS'] = genecs[list(genecs.keys())[i]].std(skipna=True).values[0]
        err.loc[mthds[m],'RMSLE'] = genermsle[list(genermsle.keys())[i]].std(skipna=True).values[0]
        m += 1
    
    fig_gw = plot4synST(data,err,mthds)
    
    # for spot-wise
    labels = ['PCC', 'CS', 'RMSLE']
    data = pd.DataFrame(index=mthds, columns=labels)
    m = 0
    for i in [1,2,0,3]:
        data.loc[mthds[m],'PCC'] = spotpcc[list(spotpcc.keys())[i]].mean(skipna=True).values[0]
        data.loc[mthds[m],'CS'] = spotcs[list(spotcs.keys())[i]].mean(skipna=True).values[0]
        data.loc[mthds[m],'RMSLE'] = spotrmsle[list(spotrmsle.keys())[i]].mean(skipna=True).values[0]
        m += 1
    
    
    err = pd.DataFrame(index=mthds, columns=labels)
    m = 0
    for i in [1,2,0,3]:
        err.loc[mthds[m],'PCC'] = spotpcc[list(spotpcc.keys())[i]].std(skipna=True).values[0]
        err.loc[mthds[m],'CS'] = spotcs[list(spotcs.keys())[i]].std(skipna=True).values[0]
        err.loc[mthds[m],'RMSLE'] = spotrmsle[list(spotrmsle.keys())[i]].std(skipna=True).values[0]
        m += 1
    
    fig_sw = plot4synST(data,err,mthds)
    
    # save figs
    fig_gw.savefig(fig_path+'Figure6a_genewise_SimulatedST.tiff', dpi=600, format="tiff")
    fig_sw.savefig(fig_path+'Figure6b_spotwise_SimulatedST.tiff', dpi=600, format="tiff")    
    
    return(fig_gw,fig_sw)

##### Suppl.Figure 2
def SupplFigD(mtrcName,cutOff,xName, fig_path):
    Gene_set = getGenes()
    Gene_set = Gene_set.drop(['s_stA','s_stB','s_stL','s_stSyn'],axis=1)
    file_path = './evaluations/'
    # mthds = ['spage', 'stplus', 'gimvi', 'tangram', '5stlearn', '30stlearn']
    
    fname = [i for i in os.listdir(file_path) if 'genewise_{}pvals_Liver'.format(mtrcName) in i]
    data=pd.read_csv(file_path+fname[0], index_col=0)
    data['sparsity_group'] = Gene_set['Groups']
    data = data.fillna(1)
    
    
    data['colors'] = ['grey']*len(data)
    data.loc[(data[xName] >= cutOff) & (data['pvals'] < 0.05),'colors'] = 'red'
    
    #plotting
    fig = plt.figure(tight_layout=True, figsize=(18,7))
    for i in range(1,5):
        ax = plt.subplot(1,4,i)
        ax.scatter(data.loc[(data['sparsity_group'] == i),xName],
                -(np.log10(data.loc[(data['sparsity_group'] == i),'pvals'])),
                  c=data.loc[(data['sparsity_group'] == i),'colors'], s=20)
        # ax.set_ylim(ymin=0)
        # ax.spines['left'].set_position('zero')
        ax.tick_params(axis='both', which='major', labelsize=18)
        plt.xlabel(mtrcName, fontsize=25)
        if i == 1:
            plt.ylabel('-log10(p)', fontsize=25)
    
    fig.savefig(fig_path+'SupplFigure2_genewise_{}pvals_Liver_gimvi.tiff'.format(mtrcName),
                dpi=600, format="tiff")    
    
    return()



########################################################################
######################## Extra in text

def below_metric(mtrcName, mthdName, cutOff):
    file_path = './evaluations/'
    i = mthdName
    mtrc_dict = all_for_metrics(file_path, mtrcName, 'genewise')
    mtrc_df = pd.concat([mtrc_dict['genewise_{}_pdacA_{}.txt'.format(mtrcName,i)],
                    mtrc_dict['genewise_{}_pdacB_{}.txt'.format(mtrcName,i)],
                    mtrc_dict['genewise_{}_Liver_{}.txt'.format(mtrcName,i)]], axis=1)
    mtrc_df = mtrc_df.mean(axis = 1, skipna = True)
    
    belows = mtrc_df[mtrc_df <= cutOff]
    
    return((len(belows)/len(mtrc_df))*100)


    




########################################################################
########################################################################
########################################## QQeval
##### Figure 3
def comparison_pl(org_data, imp_data, g, dName, fig_path,
                  mName=None, pl_g=None, fig_size=(18,4),
                  a_org=1, a_imp=0.5):
    # if there is only one value, it plots barplot otherwise filled scatter plot is done
    
    
    if pl_g is comparison_pl.__defaults__[1]:
        geneSet = list( g[i] for i in np.random.choice(np.arange(len(g)), size=4))
    else:
        geneSet = pl_g
    
    
    gene_idx_dict = {} #indices of non-zero values for given genes
    for i in geneSet:
        if i not in org_data.columns:
            print(i, '    not found in ST dataset')
            continue
        else:
            indices = [idx for idx, element in enumerate(org_data.loc[:,i]) if element != 0]
            gene_idx_dict[i] = indices    
    
    fig = plt.figure(tight_layout=True, figsize=fig_size)
    plt.subplots_adjust(wspace= 0.25, hspace= 0.25)
    for i in range(0,len(geneSet)):
        sub=fig.add_subplot(1,4,i+1)
        cmap=['brown', 'teal']
        ld = gene_idx_dict[geneSet[i]]
        if len(ld) <= 5:
            X = np.arange(len(ld))
            org = org_data.iloc[ld,org_data.columns == geneSet[i]][geneSet[i]].values
            imp = imp_data.iloc[ld,imp_data.columns == geneSet[i]][geneSet[i]].values
            sub.bar(X, org, 0.10, color='lightcoral', alpha=a_org)
            sub.bar(X, imp, 0.10, color='turquoise', alpha=a_imp)
            sub.set_title(geneSet[i], fontsize=20)
            sub.set_xlim([-0.5, len(X)+0.05])
            sub.set_xticks(np.arange(len(X)))
            sub.tick_params(axis='y', labelsize=18)
            sub.tick_params(axis='x', labelsize=18)
            sub.grid(False)
        else:
            org = org_data.iloc[ld,org_data.columns == geneSet[i]]
            imp = imp_data.iloc[ld,imp_data.columns == geneSet[i]]
            sub.plot(range(len(ld)), org,
                     'o-', c=cmap[0], markersize=0.5)
            sub.plot(range(len(ld)), imp,
                     'o-', c=cmap[1], markersize=0.5)
            sub.fill_between(range(len(ld)),
                             org.to_numpy().reshape(1,-1)[0],
                             0, color='lightcoral',       # The outline color
                             alpha=0.85)          # Transparency of the fill
            sub.fill_between(range(len(ld)),
                             imp.to_numpy().reshape(1,-1)[0],
                             0, color='turquoise',       # The outline color
                             alpha=0.2)          # Transparency of the fill
            sub.set_title(geneSet[i], fontsize=20)
            sub.tick_params(axis='y', labelsize=18)
            sub.tick_params(axis='x', labelsize=18)
            sub.set_facecolor('white')
        if i == len(geneSet)-1:
            sub.legend(['original','imputed'],bbox_to_anchor=(1.05, 1),
                       loc='upper left', fontsize=20)

    fig.savefig(fig_path+'comparison_pl_{}_{}.tiff'.format(dName,mName), dpi=600, format="tiff")
    
    return(fig)


##### Figure 4
def linearity_pl(org_data, imp_data, g, dName, fig_path,
                  mName=None, pl_g=None, fig_size=(18,4),
                  PointColor='brown', LineColor='dimgrey'):
    # if there is only one value, it plots barplot otherwise filled scatter plot is done
    
    
    if pl_g is comparison_pl.__defaults__[1]:
        geneSet = list( g[i] for i in np.random.choice(np.arange(len(g)), size=4))
    else:
        geneSet = pl_g
    
    
    idx = {} #indices of non-zero values for given genes
    for i in geneSet:
        if i not in org_data.columns:
            print(i, '    not found in ST dataset')
            continue
        else:
            indices = [idx for idx, element in enumerate(org_data.loc[:,i]) if element != 0]
            idx[i] = indices
    
    imp_data.index = org_data.index
    spots = np.array(org_data.index)        
    fig = plt.figure(tight_layout=True, figsize=fig_size)
    plt.subplots_adjust(wspace= 0.25, hspace= 0.25)
    for i in range(len(geneSet)):
            sub=fig.add_subplot(1,4,i+1)
            org = org_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
            imp = imp_data.loc[spots[idx[geneSet[i]]],geneSet[i]]
            sub.scatter(org, imp, c=PointColor)
            sub.plot(org, org, c=LineColor)
            sub.tick_params(axis='y', labelsize=15)
            sub.tick_params(axis='x', labelsize=15)
            sub.set_xlabel('original', fontsize=15)
            sub.set_ylabel('imputed', fontsize=15)
            sub.set_title(geneSet[i], fontsize=20)
    
    fig.savefig(fig_path+'linearity_pl_{}_{}.tiff'.format(dName,mName), dpi=600, format="tiff")

    return(fig)








    







