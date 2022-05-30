# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:48:36 2022

@author: Gulben AVSAR
"""
from scripts import loadData_forEval as ppD
from scripts import eval_metrics as evm
from scripts import figure_scripts as figs
from scripts import table_scripts as tabs
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

def reprodct_results():
    # define the path for figures
    fig_path = './figures/'
    fpath = Path(fig_path)
    fpath.mkdir(parents=True, exist_ok=True)
    
    # define the path for evaluation files
    eval_path = './evaluations/'
    epath = Path(eval_path)
    epath.mkdir(parents=True, exist_ok=True)
    
    # upload the preprocessed datasets
    spa_org_A, spa_imp_A, _ = ppD.load_datasets('PDAC-A','spage')
    spa_org_B, spa_imp_B, _ = ppD.load_datasets('PDAC-B','spage')
    spa_org_L, spa_imp_L, _ = ppD.load_datasets('Liver','spage')
    
    stp_org_A, stp_imp_A, _ = ppD.load_datasets('PDAC-A','stplus')
    stp_org_B, stp_imp_B, _ = ppD.load_datasets('PDAC-B','stplus')
    stp_org_L, stp_imp_L, _ = ppD.load_datasets('Liver','stplus')
    
    gim_org_A, gim_imp_A, _ = ppD.load_datasets('PDAC-A','gimvi')
    gim_org_B, gim_imp_B, _ = ppD.load_datasets('PDAC-B','gimvi')
    gim_org_L, gim_imp_L, _ = ppD.load_datasets('Liver','gimvi')
    
    tan_org_A, tan_imp_A, _ = ppD.load_datasets('PDAC-A','tangram')
    tan_org_B, tan_imp_B, _ = ppD.load_datasets('PDAC-B','tangram')
    tan_org_L, tan_imp_L, _ = ppD.load_datasets('Liver','tangram')
    
    st5_org_A, st5_imp_A, mIdx5_A = ppD.load_datasets('PDAC-A','5stlearn')
    st5_org_B, st5_imp_B, mIdx5_B = ppD.load_datasets('PDAC-B','5stlearn')
    st5_org_L, st5_imp_L, mIdx5_L = ppD.load_datasets('Liver','5stlearn')
    
    st30_org_A, st30_imp_A, mIdx30_A = ppD.load_datasets('PDAC-A','30stlearn')
    st30_org_B, st30_imp_B, mIdx30_B = ppD.load_datasets('PDAC-B','30stlearn')
    st30_org_L, st30_imp_L, mIdx30_L = ppD.load_datasets('Liver','30stlearn')
    
    
    ########################################################################
    ###### Evaluations gene-wise
    mthds = ['spage', 'stplus', 'gimvi', 'tangram','5stlearn','30stlearn']
    # for PDAC-A
    orgs = [spa_org_A,stp_org_A,gim_org_A,tan_org_A,st5_org_A,st30_org_A]
    imps = [spa_imp_A,stp_imp_A,gim_imp_A,tan_imp_A,st5_imp_A,st30_imp_A]
    mIdxs = ['nan','nan','nan','nan',mIdx5_A,mIdx30_A]
    for i in range(len(mthds)):
        gw = evm.eval_GW(orgs[i], imps[i], 'pdacA', mthds[i], mIdxs[i],eval_path)
    
    # for PDAC-B
    orgs = [spa_org_B,stp_org_B,gim_org_B,tan_org_B,st5_org_B,st30_org_B]
    imps = [spa_imp_B,stp_imp_B,gim_imp_B,tan_imp_B,st5_imp_B,st30_imp_B]
    mIdxs = ['nan','nan','nan','nan',mIdx5_B,mIdx30_B]
    for i in range(len(mthds)):
        gw = evm.eval_GW(orgs[i], imps[i], 'pdacB', mthds[i], mIdxs[i],eval_path)
    
    # for Liver
    orgs = [spa_org_L,stp_org_L,gim_org_L,tan_org_L,st5_org_L,st30_org_L]
    imps = [spa_imp_L,stp_imp_L,gim_imp_L,tan_imp_L,st5_imp_L,st30_imp_L]
    mIdxs = ['nan','nan','nan','nan',mIdx5_L,mIdx30_L]
    for i in range(len(mthds)):
        gw = evm.eval_GW(orgs[i], imps[i], 'Liver', mthds[i], mIdxs[i],eval_path)
    
    ###### Evaluations spot-wise
    mthds = ['spage', 'stplus', 'gimvi', 'tangram','5stlearn','30stlearn']
    # for PDAC-A
    orgs = [spa_org_A,stp_org_A,gim_org_A,tan_org_A,st5_org_A,st30_org_A]
    imps = [spa_imp_A,stp_imp_A,gim_imp_A,tan_imp_A,st5_imp_A,st30_imp_A]
    mIdxs = ['nan','nan','nan','nan',mIdx5_A,mIdx30_A]
    for i in range(len(mthds)):
        sw = evm.eval_SW(orgs[i], imps[i], 'pdacA', mthds[i], mIdxs[i],eval_path)
    
    # for PDAC-B
    orgs = [spa_org_B,stp_org_B,gim_org_B,tan_org_B,st5_org_B,st30_org_B]
    imps = [spa_imp_B,stp_imp_B,gim_imp_B,tan_imp_B,st5_imp_B,st30_imp_B]
    mIdxs = ['nan','nan','nan','nan',mIdx5_B,mIdx30_B]
    for i in range(len(mthds)):
        sw = evm.eval_SW(orgs[i], imps[i], 'pdacB', mthds[i], mIdxs[i],eval_path)
    
    # for Liver
    orgs = [spa_org_L,stp_org_L,gim_org_L,tan_org_L,st5_org_L,st30_org_L]
    imps = [spa_imp_L,stp_imp_L,gim_imp_L,tan_imp_L,st5_imp_L,st30_imp_L]
    mIdxs = ['nan','nan','nan','nan',mIdx5_L,mIdx30_L]
    for i in range(len(mthds)):
        sw = evm.eval_SW(orgs[i], imps[i], 'Liver', mthds[i], mIdxs[i],eval_path)
    
    
    ########################################################################
    #################################### Figure 1.a
    figAa = figs.Fig1a(fig_path)
    
    
    # #################################### Figure 1.b
    # table in the figure
    orgs = [[spa_org_A,stp_org_A,gim_org_A,tan_org_A,st5_org_A,st30_org_A],
            [spa_org_B,stp_org_B,gim_org_B,tan_org_B,st5_org_B,st30_org_B],
            [spa_org_L,stp_org_L,gim_org_L,tan_org_L,st5_org_L,st30_org_L]]
    imps = [[spa_imp_A,stp_imp_A,gim_imp_A,tan_imp_A,st5_imp_A,st30_imp_A],
            [spa_imp_B,stp_imp_B,gim_imp_B,tan_imp_B,st5_imp_B,st30_imp_B],
            [spa_imp_L,stp_imp_L,gim_imp_L,tan_imp_L,st5_imp_L,st30_imp_L]]
    mIdxs = [['nan','nan','nan','nan',mIdx5_A,mIdx30_A],
             ['nan','nan','nan','nan',mIdx5_B,mIdx30_B],
             ['nan','nan','nan','nan',mIdx5_L,mIdx30_L]]
    
    figAb = tabs.Fig1b(orgs, imps, mIdxs, fig_path)
    
    ########################################################################
    #################################### Figure 2
    figB = figs.Fig2(fig_path)
    
    #################################### Figure 3 (parts)
    mthds = ['spage', 'stplus', 'gimvi', 'tangram','5stlearn','30stlearn']
    # for PDAC-A
    orgs = [spa_org_A,stp_org_A,gim_org_A,tan_org_A,st5_org_A,st30_org_A]
    imps = [spa_imp_A,stp_imp_A,gim_imp_A,tan_imp_A,st5_imp_A,st30_imp_A]
    mIdxs = ['nan','nan','nan','nan',mIdx5_A,mIdx30_A]
    for i in range(len(mthds)):
        figC = figs.FigC(orgs[i], imps[i], mIdxs[i], 'PDAC-A', mthds[i],
                         (18,4), 1, 0.5, fig_path)
    
    # for PDAC-B
    orgs = [spa_org_B,stp_org_B,gim_org_B,tan_org_B,st5_org_B,st30_org_B]
    imps = [spa_imp_B,stp_imp_B,gim_imp_B,tan_imp_B,st5_imp_B,st30_imp_B]
    mIdxs = ['nan','nan','nan','nan',mIdx5_B,mIdx30_B]
    for i in range(len(mthds)):
        figC = figs.FigC(orgs[i], imps[i], mIdxs[i], 'PDAC-B', mthds[i],
                         (18,4), 1, 0.5, fig_path)
    
    # for Liver
    orgs = [spa_org_L,stp_org_L,gim_org_L,tan_org_L,st5_org_L,st30_org_L]
    imps = [spa_imp_L,stp_imp_L,gim_imp_L,tan_imp_L,st5_imp_L,st30_imp_L]
    mIdxs = ['nan','nan','nan','nan',mIdx5_L,mIdx30_L]
    for i in range(len(mthds)):
        figC = figs.FigC(orgs[i], imps[i], mIdxs[i], 'Liver', mthds[i],
                         (18,4), 1, 0.5, fig_path)
    
    
    #################################### Figure 4 (parts)
    mthds = ['stplus', 'gimvi']
    # for PDAC-A
    orgs = [stp_org_A,gim_org_A]
    imps = [stp_imp_A,gim_imp_A]
    mIdxs = ['nan','nan']
    for i in range(len(mthds)):
        figD = figs.FigD(orgs[i], imps[i], mIdxs[i], 'PDAC-A', mthds[i],
                         (18,4), 'brown', 'dimgrey', fig_path)
    
    # for PDAC-A
    orgs = [stp_org_B,gim_org_B]
    imps = [stp_imp_B,gim_imp_B]
    mIdxs = ['nan','nan']
    for i in range(len(mthds)):
        figD = figs.FigD(orgs[i], imps[i], mIdxs[i], 'PDAC-B', mthds[i],
                         (18,4), 'brown', 'dimgrey', fig_path)
    
    # for Liver
    orgs = [stp_org_L,gim_org_L]
    imps = [stp_imp_L,gim_imp_L]
    mIdxs = ['nan','nan']
    for i in range(len(mthds)):
        figD = figs.FigD(orgs[i], imps[i], mIdxs[i], 'Liver', mthds[i],
                         (18,4), 'brown', 'dimgrey', fig_path)
    
    
    #################################### Suppl.Figure 5 (parts)
    mthds = ['stplus','stplus','gimvi','gimvi']
    # for Liver
    orgs = [stp_org_L,stp_org_L,gim_org_L,gim_org_L]
    imps = [stp_imp_L,stp_imp_L,gim_imp_L,gim_imp_L]
    mIdxs = ['nan','nan','nan','nan']
    lims = [[1,2], [0.5,1.3], [0,0.0004], [0, 0.00012]]
    genes = ['PDIA3', 'EPN1', 'PDIA3', 'EPN1']
    for i in range(len(mthds)):
        sfigC = figs.SupplFigC(orgs[i], imps[i], mIdxs[i], genes[i], 'Liver',
                               mthds[i], lims[i], (8,6), 'brown', 'dimgrey',
                               1, fig_path)
    
    
    #################################### Figure 5 (parts)
    figE = figs.FigE(fig_path)
    
    
    #################################### Table 1
    tabA = tabs.TabA(fig_path)
    
    
    #################################### Table 2
    orgs = [[spa_org_A,stp_org_A,gim_org_A,tan_org_A,st5_org_A,st30_org_A],
            [spa_org_B,stp_org_B,gim_org_B,tan_org_B,st5_org_B,st30_org_B],
            [spa_org_L,stp_org_L,gim_org_L,tan_org_L,st5_org_L,st30_org_L]]
    imps = [[spa_imp_A,stp_imp_A,gim_imp_A,tan_imp_A,st5_imp_A,st30_imp_A],
            [spa_imp_B,stp_imp_B,gim_imp_B,tan_imp_B,st5_imp_B,st30_imp_B],
            [spa_imp_L,stp_imp_L,gim_imp_L,tan_imp_L,st5_imp_L,st30_imp_L]]
    mIdxs = [['nan','nan','nan','nan',mIdx5_A,mIdx30_A],
             ['nan','nan','nan','nan',mIdx5_B,mIdx30_B],
             ['nan','nan','nan','nan',mIdx5_L,mIdx30_L]]
    
    tabB = tabs.TabB(orgs, imps, mIdxs, fig_path)
    
    
    ########################################################################
    #################################### Extra in text
    ####### find percent HIGHER THAN cutOff in a metric for a method
    100 - figs.below_metric('PCC', 'gimvi', 0.0)
    100 - figs.below_metric('CS', 'gimvi', 0.7)
    100 - figs.below_metric('CS', 'stplus', 0.7)
    
    ####### find percent LOWER THAN cutOff in a metric for a method
    figs.below_metric('RMSLE', 'gimvi', 0.01)
    figs.below_metric('RMSLE', 'tangram', 0.01)
    
    
    #######(predictions with only at most 3% error rate)
    mthds = ['stplus', 'gimvi']
    # for Liver
    orgs = [stp_org_L,gim_org_L]
    imps = [stp_imp_L,gim_imp_L]
    
    tabs.preds_w_Eratio(orgs, imps, 0.03, ['PDIA3'])
    tabs.preds_w_Eratio(orgs, imps, 0.03, ['EPN1'])
    
    
    ########################################################################
    ################################ Evaluations with synthetic SRT dataset
    spa_org_S, spa_imp_S, _ = ppD.load_datasets('SimulatedST','spage')
    stp_org_S, stp_imp_S, _ = ppD.load_datasets('SimulatedST','stplus')
    gim_org_S, gim_imp_S, _ = ppD.load_datasets('SimulatedST','gimvi')
    tan_org_S, tan_imp_S, _ = ppD.load_datasets('SimulatedST','tangram')
    
    ###### Evaluations gene-wise
    mthds = ['spage', 'stplus', 'gimvi', 'tangram']
    # for SimulatedST
    orgs = [spa_org_S,stp_org_S,gim_org_S,tan_org_S]
    imps = [spa_imp_S,stp_imp_S,gim_imp_S,tan_imp_S]
    mIdxs = ['nan','nan','nan','nan']
    for i in range(len(mthds)):
        gw = evm.eval_GW(orgs[i], imps[i], 'SimulatedST', mthds[i], mIdxs[i], eval_path)
    
    
    ###### Evaluations spot-wise
    mthds = ['spage', 'stplus', 'gimvi', 'tangram']
    
    orgs = [spa_org_S,stp_org_S,gim_org_S,tan_org_S]
    imps = [spa_imp_S,stp_imp_S,gim_imp_S,tan_imp_S]
    mIdxs = ['nan','nan','nan','nan']
    for i in range(len(mthds)):
        sw = evm.eval_SW(orgs[i], imps[i], 'SimulatedST', mthds[i], mIdxs[i], eval_path)
    
    
    ###### Figure 6
    fgw, fsw = figs.FigF(fig_path)
    
    #######(predictions with only at most 3% error rate)
    mthds = ['stplus', 'gimvi']
    # for Liver
    orgs = [stp_org_S,gim_org_S]
    imps = [stp_imp_S,gim_imp_S]
    
    tabs.preds_w_Eratio(orgs, imps, 0.03, [])
    
    return()
