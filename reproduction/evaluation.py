"""
@author: Gulben AVSAR
Perform the qualitative and quantitative evaluation for the predictions.

Reference
-------
    Avsar G.,
"""
# from Codes_Imputation_Manuscript.scripts import loadData_forEval as ppD
from scripts import eval_metrics as evm
from scripts import figure_scripts as figs
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def QQeval(org_data, imp_data, g, dName, mName=None, pl_g=None):
    #org_data: original dataset
    #imp_data: imputed dataset
    #g: List of genes for evaluations
    #dName: name for dataset (will be used in saved file names)
    #mName: name for used method (will be used in saved file names)
    #pl_g: name of 4 genes to plot (default: will be chosen randomly)
    
    assert(type(g) != 'list'), "genes should be given in a list"
    assert(type(org_data) != 'pandas.core.frame.DataFrame'), "dataset should be given in a pandas df"
    assert(type(imp_data) != 'pandas.core.frame.DataFrame'), "dataset should be given in a pandas df"
    assert(type(dName) != 'str'), "dataset name should be given in a string"
    assert(type(mName) != 'str'), "dataset name should be given in a string"
    assert(len(pl_g) <= 4), "give at most 4 genes to plot"
    
    # define the path for figures
    fig_path = './QQeval/figures/'
    fpath = Path(fig_path)
    fpath.mkdir(parents=True, exist_ok=True)
    
    # define the path for evaluation files
    eval_path = './QQeval/evaluations/'
    epath = Path(eval_path)
    epath.mkdir(parents=True, exist_ok=True)
    
    org = org_data.copy()
    imp = imp_data.copy()
    
        
    mIdx = 'nan'
    
    ###### Evaluations gene-wise
    gw = evm.eval_GW(org, imp, dName, mName, mIdx, eval_path)
    
    ###### Evaluations spot-wise
    sw = evm.eval_SW(org, imp, dName, mName, mIdx, eval_path)
    
    ######
    figs.comparison_pl(org, imp, g, dName, fig_path,
                      mName=mName, pl_g=pl_g, fig_size=(18,4),
                      a_org=1, a_imp=0.5)
    
    ######
    figs.linearity_pl(org, imp, g, dName, fig_path,
                      mName=mName, pl_g=pl_g, fig_size=(18,4),
                      PointColor='brown', LineColor='dimgrey')

    return()

import pandas as pd
def QQtest():
    org = pd.read_csv('./reproduction/predictions/gimvi/Orig_gimvi_Liver.txt.gz', index_col=0)
    imp = pd.read_csv('./reproduction/predictions/gimvi/Imp_gimvi_Liver.txt.gz', index_col=0)
    genes = org.columns.to_list()
    return(org, imp, genes)
