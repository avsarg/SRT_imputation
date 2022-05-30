# -*- coding: utf-8 -*-
"""
@author: Gulben AVSAR
"""
import os
os.chdir("E:\\\PhD_Thesis")
os.getcwd()

from Codes_Imputation_Manuscript.scripts import loadData_forEval as ppD
from Codes_Imputation_Manuscript.scripts import eval_metrics as evm
from Codes_Imputation_Manuscript.scripts import figure_scripts as figs
from Codes_Imputation_Manuscript.scripts import table_scripts as tabs
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# define the path for figures
fig_path = './Codes_Imputation_Manuscript/figures/'
fpath = Path(fig_path)
fpath.mkdir(parents=True, exist_ok=True)


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
    gw = evm.eval_GW(orgs[i], imps[i], 'pdacA', mthds[i], mIdxs[i])

# for PDAC-B
orgs = [spa_org_B,stp_org_B,gim_org_B,tan_org_B,st5_org_B,st30_org_B]
imps = [spa_imp_B,stp_imp_B,gim_imp_B,tan_imp_B,st5_imp_B,st30_imp_B]
mIdxs = ['nan','nan','nan','nan',mIdx5_B,mIdx30_B]
for i in range(len(mthds)):
    gw = evm.eval_GW(orgs[i], imps[i], 'pdacB', mthds[i], mIdxs[i])

# for Liver
orgs = [spa_org_L,stp_org_L,gim_org_L,tan_org_L,st5_org_L,st30_org_L]
imps = [spa_imp_L,stp_imp_L,gim_imp_L,tan_imp_L,st5_imp_L,st30_imp_L]
mIdxs = ['nan','nan','nan','nan',mIdx5_L,mIdx30_L]
for i in range(len(mthds)):
    gw = evm.eval_GW(orgs[i], imps[i], 'Liver', mthds[i], mIdxs[i])

###### Evaluations spot-wise
mthds = ['spage', 'stplus', 'gimvi', 'tangram','5stlearn','30stlearn']
# for PDAC-A
orgs = [spa_org_A,stp_org_A,gim_org_A,tan_org_A,st5_org_A,st30_org_A]
imps = [spa_imp_A,stp_imp_A,gim_imp_A,tan_imp_A,st5_imp_A,st30_imp_A]
mIdxs = ['nan','nan','nan','nan',mIdx5_A,mIdx30_A]
for i in range(len(mthds)):
    sw = evm.eval_SW(orgs[i], imps[i], 'pdacA', mthds[i], mIdxs[i])

# for PDAC-B
orgs = [spa_org_B,stp_org_B,gim_org_B,tan_org_B,st5_org_B,st30_org_B]
imps = [spa_imp_B,stp_imp_B,gim_imp_B,tan_imp_B,st5_imp_B,st30_imp_B]
mIdxs = ['nan','nan','nan','nan',mIdx5_B,mIdx30_B]
for i in range(len(mthds)):
    sw = evm.eval_SW(orgs[i], imps[i], 'pdacB', mthds[i], mIdxs[i])

# for Liver
orgs = [spa_org_L,stp_org_L,gim_org_L,tan_org_L,st5_org_L,st30_org_L]
imps = [spa_imp_L,stp_imp_L,gim_imp_L,tan_imp_L,st5_imp_L,st30_imp_L]
mIdxs = ['nan','nan','nan','nan',mIdx5_L,mIdx30_L]
for i in range(len(mthds)):
    sw = evm.eval_SW(orgs[i], imps[i], 'Liver', mthds[i], mIdxs[i])


###### Gene-wise CS with p-values
mthds = ['spage', 'stplus', 'gimvi', 'tangram','5stlearn','30stlearn']
# for Liver
orgs = [spa_org_A,stp_org_A,gim_org_A,tan_org_A,st5_org_A,st30_org_A]
imps = [spa_imp_A,stp_imp_A,gim_imp_A,tan_imp_A,st5_imp_A,st30_imp_A]
mIdxs = ['nan','nan','nan','nan',mIdx5_A,mIdx30_A]
# for i in range(len(mthds)):
#     csp = evm.CS_pvals(orgs[i], imps[i], 'Liver', mthds[i], mIdxs[i])
i=2 #for gimVI
csp = evm.CS_pvals(orgs[i], imps[i], 'Liver', mthds[i], mIdxs[i])


###### Gene-wise PCC with p-values
mthds = ['spage', 'stplus', 'gimvi', 'tangram','5stlearn','30stlearn']
# for Liver
orgs = [spa_org_A,stp_org_A,gim_org_A,tan_org_A,st5_org_A,st30_org_A]
imps = [spa_imp_A,stp_imp_A,gim_imp_A,tan_imp_A,st5_imp_A,st30_imp_A]
mIdxs = ['nan','nan','nan','nan',mIdx5_A,mIdx30_A]
# for i in range(len(mthds)):
#     pccp = evm.PPC_pvals(orgs[i], imps[i], 'Liver', mthds[i], mIdxs[i])
i=2 #for gimVI
pccp = evm.PCC_pvals(orgs[i], imps[i], 'Liver', mthds[i], mIdxs[i])

########################################################################
#################################### Figure A.a
figAa = figs.FigAa()


# #################################### Figure A.b
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

figAb = tabs.FigAb(orgs, imps, mIdxs)

########################################################################
#################################### Figure B
figB = figs.FigB()

#################################### Figure C (parts)
mthds = ['spage', 'stplus', 'gimvi', 'tangram','5stlearn','30stlearn']
# for PDAC-A
orgs = [spa_org_A,stp_org_A,gim_org_A,tan_org_A,st5_org_A,st30_org_A]
imps = [spa_imp_A,stp_imp_A,gim_imp_A,tan_imp_A,st5_imp_A,st30_imp_A]
mIdxs = ['nan','nan','nan','nan',mIdx5_A,mIdx30_A]
for i in range(len(mthds)):
    figC = figs.FigC(orgs[i], imps[i], mIdxs[i], 'PDAC-A', mthds[i], (18,4), 1, 0.5)

# for PDAC-B
orgs = [spa_org_B,stp_org_B,gim_org_B,tan_org_B,st5_org_B,st30_org_B]
imps = [spa_imp_B,stp_imp_B,gim_imp_B,tan_imp_B,st5_imp_B,st30_imp_B]
mIdxs = ['nan','nan','nan','nan',mIdx5_B,mIdx30_B]
for i in range(len(mthds)):
    figC = figs.FigC(orgs[i], imps[i], mIdxs[i], 'PDAC-B', mthds[i], (18,4), 1, 0.5)

# for Liver
orgs = [spa_org_L,stp_org_L,gim_org_L,tan_org_L,st5_org_L,st30_org_L]
imps = [spa_imp_L,stp_imp_L,gim_imp_L,tan_imp_L,st5_imp_L,st30_imp_L]
mIdxs = ['nan','nan','nan','nan',mIdx5_L,mIdx30_L]
for i in range(len(mthds)):
    figC = figs.FigC(orgs[i], imps[i], mIdxs[i], 'Liver', mthds[i], (18,4), 1, 0.5)


#################################### Figure D (parts)
mthds = ['stplus', 'gimvi']
# for PDAC-A
orgs = [stp_org_A,gim_org_A]
imps = [stp_imp_A,gim_imp_A]
mIdxs = ['nan','nan']
for i in range(len(mthds)):
    figD = figs.FigD(orgs[i], imps[i], mIdxs[i], 'PDAC-A', mthds[i], (18,4), 'brown', 'dimgrey')

# for PDAC-A
orgs = [stp_org_B,gim_org_B]
imps = [stp_imp_B,gim_imp_B]
mIdxs = ['nan','nan']
for i in range(len(mthds)):
    figD = figs.FigD(orgs[i], imps[i], mIdxs[i], 'PDAC-B', mthds[i], (18,4), 'brown', 'dimgrey')

# for Liver
orgs = [stp_org_L,gim_org_L]
imps = [stp_imp_L,gim_imp_L]
mIdxs = ['nan','nan']
for i in range(len(mthds)):
    figD = figs.FigD(orgs[i], imps[i], mIdxs[i], 'Liver', mthds[i], (18,4), 'brown', 'dimgrey')


#################################### Suppl.Figure C (parts)
mthds = ['stplus','stplus','gimvi','gimvi']
# for Liver
orgs = [stp_org_L,stp_org_L,gim_org_L,gim_org_L]
imps = [stp_imp_L,stp_imp_L,gim_imp_L,gim_imp_L]
mIdxs = ['nan','nan','nan','nan']
lims = [[1,2], [0.5,1.3], [0,0.0004], [0, 0.00012]]
genes = ['PDIA3', 'EPN1', 'PDIA3', 'EPN1']
for i in range(len(mthds)):
    sfigC = figs.SupplFigC(orgs[i], imps[i], mIdxs[i], genes[i], 'Liver',
                           mthds[i], lims[i], (8,6), 'brown', 'dimgrey', 1)


#################################### Figure E (parts)
figE = figs.FigE()


#################################### Suppl.Figure D
sfigD = figs.SupplFigD('CS',0.7,'cos_similarity')
# sfigD = figs.SupplFigD('PCC',0.5, 'corr_val')

#################################### Table A
tabA = tabs.TabA()


#################################### Table A
orgs = [[spa_org_A,stp_org_A,gim_org_A,tan_org_A,st5_org_A,st30_org_A],
        [spa_org_B,stp_org_B,gim_org_B,tan_org_B,st5_org_B,st30_org_B],
        [spa_org_L,stp_org_L,gim_org_L,tan_org_L,st5_org_L,st30_org_L]]
imps = [[spa_imp_A,stp_imp_A,gim_imp_A,tan_imp_A,st5_imp_A,st30_imp_A],
        [spa_imp_B,stp_imp_B,gim_imp_B,tan_imp_B,st5_imp_B,st30_imp_B],
        [spa_imp_L,stp_imp_L,gim_imp_L,tan_imp_L,st5_imp_L,st30_imp_L]]
mIdxs = [['nan','nan','nan','nan',mIdx5_A,mIdx30_A],
         ['nan','nan','nan','nan',mIdx5_B,mIdx30_B],
         ['nan','nan','nan','nan',mIdx5_L,mIdx30_L]]

tabB = tabs.TabB(orgs, imps, mIdxs)



########################################################################
#################################### Extra in text
####### find percent HIGHER THAN cutOff in a metric for a method
100 - figs.below_metric('PCC', 'gimvi', 0.0)
100 - figs.below_metric('CS', 'gimvi', 0.7)
100 - figs.below_metric('CS', 'stplus', 0.7)

####### find percent LOWER THAN cutOff in a metric for a method
figs.below_metric('RMSLE', 'gimvi', 1)
figs.below_metric('RMSLE', 'tangram', 1)
figs.below_metric('RMSLE', 'stplus', 1)
figs.below_metric('RMSLE', 'spage', 1)

#######(predictions with only at most 3% error rate)
mthds = ['stplus', 'gimvi']
# for Liver
orgs = [stp_org_L,gim_org_L]
imps = [stp_imp_L,gim_imp_L]

tabs.preds_w_Eratio(orgs, imps, 0.03, ['PDIA3'])
tabs.preds_w_Eratio(orgs, imps, 0.03, ['EPN1'])

tabs.preds_w_Eratio(orgs, imps, 0.03, [])


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
    gw = evm.eval_GW(orgs[i], imps[i], 'SimulatedST', mthds[i], mIdxs[i])


###### Evaluations spot-wise
mthds = ['spage', 'stplus', 'gimvi', 'tangram']

orgs = [spa_org_S,stp_org_S,gim_org_S,tan_org_S]
imps = [spa_imp_S,stp_imp_S,gim_imp_S,tan_imp_S]
mIdxs = ['nan','nan','nan','nan']
for i in range(len(mthds)):
    sw = evm.eval_SW(orgs[i], imps[i], 'SimulatedST', mthds[i], mIdxs[i])


###### Figure F
fgw, fsw = figs.FigF()



#######(predictions with only at most 3% error rate)
mthds = ['stplus', 'gimvi']
# for Liver
orgs = [stp_org_S,gim_org_S]
imps = [stp_imp_S,gim_imp_S]

tabs.preds_w_Eratio(orgs, imps, 0.03, [])











