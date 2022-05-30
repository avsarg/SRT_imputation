# -*- coding: utf-8 -*-
"""
Created on Tue May 24 12:46:22 2022

@author: Gulben AVSAR

DEG analysis with SpatialDE

doi: https://doi.org/10.1038/nmeth.4636
https://github.com/Teichlab/SpatialDE
"""
######################################################################################
import os
os.chdir("E:\\\PhD_Thesis")
os.getcwd()

from scripts import loadData_forEval as ppD
from scripts import run_SpatialDE as spade
import matplotlib.pyplot as plt
from gprofiler import GProfiler
import numpy as np

######################################## Load datasets
stp_org_L, stp_imp_L, _ = ppD.load_datasets('Liver','stplus')
gim_org_L, gim_imp_L, _ = ppD.load_datasets('Liver','gimvi')

gim_org_L.index = stp_org_L.index
gim_imp_L.index = stp_org_L.index

######################################## stPlus
############# for original dataset
res_stp_org_L = spade.SpatialDE_run(stp_org_L, 'Liver')
sig_stp_org_L = spade.pull_signf_DEGs(res_stp_org_L, 'pval')

# export the result
res_stp_org_L.to_csv('./SpatialDE_stplus_org_L.txt')

#  visualize the spatially DE genes
ind = sig_stp_org_L['g'].index[0:3]
g = sig_stp_org_L.loc[ind,'g'].to_list()

fig = spade.visualise_genes(stp_org_L, 'Liver', g, 10, (18,5))

############# for imputed dataset
res_stp_imp_L = spade.SpatialDE_run(stp_imp_L, 'Liver')
sig_stp_imp_L = spade.pull_signf_DEGs(res_stp_imp_L, 'pval')

# export the result
res_stp_imp_L.to_csv('./SpatialDE_stplus_imp_L.txt')

#  visualize the spatially DE genes
ind = sig_stp_imp_L['g'].index[0:3]
g = sig_stp_imp_L.loc[ind,'g'].to_list()

fig = spade.visualise_genes(stp_imp_L, 'Liver', g, 10, (18,5))


######################################## gimVI
############# for original dataset
gim_org_L.isnull().values.any()
res_gim_org_L = spade.SpatialDE_run(gim_org_L, 'Liver')
sig_gim_org_L = spade.pull_signf_DEGs(res_gim_org_L, 'pval')

# export the result
res_gim_org_L.to_csv('./SpatialDE_gimvi_org_L.txt')

#  visualize the spatially DE genes
ind = sig_gim_org_L['g'].index[0:3]
g = sig_gim_org_L.loc[ind,'g'].to_list()

fig = spade.visualise_genes(gim_org_L, 'Liver', g, 10, (18,5))

############# for imputed dataset
res_gim_imp_L = spade.SpatialDE_run(gim_imp_L, 'Liver')
sig_gim_imp_L = spade.pull_signf_DEGs(res_gim_imp_L, 'pval')

# export the result
res_gim_imp_L.to_csv('./SpatialDE_gimvi_imp_L.txt')

#  visualize the spatially DE genes
ind = sig_gim_imp_L['g'].index[0:3]
g = sig_gim_imp_L.loc[ind,'g'].to_list()

fig = spade.visualise_genes(gim_imp_L, 'Liver', g, 10, (18,5))




######################################## plot the number of significant DEGs
sizes = [len(sig_stp_org_L), len(sig_stp_imp_L),
         len(sig_gim_org_L), len(sig_gim_imp_L)]
labels= ['stplus original', 'stplus predicted',
         'gimvi original', 'gimvi predicted']
fig, ax = plt.subplots(figsize =(12, 9))
# Horizontal Bar Plot
# c = ['#e894a7', '#c06378','#8b4e5e','#663845']
ax.barh(labels, sizes, color = 'mediumseagreen')
# Remove axes splines
for s in ['top', 'right']:
    ax.spines[s].set_visible(False)
# Add padding between axes and labels
ax.xaxis.set_tick_params(pad = 5, labelsize=18)
ax.yaxis.set_tick_params(pad = 10, labelsize=18)
# Show top values
ax.invert_yaxis()
# Add annotation to bars
for i in ax.patches:
    plt.text(i.get_width()+5.0, i.get_y()+0.5,
                  str(round((i.get_width()), 2)),
                  fontsize = 18, fontweight ='bold',
                  color ='black')
# Add Plot Title
ax.set_title('# of spatially DEGs', loc ='left', fontsize = 20)
fig.tight_layout()
fig.savefig('./figures/NumbOfSpatialDEGs.tiff', dpi=600, format="tiff")
# Show Plot
plt.show()



######################################## GO analysis
gp = GProfiler(return_dataframe=True)

# stplus
GO_stp_org = gp.profile(organism='hsapiens',
                        query=res_stp_org_L.query('qval < 0.05')['g'].sort_values().to_list(),
                        ordered=True)
GOBP_stp_org = GO_stp_org[GO_stp_org.source == 'GO:BP']

GO_stp_imp = gp.profile(organism='hsapiens',
                        query=res_stp_imp_L.query('qval < 0.05')['g'].sort_values().to_list(),
                        ordered=True)
GOBP_stp_imp = GO_stp_imp[GO_stp_imp.source == 'GO:BP']


# gimvi
GO_gim_org = gp.profile(organism='hsapiens',
                        query=res_gim_org_L.query('qval < 0.05')['g'].sort_values().to_list(),
                        ordered=True)
GOBP_gim_org = GO_gim_org[GO_gim_org.source == 'GO:BP']

GO_gim_imp = gp.profile(organism='hsapiens',
                        query=res_gim_imp_L.query('qval < 0.05')['g'].sort_values().to_list(),
                        ordered=True)
GOBP_gim_imp = GO_gim_imp[GO_gim_imp.source == 'GO:BP']



# GO for common DEGs
shrd_DEGs = list(np.intersect1d(np.intersect1d(sig_gim_org_L.g, sig_stp_org_L.g),
               np.intersect1d(sig_gim_imp_L.g, sig_stp_imp_L.g)))
gp.profile(organism='hsapiens', query=shrd_DEGs)
gp.profile(organism='hsapiens', query=shrd_DEGs)['name']



GO_all = gp.profile(organism='hsapiens', query=stp_org_L.columns.to_list(), ordered=True)
GOBP_all = GO_all[GO_all.source == 'GO:BP']



######################################## Visualise GO analysis
def go_term_plot(df, numb, color):
    # Example data
    trms = df['terms'][0:numb].to_list()
    y_pos = np.arange(len(trms))
    pvals = df['adj_p'][0:numb].to_list()
    pvals=[float(x) for x in pvals]
    
    # # error = np.random.rand(len(people))
    # fig, ax = plt.subplots()
    # ax.barh(y_pos, -np.log10(pvals), align='center', color=color)
    # ax.set_yticks(y_pos)
    # ax.set_yticklabels(trms)
    # ax.invert_yaxis()  # labels read top-to-bottom
    # ax.set_xlabel('- log10(p)')
    
    return(trms, y_pos, pvals)
    

GOBP_list = [GOBP_stp_org, GOBP_stp_imp, GOBP_gim_org, GOBP_gim_imp, GOBP_all]
title_list = ['stplus_original', 'stplus_predicted',
              'gimvi_original ', 'gimvi_predicted',
              'all 400 genes', '% of same terms']
fig, ax = plt.subplots(3,2, figsize=(20,8))
ax = ax.ravel()
for i in range(len(GOBP_list)+1):
    if i == len(GOBP_list):
        labels = ['stplus_predicted','gimvi_predicted']
        x = len(set(GOBP_all['native'].to_list()) & set(GOBP_stp_imp['native'].to_list()))
        y = len(set(GOBP_all['native'].to_list()) & set(GOBP_gim_imp['native'].to_list()))
        sizes = [x/len(GOBP_all['native']), y/len(GOBP_all['native'])]
    else:
        df = GOBP_list[i]
        labels = df.name[:5].to_list()
        # sizes = list((df.intersection_size[:4] / df.query_size[:4]).values)
        sizes = list(-np.log10(df.p_value[:5]).values)
    ax[i].barh(labels, sizes, align='center', color='lightseagreen')
    ax[i].set_yticks(labels)
    ax[i].set_yticklabels(labels, fontsize=20)
    ax[i].invert_yaxis()  # labels read top-to-bottom
    
    if i == len(GOBP_list):
        ax[i].set_xlabel('% of hit', fontsize=20)
    else:
        ax[i].set_xlabel('-log10(pval)', fontsize=20)
    ax[i].set_title(title_list[i], fontsize=20, fontweight='bold')
fig.tight_layout()
fig.subplots_adjust(top=1.5)
fig.savefig('./figures/GOterms.tiff',
            dpi=300, format="tiff", bbox_inches="tight")














