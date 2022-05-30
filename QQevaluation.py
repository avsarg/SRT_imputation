"""
@author: Gulben AVSAR
Perform the qualitative and quantitative evaluation for the predictions.

Reference
-------
    [1] Avsar G.,
"""

from reproduction import evaluation as evals

def QQevaluation(org_data, imp_data, g, dName, mName=None, pl_g=None):
    if mName is QQevaluation.__defaults__[0]:
        mName = 'QQeval'
    
    return(evals.QQeval(org_data, imp_data, g, dName, mName=mName, pl_g=pl_g))

def test():
    org, imp, genes = evals.QQtest()
    return(org, imp, genes)
