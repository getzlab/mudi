import numpy as np
import pandas as pd
import os
import sys

"""
Storey Q-Values - https://github.com/StoreyLab/qvalue
--------------------
Python Wrapper
Author: Francois Aguet
https://github.com/broadinstitute/tensorqtl/blob/master/tensorqtl/rfunc.py
"""
def qvalue(p, lambda_qvalue=None):
    """Wrapper for qvalue::qvalue"""
    import rpy2
    from rpy2.robjects.packages import importr
    from collections import Iterable

    qvalue = importr("qvalue")
    rp = rpy2.robjects.vectors.FloatVector(p)
    if lambda_qvalue is None:
        q = qvalue.qvalue(rp)
    else:
        if not isinstance(lambda_qvalue, Iterable):
            lambda_qvalue = [lambda_qvalue]
        rlambda = rpy2.robjects.vectors.FloatVector(lambda_qvalue)
        q = qvalue.qvalue(rp, **{'lambda':rlambda})
    qval = np.array(q.rx2('qvalues'))
    pi0 = np.array(q.rx2('pi0'))[0]
    return qval, pi0

def t_test(mat: pd.DataFrame, group_s: pd.Series, equal_var: bool = False) -> pd.DataFrame:
    """
    t-test
    ---------------------
    Args:
        * mat: pd.DataFrame (genes x samples)
        * group_s: series of groupings
        * equal_var: wald-ttest (False)
    """
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    mat = mat[group_s.index]

    def _collapser(x, index, columns, name):
        _df = pd.DataFrame(x, index=index, columns=columns).reset_index()
        _id = _df.columns[0]
        return pd.melt(
            pd.DataFrame(x, index=index, columns=columns).reset_index(),
            id_vars=_id,
            ).set_index(_id).rename(columns={'variable':group_s.name,'value':name})

    groups = np.array(group_s)
    X = mat.values

    n_groups = np.unique(groups).shape[0]
    n_genes = X.shape[0]

    # Init np.arrays
    t_stat = np.zeros((n_genes, n_groups))
    pval = np.zeros((n_genes, n_groups))
    pval_adj = np.zeros((n_genes, n_groups))
    qval = np.zeros((n_genes, n_groups))
    x_in = np.zeros((n_genes, n_groups))
    x_out = np.zeros((n_genes, n_groups))

    for idx,group in enumerate(np.unique(groups)):
        mask = groups==group
        if sum(mask) > 1:
            X_in = X[:,mask]
            X_out = X[:,~mask]

            t_stat[:,idx], pval[:,idx] = stats.ttest_ind(X_in, X_out, axis=1, equal_var=equal_var)

            _,pval_adj[:,idx],_,_ = multipletests(
                pval[:,idx],
                alpha=0.05,
                method='fdr_bh',
                is_sorted=False,
                returnsorted=False
            )

            qval[:,idx],_ = qvalue(pval[:,idx])
            x_in[:,idx] = np.mean(X_in,1)
            x_out[:,idx] = np.mean(X_out,1)

    # Collapse to dataframe
    de_df = pd.concat([
                _collapser(x_in, mat.index, np.unique(groups), 'x_in'),
                _collapser(x_out, mat.index, np.unique(groups), 'x_out')['x_out'],
                _collapser(t_stat, mat.index, np.unique(groups), 't')['t'],
                _collapser(pval, mat.index, np.unique(groups), 'pval')['pval'],
                _collapser(pval_adj, mat.index, np.unique(groups), 'pval_adj')['pval_adj'],
                _collapser(qval, mat.index, np.unique(groups), 'qval')['qval']
            ],1)

    # Fold-change
    de_df['diff'] = de_df['x_in'] - de_df['x_out']

    # Signed FC * -log10(qval)
    de_df['gsea_rank'] = de_df['diff'] * -np.log10(de_df['pval_adj'])

    return de_df

def mannwhitneyu(mat: pd.DataFrame, group_s: pd.Series) -> pd.DataFrame:
    """
    mannwhitneyu
    ---------------------
    Args:
        * mat: pd.DataFrame (genes x samples)
        * group_s: series of groupings
    """
    from tqdm import tqdm
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    from sys import stdout

    mat = mat[group_s.index]

    def _collapser(x, index, columns, name):
        _df = pd.DataFrame(x, index=index, columns=columns).reset_index()
        _id = _df.columns[0]
        return pd.melt(
            pd.DataFrame(x, index=index, columns=columns).reset_index(),
            id_vars=_id,
            ).set_index(_id).rename(columns={'variable':group_s.name,'value':name})

    groups = np.array(group_s)
    X = mat.values

    n_groups = np.unique(groups).shape[0]
    n_genes = X.shape[0]

    # Init np.arrays
    u_stat = np.zeros((n_genes, n_groups))
    pval = np.zeros((n_genes, n_groups))
    pval_adj = np.zeros((n_genes, n_groups))
    qval = np.zeros((n_genes, n_groups))
    x_in = np.zeros((n_genes, n_groups))
    x_out = np.zeros((n_genes, n_groups))

    for idx,group in enumerate(np.unique(groups)):
        stdout.write("\r{} of {}".format(idx+1, n_groups))
        mask = groups==group
        if sum(mask) > 1:
            X_in = X[:,mask]
            X_out = X[:,~mask]

            for gn in range(X_in.shape[0]):
                #u_stat[gn,idx], pval[gn,idx] = stats.mannwhitneyu(X_in[gn], X_out[gn])
                u_stat[gn,idx], pval[gn,idx] = stats.mannwhitneyu(X_in[gn], X_out[gn], alternative='two-sided')

            _,pval_adj[:,idx],_,_ = multipletests(
                pval[:,idx],
                alpha=0.05,
                method='fdr_bh',
                is_sorted=False,
                returnsorted=False
            )

            try:
                qval[:,idx],_ = qvalue(fgsea_df['pval'].values)
            except:
                try:
                    qval[:,idx],_ = qvalue(fgsea_df['pval'].values, lambda_qvalue=0.5)
                except:
                    qval[:,idx] = None

            x_in[:,idx] = np.mean(X_in,1)
            x_out[:,idx] = np.mean(X_out,1)

    # Collapse to dataframe
    de_df = pd.concat([
                _collapser(x_in, mat.index, np.unique(groups), 'x_in'),
                _collapser(x_out, mat.index, np.unique(groups), 'x_out')['x_out'],
                _collapser(u_stat, mat.index, np.unique(groups), 'u')['u'],
                _collapser(pval, mat.index, np.unique(groups), 'pval')['pval'],
                _collapser(pval_adj, mat.index, np.unique(groups), 'pval_adj')['pval_adj'],
                _collapser(qval, mat.index, np.unique(groups), 'qval')['qval']
            ],1)

    # Fold-change
    de_df['diff'] = de_df['x_in'] - de_df['x_out']

    # Signed FC * -log10(qval)
    de_df['gsea_rank'] = de_df['diff'] * -np.log10(de_df['pval_adj'])

    return de_df
