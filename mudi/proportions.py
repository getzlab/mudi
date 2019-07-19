import numpy as np
import sys
import os
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from anndata import AnnData
import os

from scipy.sparse import vstack
from decimal import Decimal
from tqdm import tqdm_notebook, tqdm

def plot_proportions(df, ax):
    """
    Plot proportions on a given axis.
    """
    for cell_type in list(df):
        sns.distplot(df[cell_type], label=cell_type, ax=ax)
    ax.set_xlabel('')
    ax.set_ylabel('Freq')

def sample_adata_prop(adata, key='batch', label='label', n_bootstrap='min', **kwargs):
    """
    Boot-strap sampling of cell propotions.

    Takes in an anndata object and a batch key for if different samplings are
    desired (i.e. pre-treatment vs. post-treatment).

    Returns sample_proportions.

    """
    obs_data = adata.obs.reset_index().loc[:,[key,label]].copy()

    batches = list(set(obs_data[key]))
    batch_indices = [obs_data[obs_data[key]==batch].index.values for batch in batches]

    if n_bootstrap == 'min':
        n_bootstrap = min([b.shape[0] for b in batch_indices])

    sample_proportions = {}

    for idx,batch in enumerate(batches):
        sample_proportions[batch] = categorical_sampler(obs_data[label].values, \
                                                        batch_indices[idx],
                                                        n_bootstrap, **kwargs)

    return sample_proportions

def categorical_sampler(x, indices, n_bootstrap, n_iter=1000, replace=True):
    """
    Samples by catgegory to generate distribution of proportions.
    """
    M = np.zeros((n_iter, len(set(x))))
    _map = {l:idx for idx,l in enumerate(np.unique(x))}

    for i in range(n_iter):
        names, counts = np.unique(x[np.random.choice(indices, n_bootstrap, replace=replace)], return_counts=True)
        M[i,[[_map[name] for name in names]]] += counts

    return pd.DataFrame(data = M / n_bootstrap, columns=np.unique(x))

def calc_empirical_pval(x,y,two_way=True):
    """
    Compute Empirical p-values.

    Takes two n x c matrices.
        n: number of bootstrapped random samples
        c: number of categories

    Returns a dataframe with the categories and computed p-values from the distribution.
    If conducting a one-sided test, values of signifance are for up-regulated
    cell-popluations from before treatment to after treatment.
    """
    diffs = pd.DataFrame(np.sum(x-y > 0) / x.shape[0])
    if two_way:
        diffs[1] = 1-diffs[0]
        return pd.DataFrame(data = np.array(2*(1 - diffs[[0,1]].max(1))), index=diffs.index, columns=['pval'])
    else:
        return pd.DataFrame(data = diffs.values, index=diffs.index, columns=['pval'])
