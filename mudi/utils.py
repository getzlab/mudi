import numpy as np
import h5py
import scipy
import gc
import pandas as pd
import os
import time
import pkg_resources

import scanpy as sc
import scanpy.external as sce

import sys
import scrublet as scr

# ---------------------------------
# Scanpy Helpers
# ---------------------------------
def scanpy_adata_loader(path, genome='GRCh38', verbose=True):
    """
    Loader function.
    ------------------
    Can handle lists of file/dir paths, a file (.h5) or a directory format.
    Use this to load/aggregate immediately into a scanpy object.
    """
    if isinstance(path, list):
        if verbose:
            print("Combining {} inputs.".format(len(path)))
        tmp = [scanpy_adata_loader(f, genome=genome) for f in path]
        return tmp[0].concatenate(tmp[1:])

    if os.path.isfile(path):
        return sc.read_10x_h5(path, genome=genome)
    elif os.path.isdir(path):
        return sc.read_10x_mtx(path)
    else:
        raise FileError("Provide proper path.")

def aggr_markers(adata, uns='rank_genes_groups', params=['names','scores','pvals','pvals_adj','logfoldchanges']):
    """
    Aggregate markers.
    ------------------
    Returns an easy to view marker list dataframe.
    Assumes 'rank_genes_groups' has already been called to find group markers
    in AnnData Object.
    """
    assert adata.uns['rank_genes_groups'], 'Compute differentially expressed genes first.'

    joint = pd.concat(
        [pd.DataFrame(pd.DataFrame(adata.uns[uns][param]).T.stack()).reset_index().rename(columns={'level_0':'cluster','level_1':'no.',0:param})
            for param in params],
        1)

    return joint.loc[:,~joint.columns.duplicated()]

# ---------------------------------
# Utilities
# ---------------------------------
def score_cc_genes(adata, cc_genes_file=pkg_resources.resource_filename('mudi', './ref/cell_cycle_genes/Macosko_cell_cycle_genes.txt')):
    """
    Score Cell-Cycle Genes
    ------------------------------------
    How to run:
        score_cc_genes(adata)

    Loads cell cycle genes list (ex. Macosko et al 2015) and runs cycle-scoring
    on input anndata. Does everything in place. Stores the following in .obs:
        - S_score
        - G2M_score
        - phase
    """
    cc_genes = pd.read_table(cc_genes_file, delimiter='\t')

    s_genes = cc_genes['S'].dropna()
    g2m_genes = cc_genes['G2.M'].dropna()

    s_genes_i = adata.var_names[np.in1d(adata.var_names, s_genes)]
    g2m_genes_i = adata.var_names[np.in1d(adata.var_names, g2m_genes)]

    sc.tl.score_genes_cell_cycle(adata, s_genes_i, g2m_genes_i)

def score_doublets(adata, key='batch', n_prin_comps=20, verbose=False):
    """
    Scrubber: wrapper for Scrublet.
    ------------------------------------
    How to run:
        score_doublets(adata)

    Adds the following to anndata object:
        - adata.obs['scrublet_score'] --> float (0 - 1.0)
        - adata.obs['doublet'] --> bool
    """
    doublet_scores, predicted_doublets = list(),list()

    for batch in adata.obs[key].drop_duplicates().values:
        scrub = scr.Scrublet(adata[adata.obs[key]==batch].X)
        _doublet_scores, _predicted_doublets = scrub.scrub_doublets(n_prin_comps=n_prin_comps, verbose=verbose)
        doublet_scores.append(_doublet_scores)
        predicted_doublets.append(_predicted_doublets)

    adata.obs['scrublet_score'] = np.concatenate(doublet_scores)
    adata.obs['doublet'] = np.concatenate(predicted_doublets)
