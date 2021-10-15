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
        ad = sc.read_10x_h5(path, genome=genome)
        ad.var_names_make_unique()
        return ad
    elif os.path.isdir(path):
        ad = sc.read_10x_mtx(path)
        ad.var_names_make_unique()
        return ad
    else:
        raise FileError("Provide proper path.")

def get_percent_expr(adata, groupby):
    """
    Get percent expressed & mean expression.
    ------------------------------
    Requires:
        * adata.layers['counts'] -> counting percent of cells with gene expressed
        * adata.raw.X -> for computing mean expression (log1p)
    """
    from tqdm import tqdm

    groups = list(adata.obs[groupby].cat.categories)
    res_in = pd.DataFrame(columns=adata.var_names, index=groups)
    res_out = pd.DataFrame(columns=adata.var_names, index=groups)
    res_mean_in = pd.DataFrame(columns=adata.var_names, index=groups)
    res_mean_out = pd.DataFrame(columns=adata.var_names, index=groups)

    for group in tqdm(groups, desc="Computing metrics per group"):
        res_in.loc[group] = (adata[adata.obs[groupby].isin([group]),:].layers['counts'] > 0).mean(0)
        res_out.loc[group] = (adata[~adata.obs[groupby].isin([group]),:].layers['counts'] > 0).mean(0)
        res_mean_in.loc[group] = adata[adata.obs[groupby].isin([group]),:].raw.X.mean(0)
        res_mean_out.loc[group] = adata[~adata.obs[groupby].isin([group]),:].raw.X.mean(0)

    res_in = res_in.T
    res_out = res_out.T
    res_mean_in = res_mean_in.T
    res_mean_out = res_mean_out.T

    res_in = res_in.reset_index().melt(id_vars=['index']).rename(columns={'index':'names','variable':"group", 'value':'percent_in'})
    res_out = res_out.reset_index().melt(id_vars=['index']).rename(columns={'index':'names','variable':"group", 'value':'percent_out'})

    res_mean_in = res_mean_in.reset_index().melt(id_vars=['index']).rename(columns={'index':'names','variable':"group", 'value':'mean_expr_in'})
    res_mean_out = res_mean_out.reset_index().melt(id_vars=['index']).rename(columns={'index':'names','variable':"group", 'value':'mean_expr_out'})

    return pd.merge(pd.merge(res_in, res_out), pd.merge(res_mean_in, res_mean_out))

def aggr_markers(adata, uns='rank_genes_groups', expr_metrics=True):
    """
    Aggregate markers.
    ------------------
    Returns an easy to view marker list dataframe.
    Assumes 'rank_genes_groups' has already been called to find group markers
    in AnnData Object.
        * expr_metrics -> compute percent of cells expressed & mean expression for in/out groups.
    """
    assert adata.uns[uns], 'Compute differentially expressed genes first.'

    aggr_df = sc.get.rank_genes_groups_df(adata, None)

    if expr_metrics:
        aggr_percent_expr = get_percent_expr(adata, adata.uns[uns]['params']['groupby'])
        return pd.merge(aggr_df, aggr_percent_expr)
    else:
        return aggr_df

def get_de_genes_to_plot(markers_df, lfc_thresh=1, padj_thresh=0.1, n_to_plot=5):
    """
    Top DiffExp Genes.
    Return as dict for easy plotting with sc.pl.dotplot.
    """
    markers_df = markers_df[
        (markers_df['logfoldchanges']>=lfc_thresh) &
        (markers_df['pvals_adj']<=padj_thresh)
    ].groupby("group").head(n_to_plot)

    return markers_df.groupby("group").agg(list)['names'].to_dict()

def get_uns(adata, tag):
    """
    Retrieve unstructured data stored in AnnData.
    ------------------------
    Inputs:
        - adata: AnnData Object
        - tag: name of key in adata.uns
    Outputs:
        - pd.DataFrame: formatted information in adata.uns

    """
    assert tag in adata.uns, "{} not found in adata.uns".format(tag)
    try:
        return pd.DataFrame(adata.uns[tag]['values'], index=adata.uns[tag]['rows'], columns=adata.uns[tag]['cols'])
    except:
        raise ValueError("Unable to return structured dataframe from data.uns[{}]".format(tag))

def get_a_by_b(adata, a, b, norm=False):
    """
    Get A x B.
    ----------------
    Number of each .obs b per .obs a
        returns pd.Dataframe
    """
    hm = adata.obs.groupby([a,b]).size().reset_index().set_index(a).pivot(columns=b)

    if norm:
        hm = hm.div(hm.sum(1), 0)

    hm.columns = hm.columns.droplevel()
    hm.columns.name = None
    return hm
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
