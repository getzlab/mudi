import numpy as np
import scipy
from scipy import sparse
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import scanpy as sc
import scanpy.external as sce
import sys
import scrublet as scr
from typing import Union
from anndata import AnnData
from inspect import signature

from .utils import scanpy_adata_loader
from .utils import score_doublets
from .utils import score_cc_genes

# For Numba
import warnings
warnings.filterwarnings('ignore')

def bcs_by_group(
    obs: pd.DataFrame,
    group: str = 'percent_mito',
    key: str = 'louvain',
    thresh: float = 2,
    verbose:bool = False
    ) -> list:
    """
    Barcodes by Group
    ----------------------------
    Filters out quality thresholds for value of "group" variable.

    Inputs:
        - obs: observations dataframe from AnnData object
        - group: group to filter out barcodes by
        - key: what sub-groups to separately consider when doing filtering
        - thresh: std. err threshold cutoff for removing barcodes
        - verbose: print out filtering

    Outputs:
        - bcs: selected barcodes
    """
    bcs = list()

    for opt in set(obs[key]):
        obs_sub = obs[obs[key]==opt]
        d = obs_sub[group].values
        filt_thresh = np.mean(d) + thresh * np.std(d)
        obs_sub = obs_sub[obs_sub[group] < filt_thresh]

        if verbose:
            print("\tFiltering {} group {} - {} < {}".format(key, opt, group, filt_thresh))

        bcs += list(obs_sub.index)

    if verbose: print("Filtered {} / {} barcodes by {}.".format(obs.shape[0]-len(bcs), obs.shape[0], group))

    return bcs

def filter_upper(
    adata: AnnData,
    groups: list,
    verbose: bool,
    **kwargs
    ) -> list:
    """
    Filter Upper
    ----------------------------
    Performs downstream clustering on unfiltered data to estimate thresholds
    for filtering out diseased cells.

    Inputs:
        - adata: AnnData object
        - groups:  groups to do separate filtering on
        - **kwargs: key-word args for bcs_by_group

    Outpts:
        - final list of selected barcodes
    """
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=10000)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    adata = adata[:,adata.var['highly_variable']]
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, knn=True, n_pcs=40)
    sc.tl.louvain(adata, resolution=1)

    return list(
        set.intersection(*[set(bcs_by_group(adata.obs, group=g, verbose=verbose, **kwargs)) for g in groups])
        )

def recipe(
    file_name: Union[AnnData, str, list], \
    min_genes: int = 200, \
    min_cells: int = 3, \
    thresh: float = 1.25, \
    mito_thresh: Union[None,float] = None, \
    groups: Union[None, str] = None, \
    genome: Union[None, str] = None, \
    norm: str = "library", \
    scran_key: Union[None, str] = None, \
    scran_res: float = 0.5, \
    regress_vars: Union[None, str] = None, \
    regress_jobs: int = 1, \
    compute_doublets: bool = True, \
    remove_doublets: bool = False, \
    scrublet_key: str = 'batch', \
    hvg: Union[None, dict] = None, \
    qc: bool = False,
    downstream: bool = True,
    bbknn: Union[None, str] = None,
    issparse: bool = False,
    verbose: bool = False,
    make_sparse: bool = True,
    **kwargs
    ) -> AnnData:
    """
    Recipe for single-cell processing.
    ----------------------------
    Inputs:
        - min_genes: minimum number of genes required for a cell
        - min_cells: minimum number of cells for expression of a gene
        - thresh: estimated threshold for qc-filtering
        - groups: qc-metrics to threshold over (default is percent_mito)
        - genome: genome build for loading scanpy object (usually not needed)
        - norm: normalization method
            * library: default
            * scran: scran
        - scran_key: scran key to normalize on - if not provided, will do course clustering
        - scran_res: scran resolution for course clustering
        - regress_vars: variables to regress over
            *** Note: this will subset adata.X to highly variable genes
            *** Note: raw data will still be stores in adata.raw and adata.layers['counts']
        - regress_jobs: n_jobs to use for regression if specified
        - compute_doublets: run scrublet
        - remove_doublets: remove doublets before downstream processing
        - hvg: dictionary of keys specifying highly variable gene selection
        - qc: if True, returns adata object pre-filtering before downstream processing
            *** Note: good for easily computing qc-metrics
        - downstream: if True, continues with downstream processing
        - bbknn: if specified, performs light-weight batch correction on the provided
            variable
        - issparse: if provided an AnnData object, if hte X variable is alreayd in sparse format
        - make_sparse: if True, converts the hte X variable of the AnnData to be returned to sparse format

    Outputs:
        - adata: AnnData Object
    """
    # ---------------------------------
    # Argument Helpers
    # ---------------------------------
    if qc:
        min_cells = 0
        min_genes = 0

    if hvg is None:
        hvg = {'min_mean':0.0125, 'max_mean':3, 'min_disp':0.5}

    for key in signature(sc.pp.highly_variable_genes).parameters:
        if key in kwargs:
            raise KeyError("Please place {} in hvg.".format(key))

    if groups is None:
        groups = ['percent_mito']

    if remove_doublets:
        assert compute_doublets, 'Doublet removal specified but doublets are not being computed.'

    # ---------------------------------
    # Pipeline
    # ---------------------------------
    if isinstance(file_name, AnnData):
        adata = file_name
        if not issparse:
            adata.X = sparse.csr_matrix(adata.X)
    else:
        adata = scanpy_adata_loader(file_name, genome=genome)

    adata.var_names_make_unique()

    # General QC
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    if 'batch' not in list(adata.obs):
        adata.obs['batch'] = '1'

    # ---------------------------------
    # Compute Doublets
    # ---------------------------------
    if compute_doublets:
        score_doublets(adata, key=scrublet_key)

        if remove_doublets:
            print("Dropping {} doublets".format(sum(adata.obs['doublet'])))
            adata = adata[~adata.obs['doublet']]

    # ---------------------------------
    # Compute Mitochondrial + Ribo Genes
    # ---------------------------------
    mito_genes = adata.var_names.str.startswith('MT-')
    ribo_genes = adata.var_names.str.startswith('RPS') + adata.var_names.str.startswith('RPL')

    # Filter by grouped mito-thresholds
    adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    adata.obs['percent_ribo'] = np.sum(adata[:, ribo_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    m = adata.obs['percent_mito']

    if qc:
        return adata

    # ---------------------------------
    # Filter out barcodes
    # ---------------------------------
    if thresh is not None:
        adata = adata[adata.obs.index.isin(filter_upper(adata.copy(), groups, thresh=thresh, verbose=verbose, **kwargs))]

    # Mito Threshold
    if mito_thresh is not None:
        if mito_thresh == 'auto':
            mito_thresh = np.std(m)*1.5 + np.mean(m)
        print("\tFiltering mitochondiral genes over {}.".format(mito_thresh))
        prev_cells = adata.shape[0]
        adata = adata[adata.obs['percent_mito'] < mito_thresh]
        print("Filtered {} / {} barcodes.".format(prev_cells-adata.shape[0], prev_cells))

    adata.layers['counts'] = adata.X.copy()

    # ---------------------------------
    # Normalization
    # ---------------------------------
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=10000)
    sc.pp.log1p(adata)

    # Save Raw
    adata.raw = adata.copy()

    if norm=='library':
        pass
    elif norm=='scran':
        from .norm import pyscran

        adata = pyscran(adata, resolution=scran_res, scran_key=scran_key, hvg=hvg, log_norm=True)

    # ---------------------------------
    # Score cell-cycle genes
    # ---------------------------------
    try:
        score_cc_genes(adata)
    except:
        if verbose: print("Unable to compute cell-cycle genes.")
        pass

    # ---------------------------------
    # Downstream processing
    # ---------------------------------
    if downstream:
        sc.pp.highly_variable_genes(adata, **hvg)

        # Regress out vars if provided
        if regress_vars is not None:
            adata = adata[:, adata.var['highly_variable']]
            sc.pp.regress_out(adata, regress_vars, n_jobs=regress_jobs)

        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True)

        if bbknn is not None:
            # Light-weight batch correction
            sce.pp.bbknn(adata, batch_key=bbknn)
        else:
            sc.pp.neighbors(adata)

        sc.tl.louvain(adata, resolution=1)
        sc.tl.umap(adata)

    # Convert back to sparse if need-be
    if make_sparse:
        try:
            adata.X = sparse.csr_matrix(adata.X)
        except:
            pass

    return adata
