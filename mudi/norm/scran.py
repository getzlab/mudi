import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from ..repos.anndata2ri import anndata2ri

pandas2ri.activate()
anndata2ri.activate()

import scanpy as sc
import scanpy.external as sce
import numpy as np
from typing import Union

scran = robjects.r('''
        scran <- function(M, groups, min.mean=0.1) {
            require(scran)
            return(computeSumFactors(M, clusters=groups, min.mean=min.mean))
        }
        ''')

def pyscran(
    adata_i: AnnData,
    resolution: float=0.5,
    hvg: Union[dict, None] = None,
    log_norm: bool = True,
    scran_key: Union[str, None] = None):
    """
    Wraps scran.
        See here for paper: https://doi.org/10.1186/s13059-016-0947-7
    ------------------------------------
    Requires raw count data to be stored in layers.

    Uses scran to compute size factors, store these size factors as an observation for each cell.
    Scales the data by the size factors.

    Args:
        * adata_i: input AnnData object
        * resolution: float describing the lovain clustering to do on the raw counts data for
            selecting of scran groups
        * hvg: dictionary with parameters for selecting highly variable genes; defaults
            to Seurat settings
        * log_norm: whether or not to automatically log-norm data
        * scran_key: found in .obs; this will be used for scran groups if provided

    Returns:
        * adata: AnnData object with scran-normalized counts
            Automatically returns log-normed w/ HVGs; if 'log_norm = False', will return
            only the scran-normalized counts
    """
    # assert adata_i.layers['counts'], "Please provide raw counts in AnnData Object."

    if hvg is None:
        hvg = {'min_mean':0.0125, 'max_mean':3, 'min_disp':0.5}

    adata = adata_i.copy()

    if scran_key is None:
        sc.pp.neighbors(adata)
        sc.tl.louvain(adata, key_added='scran_groups', resolution=resolution)
    else:
        adata.obs['scran_groups'] = adata.obs[scran_key]

    print(adata.obs.groupby('scran_groups').size())



    adata.obs['size_factors'] = scran(
        adata.layers['counts'].toarray().astype(int).T,
        adata.obs['scran_groups']
    )

    adata.X = np.array(adata.layers['counts'] / adata.obs['size_factors'].values[:,np.newaxis])

    if log_norm:
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, **hvg)

    return adata
