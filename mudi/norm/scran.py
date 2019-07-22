import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from ..repos.anndata2ri import anndata2ri

pandas2ri.activate()
anndata2ri.activate()

import scanpy as sc
import scanpy.external as sce

scran = robjects.r('''
        scran <- function(M, groups, min.mean=0.1) {
            require(scran)
            return(computeSumFactors(M, clusters=groups, min.mean=min.mean))
        }
        ''')

def pyscran(adata_i, resolution=0.5, hvg=None):
    """
    Wraps scran.
        See here for paper: https://doi.org/10.1186/s13059-016-0947-7
    ------------------------------------
    Requires raw count data to be stored in layers.

    Uses scran to compute size factors, store these size factors as an observation for each cell.
    Scales the data by the size factors.
    Log norms data and computes highly variable genes.

    Returns the anndata object.
    """
    # assert adata_i.layers['counts'], "Please provide raw counts in AnnData Object."

    if hvg is None:
        hvg = {'min_mean':0.0125, 'max_mean':3, 'min_disp':0.5}

    adata = adata_i.copy()
    sc.tl.louvain(adata, key_added='scran_groups', resolution=resolution)

    adata.obs['size_factors'] = scran(
        adata.layers['counts'].toarray().astype(int).T,
        adata.obs['scran_groups']
    )

    adata.X = np.array(adata.layers['counts'] / adata.obs['size_factors'].values[:,np.newaxis])
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, **hvg)

    return adata
