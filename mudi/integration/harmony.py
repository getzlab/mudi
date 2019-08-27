import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from ..repos.anndata2ri import anndata2ri

pandas2ri.activate()
anndata2ri.activate()

import scanpy as sc
import scanpy.external as sce

harmony = robjects.r('''
        harmony <- function(pca,batch,theta=2,lam=1) {
            library(harmony)
            library(magrittr)

            hem <- HarmonyMatrix(pca, batch, theta=theta, lambda=lam, npcs = 20)
            hem = data.frame(hem)
            return(hem)
        }
        ''')

def pyharmony(adata_i, batch_key='batch', theta=2, lam=1):
    """
    Python wrapper for Harmony method for scRNA-seq data integration.

    Inputs:
        - adata_i: input Anndata object; pre-compute PCA of samples or it will automatically
                    compute this
        - batch_key: key to integrate data on
        - theta: harmony parameter (theta)
        - lam: harmony parameter (lambda)

    Outputs:
        - adata: returns an Anndata object with corrected PCA and re-computed UMAP based on this
                  PCA representation.

    """
    adata = adata_i.copy()
    try:
        pca = adata.obsm['X_pca']
    except:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True)

    batch = adata.obs[batch_key]
    hem = harmony(pca, batch, theta=theta, lam=lam)
    adata.obsm['X_pca'] = hem.values
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    return adata
