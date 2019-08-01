import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from ..repos.anndata2ri import anndata2ri

import scanpy as sc
import scanpy.external as sce
import numpy as np

pandas2ri.activate()
anndata2ri.activate()

sct = robjects.r('''
        sct <- function(adata) {
            require(Seurat)
            require(sctransform)

            dat <- as.Seurat(adata, counts='X', data='X')

            # run sctransform & regress out percent mitochondrial reads
            dat <- SCTransform(dat, verbose=TRUE)

            # Corrected log-transformed UMI
            genes <- dat$SCT@data@Dimnames[[1]]
            bcs <- dat$SCT@data@Dimnames[[2]]
            nmat <- dat$SCT@data

            genes_dataframe <- dat$SCT[[]]
            X <- dat$SCT@scale.data

            return(list(nmat,genes,bcs,genes_dataframe,X))
        }
        ''')

def pysctform(adata_i, num_hvg=3000):
    """
    Wraps sctransform.
        See here for preprint: https://www.biorxiv.org/content/10.1101/576827v2
        See here for tutorial: https://rawgit.com/ChristophH/sctransform/master/inst/doc/seurat.html
    ------------------------------------
    Requires raw count data to be stored in layers.

    Uses sctransform to compute new normalized count matrix.
        - Replaces normalize data
        - Replaces scaling data
        - Replaces finding highly variable features

    The rpy2 function sct returns a list:
        - results[0]: corrected, log-normalized UMI matrix (data)
        - results[1]: genes
        - results[2]: barcodes
        - results[3]: genes dataframe with information about residual variance
        - results[4]: residuals (normalized values) - non-sparse (scale.data)

    Following this transform, all genes stored in the residuals are highly-variable.
    Information about the variance contributed by each gene is stored in:
        - adata.uns['sct']
        - array with the following column names:
            detection_rate	gmean	variance	residual_mean	residual_variance

    Corrected, log-normalized UMI matrix:
        - adata.X

    Normalized residuals (to use for PCA):
        - adata.obsm['X_sct']

    ------------------------------------
    Note from the Seurat website:
        You can use the corrected log-normalized counts for differential
        expression and integration. However, in principle, it would be most
        optimal to perform these calculations directly on the residuals
        (stored in the scale.data slot) themselves.

    """
    adata = adata_i.copy()

    # Hacky fix
    # ScTransform does not like variables with '__'
    # TODO: more robust solution for this
    adata = adata[:,np.array([x for x in adata.var_names if '_' not in x])]

    # To np.array for easy conversion to rpy2
    adata.X = np.asarray(adata.layers['counts'].todense(),dtype=int)
    del adata.layers['counts']
    adata.raw = None

    print("Running scTransform...")
    results = sct(adata)
    print("Finishing running scTransform...")

    adata_i = adata_i[:,results[1]]
    adata_i.X = results[0].T
    adata_i.uns['sct'] = results[3].sort_values('residual_variance', ascending=False)
    adata_i.obsm['X_sct'] = results[4].T

    # Add highly varible genes selected from ScTransform
    adata_i.var['highly_variable'] = adata_i.var_names.isin( \
        [x for x in adata_i.uns['sct'].index[:num_hvg] if \
         not x.startswith('MT-') and \
         not x.startswith('RPS') and \
         not x.startswith('RPL')] \
        )

    adata_i.raw = adata_i

    # Downstream
    print("Downstream processing...")
    sc.tl.pca(adata_i, use_highly_variable=True)
    sc.pp.neighbors(adata_i)
    sc.tl.louvain(adata_i, resolution=1.0)
    sc.tl.umap(adata_i)

    return adata_i
