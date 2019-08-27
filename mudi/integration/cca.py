import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from ..repos.anndata2ri import anndata2ri

pandas2ri.activate()
anndata2ri.activate()

import scanpy as sc
import scanpy.external as sce

cca = robjects.r('''
    cca <- function(df, batch_ids) {
        suppressMessages(library(Seurat))

        print("Building Seurat Object...")
        # Create Seurat Object
        df = data.matrix(df)
        sdf = CreateSeuratObject(df)
        sdf@meta.data$batch = batch_ids

        # Split object into batches
        sdf.list <- SplitObject(sdf, split.by = "batch")

        # Find variable features
        for (i in 1:length(sdf.list)) {
            sdf.list[[i]] <- FindVariableFeatures(sdf.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
        }

        # Find integration anchors
        batch.anchors <- FindIntegrationAnchors(object.list = sdf.list, dims=1:30)

        # Integrate data
        batch.combined <- IntegrateData(anchorset = batch.anchors, dims = 1:30)

        # Return integrated expression matrix
        return(data.frame(batch.combined$integrated@data))
    }
    ''')

def pycca(adata_i, batch_key='batch'):
    """
    Python wrapper for Seurat's CCA.
    ------------------------------------
    Inputs:
        - adata: input adata object
        - batch_key: batch key for integration over

    Outputs:
        - Returns an adata object with:
            - Original counts stored as the .raw attribute
            - adata.X --> corrected matrix values
    """
    adata = adata_i.copy()

    df = pd.DataFrame(data=adata.X.todense().transpose(), \
                      index=adata.var_names, \
                      columns=adata.obs_names)

    batch_ids = np.array(adata.obs[batch_key].values, dtype=str)
    cca_X = cca(df, batch_ids)

    adata = adata[:,cca_X.index]
    adata.X = cca_X.T.values

    adata.raw = adata_i
    sc.pp.scale(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.louvain(adata)

    return adata
