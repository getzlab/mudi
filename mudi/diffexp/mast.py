import rpy2.robjects as robjects
from anndata import AnnData
import anndata2ri

import rpy2.rinterface_lib.callbacks
from rpy2.robjects import pandas2ri
pandas2ri.activate()
anndata2ri.activate()

mast = robjects.r('''
        mast <- function(adata, latent.vars=c('n_counts'),  stim='condition', subset=NULL, rename='Group2') {
            require(MAST)
            # Convert scanpy adata to a single-cell assay required by MAST
            sca <- SceToSingleCellAssay(adata, class = "SingleCellAssay")
            colData(sca)$n_genes = scale(colData(sca)$n_genes)

            # Subset data if provided
            if (!is.null(subset)) {
                 #print('subsetting')
                 sca <- subset(sca, with(colData(sca), label==subset))
            }

            # Latent variables
            latent.vars <- c(stim, latent.vars)

            # Formula for zlm model
            fmla <- as.formula(
                object = paste0(" ~ ", paste(latent.vars, collapse = "+"))
            )

            # Zlm model
            zlmCond <- zlm(formula = fmla, sca = sca)
            summaryCond <- summary(zlmCond, doLRT=paste(c(stim,rename),collapse=''))
            summaryDt <- summaryCond$datatable

            result <- merge(summaryDt[contrast==paste(c(stim,rename),collapse='') & component=='H',.(primerid, `Pr(>Chisq)`)], #P-vals
                 summaryDt[contrast==paste(c(stim,rename),collapse='') & component=='logFC', .(primerid, coef)],
                 by='primerid') #logFC coefficients

            result[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
            de = result[result$FDR<1,, drop=F]
            de = de[order(de$FDR),]

            return(de)
        }
        ''')

def set_var(x,y):
    if x == y:
        return "Group2"
    else:
        return "Group1"

def pymast(adata_i, var, var_set, new='condition', latent_vars=['n_counts'], subset=None):
    """
    Wrapper for MAST single-cell differential expression analysis.
    Wraps R kernel.

    Inputs:
        adata_i: provide an anndata object
        var: variable of interest (i.e. cell-type)
        subset: subset key (i.e. monocytes)

    Outputs:
        Padnas Dataframe of de genes.

    """
    adata = adata_i.copy()
    adata.obs[new] = adata.obs[var].apply(lambda x: set_var(x,var_set))

    if subset is None:
        de = mast(adata, latent_vars, stim=new)
    else:
        de = mast(adata, latent_vars, stim=new, subset=subset)

    de['coef'] = de['coef']*-1
    return de
