from pyscenic.aucell import aucell
from pyscenic.genesig import GeneSignature
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def create_gene_signatures(W_df: pd.DataFrame, subtype_idx: str = 'Subgroup', weight_idx: str = 'max_norm'):
    """
    Create gene signatures.
    ---------------------------
    Args:
        * W_df: weights dataframe
            --> (genes x weights,subtype)
        * subtype_idx: str column in W_df subtype to make individual gene signatures for
        * weight_idx: str column in W_df weight to use
    Returns:
        * Dictionary
            * subtype --> pyscenic.genesig.GeneSignature
    """
    return [GeneSignature(
                sig,
                W_df[W_df[subtype_idx]==sig][[weight_idx]].to_dict()[weight_idx]
            ) for sig in np.unique(W_df[subtype_idx])]

def auc(
    adata,
    markers_df: pd.DataFrame,
    tag: str = 'enrichment',
    num_workers: int = 1,
    normalize: bool = True,
    inplace: bool = True,
    subtype_idx: str = 'Subgroup',
    weight_idx: str = 'max_norm',
    **kwargs):
    """
    Run AUCell on AnnData Object.
    --------------------------
    Args:
        * adata: Scanpy AnnData Object
        * markers_df: pd.DataFrame of markers (genes x weights, subtype)
        * tag: name of enrichment
        * num_workers: how many threads to use for computing enrichments
        * normalize: normalize enrichment score
        * inplace: add results directly to adata.obs
            *  False: return  pd.DataFrame of enrichments
        * subtype_idx: str column in markers_df subtype to make individual gene signatures for
        * weight_idx: str column in markers_df weight to use
        **kwargs: passed to aucell

    Returns:
        * None or pd.DataFrame of enrichments
    """
    # Drop NA markers
    markers_df = markers_df.loc[markers_df.index.dropna()]

    #Create gene signatures
    gene_sigs = create_gene_signatures(markers_df, subtype_idx=subtype_idx, weight_idx=weight_idx)

    # Run AUCell
    exp_mtx = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
    e_df = aucell(exp_mtx, gene_sigs, normalize=normalize, num_workers=num_workers, **kwargs)
    e_df = e_df.rename(columns={x:tag+"_"+x for x in e_df.columns})
    e_df[tag+'_subtype'] = e_df.idxmax(1).apply(lambda x: x.split("_")[-1])

    # Attach to  adata object
    if inplace:
        adata.obs = adata.obs.join(e_df)
        return
    else:
        return e_df

def assign_bootstrap(filepath, n=100, norm=False):
    """
    Assign Bootstrap Labels
    ---------------------
    Args:
        * filepath: path to h5 file from bootstrap runs
        * n: number of boostraps
        * norm: whether to normalize results
    """
    def empirical_pval(a,b):
        return 1-np.sum(a-b > 0) / a.shape[0]

    perm = np.array([pd.read_hdf(filepath,'perm{}'.format(x)).values for x in range(n)])

    if norm:
        perm = perm / np.sum(perm, axis=2)[:,:,np.newaxis]

    groups = pd.read_hdf(filepath,'perm0').columns
    barcodes = pd.read_hdf(filepath,'perm0').index

    perm_mean = np.mean(perm,0)
    perm_std = np.std(perm,0)

    max_one = np.argsort(perm_mean,axis=1)[:,-1]
    max_two = np.argsort(perm_mean,axis=1)[:,-2]

    pvals = np.array([empirical_pval(perm[:,i,max_one[i]],perm[:,i,max_two[i]]) for i in range(perm.shape[1])])
    assignments = np.array(groups[max_one])

    X = np.concatenate((assignments[:,np.newaxis],pvals[:,np.newaxis],perm_mean),axis=1)
    cols = np.hstack((np.array(['aucell','pval']),groups))

    return pd.DataFrame(X, index=barcodes, columns=cols)
