from anndata import AnnData
import scanpy as sc
import numpy as np
import pandas as pd
import sys
import os

import matplotlib.pyplot as plt

from signatureanalyzer import ardnmf
from signatureanalyzer.utils import postprocess_msigs, get_nlogs_from_output, select_markers
from signatureanalyzer.consensus import consensus_cluster
from signatureanalyzer.plotting import k_dist, consensus_matrix
from signatureanalyzer.plotting import marker_heatmap

"""
ARD-NMF Wrapper.

See https://github.com/broadinstitute/getzlab-SignatureAnalyzer
For more details about the method and implementation.
"""

def join_nmf_output_to_anndata(adata: AnnData, filepath: str, cut_norm:float = 0, cut_diff:float = 0.1):
    """
    Join ARD-NMF output to AnnData
    -----------------------
    Args:
        * adata: AnnData object to join to
        * filepath: output file path from signatureanalyzer
        * cut_norm: marker selection param (signatureanalyzer.utils.select_markers)
        * cut_diff: marker selection param (signatureanalyzer.utils.select_markers)

    Returns:
        * None (edits adata directly)
    """
    # Load Results
    H = pd.read_hdf(filepath,"H")
    X = pd.read_hdf(filepath,"X")
    W = pd.read_hdf(filepath,"W")

    markers,signatures = select_markers(X, W, H, cut_norm=cut_norm, cut_diff=cut_diff)

    # Join to Obs
    adata.obs = adata.obs.join(H)
    adata.obs['max_id'] = adata.obs['max_id'].astype('category')

    try:
        adata.obs = adata.obs.join(pd.read_hdf(filepath,"consensus"))
        adata.obs['clusters'] = adata.obs['clusters'].astype('category')
        adata.obs = adata.obs.rename(columns={'clusters':'consensus_cluster'})
    except:
        pass

    adata.obs = adata.obs.rename(columns={'max_id':'nmf_id'})

    adata.uns['signatures'] = list(H)[:-3]
    adata.obsm['X_nmf'] = H.iloc[:,:-3].values

    # W matrix
    _nmf_genes = {}
    _nmf_genes['cols'] = list(signatures)
    _nmf_genes['rows'] = signatures.index
    _nmf_genes['values'] = signatures.values
    adata.uns['nmf_genes'] = _nmf_genes

    # Marker genes
    _nmf_markers = {}
    _nmf_markers['cols'] = list(markers)
    _nmf_markers['rows'] = markers.index
    _nmf_markers['values'] = markers.values
    adata.uns['nmf_markers'] = _nmf_markers

def NMF(
    adata: AnnData,
    nruns: int = 10,
    outdir: str = '.',
    input_type: str = 'raw',
    use_highly_variable: bool = True,
    filter_ribo: bool = True,
    filter_mito: bool = True,
    inplace: bool = True,
    verbose: bool = False,
    consensus_cluster_results: bool = False,
    **nmf_kwargs
    ):
    """
    AnnData wrapper for ARD-NMF.
    ------------------------
    Inputs:
        * adata: AnnData object
        * nruns: number of runs of NMF to complete
        * outdir: path (str) where to output ARD-NMF results
        * input_type: what AnnData input to use
            * 'raw': normalized in adata.raw
            * 'X': X
            * 'counts': raw counts saved in adata.layers['counts']
                * any input saved in adata.layers
        * use_highly_variable: whether or not to use highly variable genes
            (recommended)
        * filter_ribo: exclude ribosomal genes (default: True)
        * filter_mito: exclude mitochondrial genes (default: True)
        * inplace: whether or not to edit the AnnData object directly
            * False: returns H,W,markers,signatures (4 pd.DataFrame objects)
        * verbose: verbosity
        * consensus_cluster_results: run consensus clustering (default: False)
            * Very slow for 20K+ cells
        * nmf_kwargs: passed to signatureanalyzer.ardnmf

    Outputs:
        * Depends

    _ note _
    It is highly reccomended to use only highly variable genes for factorization.
    Default parameters are to use highly variable genes and filter out mitochondrial
    or ribosomal genes.
    """
    if outdir is not ".":
        print("   * Creating output dir at {}".format(outdir))
        os.makedirs(outdir, exist_ok=True)

    # ---------------------------------
    # Select Genes for NMF
    # ---------------------------------
    if use_highly_variable:
        assert 'highly_variable' in list(adata.var), 'Please compute highly variable genes.'
        print("   * Using {} highly variable genes.".format(sum(adata.var.highly_variable)))
        genes_for_nmf = set(adata.var[adata.var['highly_variable']].index)
    else:
        genes_for_nmf = set(adata.var.index)

    if filter_mito:
        mito_genes = {x for x in genes_for_nmf if x.startswith('MT-')}
        print("   * Filtering {} mito genes.".format(len(mito_genes)))
        genes_for_nmf -= mito_genes
    if filter_ribo:
        ribo_genes = {x for x in genes_for_nmf if x.startswith('RPS') or x.startswith('RPL')}
        print("   * Filtering {} ribo genes.".format(len(ribo_genes)))
        genes_for_nmf -= ribo_genes

    # Input type
    if input_type == 'raw':
        df = pd.DataFrame(data=adata.raw.X.toarray(), index=adata.obs_names, columns=adata.raw.var_names).T
    elif input_type == 'X':
        df = pd.DataFrame(data=adata.X.toarray(), index=adata.obs_names, columns=adata.var_names).T
    else:
        assert input_type in adata.layers, "Please save input in adata.layers['{}']".format(input_type)
        df = pd.DataFrame(data=adata.layers['counts'].toarray(), index=adata.obs_names, columns=adata.var_names).T

    matrix = df.loc[genes_for_nmf]

    # ---------------------------------
    # ARD-NMF
    #
    # Code chunk copied from
    # https://github.com/broadinstitute/getzlab-SignatureAnalyzer/blob/master/signatureanalyzer/signatureanalyzer.py
    # ------------------------------------------------------------------{
    print("   * Saving ARD-NMF outputs to {}".format(os.path.join(outdir,'nmf_output.h5')))
    store = pd.HDFStore(os.path.join(outdir,'nmf_output.h5'),'w')

    print("   * Running ARD-NMF...")
    for n_iter in range(nruns):
        store['X'] = matrix

        res = ardnmf(
            matrix,
            tag="\t{}/{}: ".format(n_iter,nruns-1),
            verbose=verbose,
            **nmf_kwargs
        )

        lam = pd.DataFrame(data=res["lam"], columns=["lam"])
        lam.index.name = "K0"

        store["run{}/H".format(n_iter)] = res["H"]
        store["run{}/W".format(n_iter)] = res["W"]
        store["run{}/lam".format(n_iter)] = lam
        store["run{}/Hraw".format(n_iter)] = res["Hraw"]
        store["run{}/Wraw".format(n_iter)] = res["Wraw"]
        store["run{}/markers".format(n_iter)] = res["markers"]
        store["run{}/signatures".format(n_iter)] = res["signatures"]
        store["run{}/log".format(n_iter)] = res["log"]

    store.close()

    # Select Best Result
    aggr = get_nlogs_from_output(os.path.join(outdir,'nmf_output.h5'))
    max_k = aggr.groupby("K").size().idxmax()
    max_k_iter = aggr[aggr['K']==max_k].shape[0]
    best_run = int(aggr[aggr['K']==max_k].obj.idxmin())
    print("   * Run {} had lowest objective with mode (n={:g}) K = {:g}.".format(best_run, max_k_iter, aggr.loc[best_run]['K']))

    store = pd.HDFStore(os.path.join(outdir,'nmf_output.h5'),'a')
    store["H"] = store["run{}/H".format(best_run)]
    store["W"] = store["run{}/W".format(best_run)]
    store["lam"] = store["run{}/lam".format(best_run)]
    store["Hraw"] = store["run{}/Hraw".format(best_run)]
    store["Wraw"] = store["run{}/Wraw".format(best_run)]
    store["markers"] = store["run{}/markers".format(best_run)]
    store["signatures"] = store["run{}/signatures".format(best_run)]
    store["log"] = store["run{}/log".format(best_run)]
    store["aggr"] = aggr
    store.close()

    # Plots
    print("   * Saving report plots to {}".format(outdir))

    H = pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "H")
    X = pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "X")
    signatures = pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "signatures")

    _ = k_dist(np.array(aggr.K, dtype=int))
    plt.savefig(os.path.join(outdir, "k_dist.pdf"), dpi=100, bbox_inches='tight')

    _ = marker_heatmap(X, signatures, H.sort_values('max_id').max_id)
    plt.savefig(os.path.join(outdir, "marker_heatmap.pdf"), dpi=100, bbox_inches='tight')

    if consensus_cluster_results:
        print("   * Computing consensus matrix")
        cmatrix, _ = consensus_cluster(os.path.join(outdir, 'nmf_output.h5'))
        f,d = consensus_matrix(cmatrix, n_clusters=max_k_iter)

        cmatrix.to_csv(os.path.join(outdir, 'consensus_matrix.tsv'), sep='\t')
        d.to_csv(os.path.join(outdir, 'consensus_assign.tsv'), sep='\t')
        plt.savefig(os.path.join(outdir, 'consensus_matrix.pdf'), dpi=100, bbox_inches='tight')

        store = pd.HDFStore(os.path.join(outdir,'nmf_output.h5'),'a')
        store['consensus'] = d
        store.close()
    # ------------------------------------------------------------------}

    if inplace:
        join_nmf_output_to_anndata(adata, os.path.join(outdir,'nmf_output.h5'))
        return
    else:
        W = pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "W")
        markers = pd.read_hdf(os.path.join(outdir,'nmf_output.h5'), "markers")
        return H,W,markers,signatures
