import numpy as np
import sys
import os
import pandas as pd
import glob
import subprocess
import scanpy as sc
from agutil import parallel

def prep_inputs(adata, diffexp_dir, meta_vars=None, genes_to_use=None, n_threads=15):
    """
    Prepare Inputs for Differential Expression
    -------------------------------
    Args:
        * adata
        * diffexp_dir
        * meta_vars: list of variables in metadata

    Returns:
        None
    """
    if meta_vars is None:
        meta_vars = list(adata.obs)
    if genes_to_use is None:
        genes_to_use = adata.var_names

    for v in meta_vars:
        assert v in list(adata.obs), print("{} not in adata.obs.".format(v))

    os.makedirs(diffexp_dir, exist_ok=True)
    os.makedirs(os.path.join(diffexp_dir, "inputs"), exist_ok=True)
    os.makedirs(os.path.join(diffexp_dir, "inputs", "genes"), exist_ok=True)

    # --------------------------
    # Create counts dataframe
    # --------------------------
    counts_df = pd.DataFrame(
        adata.layers['counts'].todense(),
        index=adata.obs_names,
        columns=adata.var_names,
        dtype=int
    )[genes_to_use]

    print("{} barcodes".format(counts_df.shape[1]))
    print("{} genes".format(counts_df.shape[0]))

    counts_df.T.to_csv(os.path.join(diffexp_dir, "inputs", "raw_counts.csv"))
    counts_df.T.to_parquet(os.path.join(diffexp_dir, "inputs", "raw_counts.parquet"))

    adata.obs[meta_vars].rename(
        columns={'log1p_n_genes_by_counts':'loggenesxcounts'}
        ).to_csv(os.path.join(diffexp_dir, "inputs", "meta.csv"))

    @parallel.parallelize(maximum=n_threads)
    def save_slice(gene_i):
        gene = gene_i.replace("/",'_')
        try:
            counts_df[gene_i].to_csv(os.path.join(diffexp_dir, "inputs", "genes", "{}.csv".format(gene)))
        except:
            print("Error with {}".format(gene))

    print("Saving {} genes...".format(counts_df.shape[0]))
    _ = [x for x in save_slice(counts_df)]
