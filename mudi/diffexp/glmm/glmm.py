import numpy as np
import sys
import os
import pandas as pd
import glob
import subprocess
import scanpy as sc
from agutil import parallel
from anndata import AnnData
from typing import Union
from tqdm import tqdm
import pkg_resources
from typing import Union

import rpy2
from rpy2.robjects.packages import importr

"""
Author: Francois Aguet
https://github.com/broadinstitute/tensorqtl/blob/master/tensorqtl/rfunc.py
"""
def qvalue(p, lambda_qvalue=None):
    """Wrapper for qvalue::qvalue"""
    qvalue = importr("qvalue")
    rp = rpy2.robjects.vectors.FloatVector(p)
    if lambda_qvalue is None:
        q = qvalue.qvalue(rp)
    else:
        if not isinstance(lambda_qvalue, Iterable):
            lambda_qvalue = [lambda_qvalue]
        rlambda = rpy2.robjects.vectors.FloatVector(lambda_qvalue)
        q = qvalue.qvalue(rp, **{'lambda':rlambda})
    qval = np.array(q.rx2('qvalues'))
    pi0 = np.array(q.rx2('pi0'))[0]
    return qval, pi0

"""
GLMM Pipeline Functions
------------
* prep_inputs --> prepare inputs from AnnData for lme4 tests
* dispatch --> dispatch jobs on Canine
* compile_de_result --> compiles result from a single diffexp canine output
* compile_all_results --> compiles all outputs

"""
def prep_inputs(
    adata: AnnData,
    diffexp_dir: str,
    groupby: str,
    meta_vars: Union[None,list] = None,
    genes_to_use: Union[None, list] = None,
    n_threads: int = 15
    ):
    """
    Prepare Inputs for Differential Expression
    -------------------------------
    Args:
        * adata: AnnData file
        * diffexp_dir: output directory to write out inptus for lme4 to
        * groupby: variable to groupy
        * meta_vars: list of variables in metadata
        * genes_to_use: list of genes to use for differential expression tests
        * n_threads: nthreads for saving gene counts

    Returns:
        None
    """
    if meta_vars is None:
        meta_vars = list(adata.obs)
    if genes_to_use is None:
        genes_to_use = adata.var_names

    for v in meta_vars:
        assert v in list(adata.obs), print("{} not in adata.obs.".format(v))

    assert groupby in meta_vars, "{} not in metavars".format(groupby)

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

    meta_df = adata.obs[meta_vars].rename(
        columns={'log1p_n_genes_by_counts':'loggenesxcounts'}
    )

    dummy_df = pd.get_dummies(meta_df[groupby])
    dummy_df = dummy_df.rename(columns={x:"groupby_"+x for x in dummy_df})
    meta_df.join(dummy_df).to_csv(os.path.join(diffexp_dir, "inputs", "meta.csv"))

    @parallel.parallelize(maximum=n_threads)
    def save_slice(gene_i):
        gene = gene_i.replace("/",'_')
        try:
            counts_df[gene_i].to_csv(os.path.join(diffexp_dir, "inputs", "genes", "{}.csv".format(gene)))
        except:
            print("Error with {}".format(gene))

    print("Saving {} genes...".format(counts_df.shape[0]))
    _ = [x for x in save_slice(counts_df)]

def dispatch(
    meta: str,
    genes_dir: str,
    dispersions: str,
    model: str,
    transfer_bucket: str,
    n_nodes: int = 1,
    worker_type: str = "n1-highmem-4",
    verbose:bool = False,
    ):
    """
    Dispatch.
    ----------------------
    Args:
        * meta: path to metadata dataframe (barcodes x covariates) processed by prep.py
        * genes_dir: directory of all genes
        * dispersions:
        * model:
        * transfer_bucket:
        * n_nodes:
        * worker_type

    """
    meta = os.path.abspath(meta)
    genes_dir = os.path.abspath(genes_dir)
    dispersions = os.path.abspath(dispersions)

    genes = [x.split(".csv")[0] for x in os.listdir(genes_dir)]
    print("Loaded {} genes.".format(len(genes)))

    conf = dict()

    conf["name"] = "lme4_diffexp"
    conf["backend"] = {
        "type": "TransientGCP",
        "compute_zone": "us-central1-a",
        "controller_type": "n1-highmem-16",
        "secondary_disk_size": 100,
        "worker_type": worker_type,
        "max_node_count": n_nodes
    }

    conf["localization"] = {
        "staging_dir": "/mnt/disks/sec/canine",
        "transfer_bucket": transfer_bucket
    }

    conf["inputs"] = dict()
    conf["inputs"]["SCRIPT"] = pkg_resources.resource_filename('mudi', 'diffexp/glmm/lme4_lrt.R')
    conf["inputs"]["METADATA"] = meta
    conf["inputs"]["DISP"] = dispersions
    conf["inputs"]["COUNTS"] = [os.path.join(genes_dir, "{}.csv".format(gene)) for gene in genes]

    conf["outputs"] = dict()
    conf["outputs"]["output"] = "*.tsv"

    conf["resources"] = {
        "cpus-per-task": 1,
        "mem-per-cpu": 6144
    }

    conf["script"] = ["set -e -o pipefail"]
    conf["script"].extend([
        "sudo docker run -v $CANINE_ROOT:$CANINE_ROOT --rm \
                    --cpus $SLURM_CPUS_PER_TASK \
                    --memory $(expr $SLURM_CPUS_PER_TASK '*' $SLURM_MEM_PER_CPU)MB \
                    -t gcr.io/broad-cga-sanand-gtex/r36:latest Rscript --vanilla $SCRIPT -i $COUNTS -m $METADATA -d $DISP -f '{}' -o $CANINE_JOB_ROOT".format(model)
    ])

    if verbose:
        print("Running:")
        print(conf["script"])

    orch = canine.Orchestrator(config=conf)
    R = orch.run_pipeline()

def compile_de_result(out_dir: str):
    """
    Compile gene-level differential expression result
    ------------------------
    Args:
        * out_dir: ingividual output directory per gene
    """
    summary_tsv = glob.glob(os.path.join(out_dir, "*.summary.tsv"))[0]
    lrt_tsv = glob.glob(os.path.join(out_dir, "*.lrt.tsv"))[0]

    gene = summary_tsv.split("/")[-1].split(".")[0]

    summary_df = pd.read_csv(summary_tsv, sep='\t')

    pval_df = summary_df.join(pd.read_csv(lrt_tsv, sep='\t')).dropna().rename(columns={"Pr(>Chisq)":"lrt_p_val"}).drop(columns="celltype")
    pval_df = pval_df.rename(index={x:x.split("groupby_")[-1] for x in pval_df.index})
    pval_df['gene'] = gene

    summary_df['gene'] = gene
    summary_df = summary_df.rename(columns={'celltype':'groupby'})
    summary_df['groupby'] = summary_df['groupby'].apply(lambda x: x.split("groupby_")[-1])

    return pval_df, summary_df

def compile_all_results(output_dirs: list, canine_dir_name: str = "canine_output"):
    """
    Compile all results.
    ------------------------
    Args:
        * output_dirs: list of directories where canine_outputs are located
    """
    pvals = list()
    sums = list()

    for output_dir in output_dirs:
        file_list = glob.glob(os.path.join(output_dir, canine_dir_name, "**/output"))

        for x in tqdm(file_list, desc=output_dir):
            p,s = compile_de_result(x)
            pvals.append(p)
            sums.append(s)

    # ------------------------
    # Q-values for lrt-pvalues
    # ------------------------
    p_df = pd.concat(pvals)
    _p_df = list()

    for group in np.unique(p_df.index):
        group_df = p_df.loc[group].sort_values('lrt_p_val')
        group_df['qval'], pi0 = qvalue(np.array(group_df['lrt_p_val']))
        _p_df.append(group_df)

    p_df = pd.concat(_p_df)
    p_df = p_df.loc[:,['gene','lrt_p_val','qval','coef','stderr','z','p_val','Df','AIC','BIC','logLik','deviance','Chisq','Chi Df']].rename(columns={'coef':'logFC'})

    # ------------------------
    # Q-values for all covariates
    # ------------------------
    s_df = pd.concat(sums)
    _s_df = list()

    vals = s_df.iloc[0].loc[['groupby','gene']]
    labs = s_df[(s_df['groupby']==vals[0]) & (s_df['gene']==vals[1])].index


    for idx,lab in enumerate(labs):
        group_df = s_df.iloc[idx::4,:].sort_values('p_val')
        group_df['label'] = lab
        group_df['qval'], pi0 = md.de.glmm.qvalue(np.array(group_df['p_val']))

        _s_df.append(group_df)

    s_df = pd.concat(_s_df).set_index('label')
    s_df = s_df.loc[:,["gene","groupby","coef","p_val","qval","z","stderr"]].rename(columns={"coef":"logFC"})

    return p_df, s_df
