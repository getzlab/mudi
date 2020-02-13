import canine
import os
import numpy as np
import glob
import pkg_resources
from typing import Union

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
