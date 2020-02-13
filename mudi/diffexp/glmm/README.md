### Differential Expression in Single-cell

**Author**: Shankara Anand

**Date**: 02/12/20

_**Note**: this is intended for 10x Data which is assumed to be subject to overdispersion. A natural distribution to fit expression profiles to is negative-binomial._

Docker: `gcr.io/broad-cga-sanand-gtex/r36:latest`

---

### Pipeline

####  1. Prepare Inputs

We have the following needs before running the computation for differential expression:
* Take an AnnData object and save its count matrices
* Save the relevant AnnData object metadata to use as covariates when fitting GLMs
* Save individual counts file for each gene

This is all accessible using the following function:

```{python}
import mudi as md
import scanpy as sc

# Load AnnData
adata = sc.read("anndatas/adata.h5")

# Args
diff_exp_dir = "./de/diffexp"
var_to_groupby = "molecular_subtype"
effects_to_include = ["sample","age","sex"]

md.de.prep_inputs(adata, diff_exp_dir, var_to_groupby, effects_to_include)
```

This will result in the following directory structure:
```
./de/diffexp/inputs
  * genes
    * GENE.csv ... etc.
  * meta.csv
  * raw_counts.csv
  * raw_counts.parquet
```

#### 2. Estimating Dispersions

If we take a page from Bulk RNA-Seq differential expression, a core step before fitting negative binomial GLMs for each gene is estimating the appropriate dispersions. A GLM for negative binomial distribution is only defined with a stated dispersion. The tool we will be using in later steps is `lme4` - our statistical hammer to fit these models - can estimate dispersions per gene fit. However, this can by noisy and not converge with sparse-single-cell data. Thus we will use, `edgeR`' to estimate dispersions using `edgeR::estimateDisp`.

A script is provided to compute this:

```
Rscript compute_disp.R -i inputs/raw_counts.csv \
                       -m ./de/diffexp/inputs/meta.csv \
                       -o ./de/diffexp/inputs/dispersions.tsv
```

_**Note**: this is highly memory intensive and might take a few hours given the size of your input data._

_**Note**: run with_ `export OPENBLAS_NUM_THREADS=1` _to prevent overuse of unncessary threads._

#### 3. Running lme4 Tests for Each Gene

Now, we need to fit a negative binomial model for each gene, specify random and fixed effects, and test for a group of interest to generate our differential expression results. This would take a long time to loop through given each gene requires `n_groups`+1 (null) models to fit. This roughly takes minimum 5-10 minutes for 20K genes which could take up to 2.22 days minimum. Thes, we use `canine` to dynamically split up google compute nodes, `slurm` to dispatch the jobs and rapidly run all of these analyses.

This is all packaged together in mudi as follows where you can specify the model of your choosing:

```{python}
import mudi as md

md.de.glmm.dispatch(
  "./de/diffexp/inputs/meta.csv",
  "./de/diffexp/inputs/genes",
  "./de/diffexp/inputs/dispersions.tsv",
  "Y ~ 1 + (1|sample) + age + sex + offset(loggenesxcounts)",
  "sa_transfer",
  n_nodes=1,
  worker_type='n1-highmem-64'
)
```

* `n_nodes`: the number of nodes `canine` is allowed to spin up
* `worker_type`: gcp worker type

The model formula must be specified with the syntax `Y ~ `; aside from that, you are only limited by what covariates you would like to include. **This is only the null model.** This means that any groups that you are using to do your one vs. rest differential expression should NOT be represented here, only covariates (random and fixed effects) you would like to model.

The transfer bucket must be accessible by your default account.
