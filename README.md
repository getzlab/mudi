# mudi
A collection of python-based single-cell utilities and pipelines.

###### **Author**: Shankara Anand

###### **Email**: sanand@broadinstitute.edu

![](./mudi.jpg)

_For more information about mudis, the amazing dog breed this package is named after, please visit_ https://en.wikipedia.org/wiki/Mudi.

---

### Requirements
  * python>=3.6
  * R>=3.5.1
    * Seurat>=3.0.0
    * scran>=1.10.1

### Installation
```{bash}
git clone https://github.com/broadinstitute/mudi.git

cd mudi
pip install -e .
```

_Support for PIP coming soon._

### Package Details

###### Includes
* wrappers for sequence level QC (ie `FastQC`)
* pre-processing recipes for raw count data with qc & filtering helpers
* wrappers for single-cell normalization (i.e. `scran`)
* wrapper for `signatureanalyzer` (bayesian NMF)
* wrapper for `gprofiler`

###### Examples: see `./examples`
* `normalization`: general processing & normalization
* `cell_types`: assigning cell-types
* `ardnmf`: running ard-nmf via `signatureanalyzer`

### Example Usage

```{python}
import mudi as md

adata = md.recipe(
  list_of_10x_files,
  min_genes=200,
  min_cells=3,
  compute_doublets=True,
  remove_doublets=False,
  verbose=True,
  norm='scran'
)
```

Examples and testing coming soon...

- [ ] Batch Integration
- [ ] Differential Expression
- [ ] Proportion Tests
- [ ] Plotting features
- [ ] QC plotting
- [ ] Integration with SCRINVEX
