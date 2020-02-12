# mudi
A collection of python-based single-cell utilities and pipelines.

Author: Shankara Anand
Email: sanand@broadinstitute.edu

![](./mudi.jpg)

For more information about mudis, the amazing dog breed this package is named after, please visit https://en.wikipedia.org/wiki/Mudi.
---

### Requirements
  * python>=3.6
  * R>=3.5.1
    * Seurat>=3.0.0
    * scran>=1.10.1

#### Submodules

`git clone --recursive https://github.com/broadinstitute/mudi.git`

or if you have already cloned:

`git submodule update --init --recursive`

Then, `cd` into `mudi` and:

`pip install -e .`
---

### Usage

#### Raw Pre-processing of Count Data

```{python}
import mudi as md
from mudi.process import recipe

adata = recipe(
  list_of_10x_files,
  min_genes=200,
  min_cells=3,
  compute_doublets=True,
  remove_doublets=False,
  verbose=True
)

```

#### Cell-Type Identification

```python
from mudi.markers import build_marker_set, sub_cluster_and_rename

# Find marker genes, label clusters in each "groupby" group based on markers annotation
# Stored in adata.obs['cell_type']
scores, aggr, labels = build_marker_set(adata, heart_markers, groupby='louvain', key_added='cell_type' thresh=1e-2)

# Sub-clustering can be done on broader categories
# ex. this sub-clusters the "Myocardial Pericyte" cells and renames the newly created clusters
_ = sub_cluster_and_rename(adata, 'cell_type', ['Myocardial Pericyte'], heart_markers)
```
---

Examples and testing coming soon...

- [ ] Batch Integration
- [ ] Differential Expression
- [ ] Proportion Tests
- [ ] Plotting features
- [ ] QC plotting
- [ ] Integration with SCRINVEX
