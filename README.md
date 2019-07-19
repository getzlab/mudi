# mudi
A collection of single-cell utilities and pipelines for the Getz Lab.


![](./mudi.jpg)

---

### Requirements
  * python>=3.6
  * R>=3.5.1
    * Seurat>=3.0.0
    * scran>=1.10.1

---

### Usage

#### General Processing

```python

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
from mudi.markers import build_marker_set
from mudi.markers import labeler

# Find marker genes
sc.tl.rank_genes_groups(adata, groupby='louvain')

# Build marker set
scores, aggr, labels = build_marker_set(adata, heart_markers, thresh=1e-2)

# Assign cell-type
adata.obs['cell_type'] = adata.obs['louvain'].apply(lambda x: labeler(labels,x))

```
---

Examples and testing coming soon...

- [ ] Batch Integration
- [ ] Differential Expression
- [ ] Proportion Tests
- [ ] Plotting features
- [ ] QC plotting
- [ ] Integration with SCRINVEX
