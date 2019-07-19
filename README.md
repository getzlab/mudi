# mudi
A collection of single-cell utilities and pipelines for the Getz Lab.

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
---

Examples and testing coming soon...

- [ ] Batch Integration
- [ ] Differential Expression
- [ ] Proportion Tests
- [ ] Plotting features
- [ ] QC plotting
- [ ] Integration with SCRINVEX
