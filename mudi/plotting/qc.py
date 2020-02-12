import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

def adata_qc(adata, color_by='batch', figsize=(10,6)):
    """
    Plot AnnData QC
    ------------------
    Plot some useful QC plots.
    Requires that sc.compute_qc_metrics is run first.
    """
    fig, axes = plt.subplots(2,2,figsize=figsize)

    sc.pl.scatter(adata, x='n_counts', y='percent_mito',color=color_by,ax=axes[0,0], show=False)
    sc.pl.scatter(adata, x='n_counts', y='n_genes',color=color_by,ax=axes[0,1], show=False)

    sns.violinplot(data=adata.obs, y='log1p_n_genes_by_counts', x=color_by,ax=axes[1,0])
    sns.violinplot(data=adata.obs, y='percent_mito', x=color_by,ax=axes[1,1])

    plt.tight_layout()
