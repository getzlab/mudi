import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import numpy as np

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

def _plot_scatter(ad, x='log1p_total_counts', y='n_genes', c=None, ax=None, vmin=None, vmax=None, **kwargs):
    if ax is None:
        fix,ax= plt.subplots(figsize=(5,5))

    ad_df = ad.obs

    if c is not None:
        ad_df = ad_df.sort_values(c, ascending=False)
        _color = ad_df[c]

    ax.scatter(ad_df[x], ad_df[y], c=_color, vmin=vmin, vmax=vmax, **kwargs)
    ax.grid(False)

def adata_qc_grid(adata, batch, c='percent_mito', vmin=0, vmax=1, s=10, alpha=0.6, **kwargs):
    """
    Adata QC Grid.
    """
    samples_to_plot = adata.obs[batch].unique()
    n_samples = len(samples_to_plot)

    print("   * {} batches to plot".format(str(n_samples)))

    if len(samples_to_plot) <= 2:
        return
    else:
        fig,axes = plt.subplots(
            int(np.ceil(n_samples/2)),
            2,
            figsize=(8, 4 * int(np.ceil(n_samples/2))),
            sharex=True,
            sharey=True
        )

        c=0
        for i in range(int(np.ceil(n_samples/2))):
            for j in range(2):
                try:
                    _plot_scatter(
                        adata[adata.obs['sample']==samples_to_plot[c]],
                        c='percent_mito',
                        ax=axes[j,i],
                        vmin=vmin,
                        vmax=vmax,
                        s=s,
                        alpha=alpha,
                        rasterized=True,
                        **kwargs
                    )
                    axes[j,i].grid(False)
                    axes[j,i].set_xlabel("log UMI", fontsize=14)
                    axes[j,i].set_ylabel("# Genes", fontsize=14)
                    axes[j,i].set_title(samples_to_plot[c], fontsize=16)
                except:
                    axes[j,i].axis('off')
                c+=1

    plt.tight_layout()
    im = plt.gca().get_children()[0]
    cax = fig.add_axes([1,.1,0.03,.8])
    cax.set_title("% Mito")
    fig.colorbar(im, cax=cax)
