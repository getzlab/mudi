import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc

def plot_adata_qc(adata, color_by='batch', figsize=(10,6)):
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

def volcano_plot(de, thresh=20, ax=None):
    """
    Plot a volcano plot of differentially expressed genes.
    """
    de['-logQ'] = -np.log10(de['FDR'])
    lowqval_de = de.loc[de['-logQ'] > thresh]
    other_de = de.loc[de['-logQ'] < thresh]

    if ax is None:
        fig, ax = plt.subplots()

    sns.regplot(other_de['coef'], other_de['-logQ'], fit_reg=False, scatter_kws={'s':12, 'alpha':0.5,'color':'b'},ax=ax)
    sns.regplot(lowqval_de['coef'], lowqval_de['-logQ'], fit_reg=False, scatter_kws={'s':12, 'alpha':0.5,'color':'r'},ax=ax)
    ax.set_xlabel("$\log_2$FC", fontsize=20)
    ax.set_ylabel("-$\log_{10}$Q", fontsize=20)
    ax.tick_params(labelsize=15)

    # Label names and positions
    lowqval_de = lowqval_de.dropna()
    x = [i for i in lowqval_de['coef']]
    y = [i for i in lowqval_de['-logQ']]
    labels = lowqval_de['primerid']

    # Show only some labels to avoid overcrowding the figure
    to_remove = np.where([i < thresh for i in y])[0]
    labels = ["" if i in to_remove else lab for i,lab in enumerate(labels) ]

    for i,txt in enumerate(labels):
        ax.annotate(txt, (x[i], y[i]))
    return ax
