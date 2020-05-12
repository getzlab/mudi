import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import numpy as np

def volcano_plot(
    de,
    thresh=20,
    ax=None,
    xlim=None,
    yax='qval',
    xax='logFC',
    xlabel=None,
    ylabel=None,
    gene_id='genes',
    ptm_id=None,
    filter_noise=False,
    fix_extremes=True,
    arrow=False,
    label_percentile=99.5,
    shuffle=True,
    gene_fontsize=10,
    only_plot_gene_symbol=True,
    genes_to_label=None
    ):
    """
    Differential Expression Volcano Plot
    ---------------------
    Args:
        * de: pd.DataFrame of differentail expression results
        * thresh: adj. pval threshold for differential expression (color)
        * ax: matplotlib.axis
        * xlim: x-limits for the volcano plot; automatically determined / centered
            if not-provided, takes symmetric ceil of x-values
        * yax: column in dataframe for y-axis ('qval')
        * xax: column in dataframe for x-axis (usually log-fold-change)
        * xlabel: xlabel
        * ylabel: ylabel
        * gene_id: column in dataframe for gene-labels
        * filter_noise: remove noise
        * fix_extremes: fix-extremes
        * arrow: add arrows to each gene label
        * label_percentile: percentile of genes on each neg/pos lfc to label
        * shuffle: shuffle points
        * gene_fontsize: size of gene-labels
        * genes_to_label: list of genes you can provide to label with

    Returns:
        * returns fig
    """
    from numpy import inf
    from adjustText import adjust_text

    # Subset relevant Columns
    if ptm_id in de:
        de = de.loc[:,[xax, yax, gene_id, ptm_id, 'variableSites']]
    else:
        de = de.loc[:,[xax, yax, gene_id]]

    if xlim is None:
        _x = np.ceil(np.max(np.abs(de[xax])))
        xlim = (-_x, _x)

    # Shuffle columns
    if shuffle:
        de = de.sample(frac=1)

    # Create -log10q-value
    de['-logQ'] = -np.log10(de[yax])

    # Get percentile threshold for each
    de_up = de[de[xax]>0]
    de_down = de[de[xax]<0]

    # Filter extremes
    if fix_extremes:
        def _fix(x,v):
            if x==inf:
                return int(v)
            else:
                return x

        val = int(np.max(de[de['-logQ'] != inf]['-logQ'])+1)
        de['-logQ'] = de['-logQ'].apply(lambda x: _fix(x,val))

    if filter_noise:
        de[(de[xax] < -20) & (de[yax] < 1e-3)]

        noisy_down_genes = de[(de[xax]<-5) & (de[yax]>1e-3)][gene_id]
        noisy_up_genes = de[(de[xax]>5) & (de[yax]>1e-3)][gene_id]
        noisy_genes = set(noisy_up_genes) | set(noisy_down_genes)

        de = de[~de[gene_id].isin(noisy_genes)]

    lowqval_de = de.loc[de['-logQ'] > thresh]
    other_de = de.loc[de['-logQ'] < thresh]

    if ax is None:
        fig, ax = plt.subplots(figsize=(6,6))

    # Below Threshold
    sns.regplot(
        other_de[xax],
        np.abs(other_de['-logQ']),
        fit_reg=False,
        scatter_kws={'s':12, 'alpha':0.5, 'color':'k','rasterized':True},
        ax=ax,
    )

    sns.regplot(
        lowqval_de[xax],
        np.abs(lowqval_de['-logQ']),
        fit_reg=False,
        scatter_kws={'s':12, 'alpha':0.5,'color':'r','rasterized':True},
        ax=ax,
    )

    # Axis labels
    if xax=='diff':
        xax_label = "Difference in Means"
    else:
        xax_label = xax.capitalize()

    ax.set_xlabel(xax_label, fontsize=14)
    ax.set_ylabel("-log10 {}".format(yax.capitalize()), fontsize=14)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, min_n_ticks=5, nbins=4))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, min_n_ticks=5, nbins=4))

    # Label names and positions
    _up_thresh = np.percentile(de_up['-logQ'], label_percentile)
    _down_thresh = np.percentile(de_down['-logQ'], label_percentile)

    x = lowqval_de[xax].values
    y = lowqval_de['-logQ'].values
    labels = lowqval_de[gene_id].values

    to_keep = []
    for idx,x_val in enumerate(x):
        if x[idx] < 0 and y[idx] > _down_thresh:
            to_keep.append(idx)
        elif x[idx] > 0 and y[idx] > _up_thresh:
            to_keep.append(idx)

    x = x[to_keep]
    y = y[to_keep]
    labels = labels[to_keep]

    # Additional genes to label
    if genes_to_label is not None:
        _de_to_label = de[de[gene_id].isin(genes_to_label)]
        x = np.concatenate((x,_de_to_label[xax].values))
        y = np.concatenate((y,_de_to_label['-logQ'].values))
        labels = np.concatenate((labels,_de_to_label[gene_id].values))

    ax.axvline(0, color='black', alpha=0.2, rasterized=True)
    ax.axhline(0, color='black', alpha=0.2, rasterized=True)

    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=14)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=14)

    # Set xlabel
    if xlim is not None:
        ax.set_xlim(xlim)
    else:
        xlim = ax.get_xlim()

    texts = list()
    for i,txt in enumerate(labels):
        if x[i] > xlim[0] and x[i] < xlim[1]:
            texts.append(ax.text(x[i], y[i], txt, ha='center', va='center', fontsize=gene_fontsize))
    if arrow:
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='blue'))
    else:
        adjust_text(texts)

    return ax
