import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def pw_enrichment(
    df,
    xax='ES',
    padj_idx='padj',
    ax=None,
    padj_thresh=1e-2,
    figsize=(6,6),
    size_min=5,
    labelsize=9,
    cmap='viridis',
    title=None,
    cbar_ax_title='Adj P-value',
    sort_by=None,
    sort_by_ascending=True,
    xlim=None,
    n_to_plot=15,
    size_multi=1.25
    ):
    """
    Plot GSEA results
    -----------------------
    Args:
        * df: pd.DataFrame
        * xax: value to plot on x-axis (normally enrichment score)
        * padj_idx: what column in dataframe to use for plotting
        * ax: matplotli.pyplot.axis
        * padj_thresh: threshold value to plot for `padj`; pathways above this threshold will
            not be plotted
        * figsize: tuple of figure size
        * size_min: minimum gene overlap in pathway to consider for plotting
        * labelsize: int size of pathway labels
        * cmap: colormap for p-value coloring of circles
        * title: title of plot
        * cbar_ax_title: color-bar axis title
        * sort_by: how to sort entries for plotting
            * default is to use the enrichment (or normalized) enrichment score; if this is
              selected, then it wil take the top highest and lowest "n_to_plot" pathways;
              if sorted by the padj_idx, then it will take 2*n_to_plot index
        * sort_by_ascending: how to sort the "sort_by" option
        * xlim: custom x-limits; else automatically determined
        * n_to_plot: number of pathways to plot
        * size_multi: multiplier of circle-sizes

    Returns:
        * pd.DataFrame of data plotted
    """
    from matplotlib import rcParams
    from matplotlib import colors

    if ax == None:
        fig,ax = plt.subplots(figsize=figsize)

    if sort_by is None:
        sort_by = padj_idx
        sort_by_ascending = True

    data_to_plot = df[(df[padj_idx] < padj_thresh) & (df['size']>size_min)].copy()
    data_to_plot = data_to_plot.sort_values(sort_by, ascending=sort_by_ascending)

    if sort_by == padj_idx:
        # Plot top by significance
        data_to_plot = data_to_plot.head(n_to_plot*2)
    elif sort_by in ("NES", "ES"):
        # Plot top / bottom by ES
        data_to_plot = pd.concat((data_to_plot.head(n_to_plot), data_to_plot.tail(n_to_plot))).drop_duplicates(subset='pathway')

    xlim_bound = max(np.max(data_to_plot[xax]), np.abs(np.min(data_to_plot[xax])))

    norm = colors.LogNorm(data_to_plot[padj_idx].min(), data_to_plot[padj_idx].max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    rcParams.update({'font.size': 14, 'font.weight': 'normal'})

    ax.set_axisbelow(True)

    if data_to_plot.shape[0]==0:
        print("   * no significant pathways")
        return
    else:
        print("   * plotting {} pathways".format(data_to_plot.shape[0]))

    path = ax.scatter(
        x=xax,
        y='pathway',
        c=padj_idx,
        cmap=cmap,
        norm=norm,
        data=data_to_plot,
        linewidth=1,
        edgecolor="black",
        s=[(i+10)**size_multi for i in data_to_plot['size']],
        zorder=10
    )

    # Plot Bars
    bar_colors = ['pink' if x > 0 else 'lightblue' for x in data_to_plot[xax]]
    data_to_plot.plot(
        kind='barh',
        y=xax,
        x='pathway',
        legend=None,
        ax=ax,
        alpha=0.5,
        color=bar_colors,
        width=0.7,
        zorder=1,
    )


    # Lower p-value on top
    ax.invert_yaxis()
    ax.tick_params(axis='y', which='major', labelsize=labelsize)
    ax.set_ylabel('')
    ax.set_xlabel(xax, fontsize=14, fontweight='normal')
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

    if xlim is None:
        ax.set_xlim(-xlim_bound-.1, xlim_bound+.1)
    else:
        ax.set_xlim(xlim)

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # -------------------------
    # p-value Colorbar
    # -------------------------
    cbaxes = ax.figure.add_axes([0.8, 0.125, 0.03, 0.2])
    cbar = ax.figure.colorbar(sm, shrink=0.4, cax=cbaxes)
    cbar.ax.invert_yaxis()

    cbar.ax.set_title(cbar_ax_title, fontweight='normal', ha='center', x=2, y=1.15, fontsize=12)
    cbar.ax.tick_params(axis='y', which='major', labelsize=8)
    cbar.ax.tick_params(axis='y', which='minor', labelsize=8)

    # -------------------------
    # Size Legend
    # -------------------------
    min_bound = int(np.floor(min(data_to_plot['size']) / 10.0)) * 10
    max_bound = int(np.ceil(max(data_to_plot['size']) / 10.0)) * 10
    med_bound = int((min_bound + max_bound) / 2)

    l1 = ax.scatter([],[], s=(min_bound+10)**size_multi, edgecolors='none', color='black')
    l2 = ax.scatter([],[], s=(med_bound+10)**size_multi, edgecolors='none', color='black')
    l3 = ax.scatter([],[], s=(max_bound+10)**size_multi, edgecolors='none', color='black')
    labels = [min_bound, med_bound, max_bound]

    leg = ax.legend(
        [l1, l2, l3],
        labels,
        ncol=1,
        frameon=False,
        fontsize=12,
        handlelength=2.5,
        loc = 'upper right',
        labelspacing = 2,
        bbox_to_anchor=(1.45,1),
        title='Gene Overlap',
        scatterpoints = 1,
        facecolor='black'
    )

    ax.axvline(0, color='black', linewidth=1, alpha=0.5)
    ax.set_title(title)

    return data_to_plot
