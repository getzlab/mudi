import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

def plot_barplot_a_by_b(hm, palette=None, ax=None, figsize=(5,2), width=0.8, stacked=True):
    """
    Plot the barplot a by b.
    ---------------------
    Args:
        * hm: pd.DataFrame
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if palette is not None:
        c = [palette[x] for x in hm.columns]
    else:
        c = None

    hm.plot(
        kind='barh',
        ax=ax,
        stacked=stacked,
        width=width,
        edgecolor='black',
        rasterized=True,
        color=c
    )

    ax.legend(frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylim(ax.get_ylim()[0]+0.25, ax.get_ylim()[1]-0.25)
    ax.set_xlabel("Subtype Proportion", fontsize=14)
    ax.set_ylabel("")
    ax.grid(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
