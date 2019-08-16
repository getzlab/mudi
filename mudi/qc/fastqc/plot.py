import os
from io import StringIO
import subprocess
import zipfile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_fq_metric(df, hue='sample', ax=None, title=None, legend=True, **kwargs):
    """
    Takes output from the aggr_metric function and plots.
    """
    col_names = list(df)

    if ax is None:
        fig,ax = plt.subplots(figsize=(10,8))

    sns.lineplot(data=df, x=col_names[0], y=col_names[1], ax=ax, sort=False, hue=hue)

    if legend:
        ax.legend(**kwargs)
    else:
        ax.get_legend().remove()

    ax.set_title(title)

    for tick in ax.get_xticklabels():
        tick.set_rotation(90)
