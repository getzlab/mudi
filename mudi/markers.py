import gc
import pandas as pd
from tqdm import tqdm
import os
import scanpy as sc
import scanpy.external as sce
import sys

from .utils import aggr_markers

# ---------------------------------
# Cell-type Calling
# ---------------------------------
def labeler(d,x):
    try:
        return d[x]
    except:
        return 'n/a'

def build_marker_set(adata, markers_df, metric='sum', thresh=1e-5, **kwargs):
    """
    Build Marker set.
    ---------------------
    This takes in a spy object, a dictionary mapping from
    cell-types to markers of interest.

    Inputs:
        - adata: scanpy anndata object
        - markers_df: markers dataframe
        - metric: metric for assigning cell-type
        - thresh: threshold for markers to consider
        - kwargs: inputs for aggr_markers

    Outputs:
        - scores: dataframe of scores that have the sum z-scores for each cell-type
            for a given cluster; clusters with no marker genes are assigned -1
        - aggr: aggregate list of marker genes by cell type; includes pvals, z-scores,
            initial cluster assignment, gene name, what cell-type the marker comes from,
            and what cell type it was ultimately labeled
        - labels: dictionary mapping from cluster to labeled cell-type
    """
    markers = aggr_markers(adata, **kwargs)
    markers = markers[markers['pvals_adj']<thresh]

    d = {}
    dfs = list()
    for cell_type in set(markers_df['Cell-Type']):
        filt = markers[markers.names.isin(markers_df[markers_df['Cell-Type']==cell_type]['Gene'].values)]
        filt['cell_type'] = cell_type #heck
        dfs.append(filt)

        if metric == 'sum':
            d[cell_type] = dict(filt.groupby('cluster').sum().scores)
        elif metric == 'mean':
            d[cell_type] = dict(filt.groupby('cluster').mean().scores)
        else:
            raise ValueError("Not yet implemented.")

    # Compile scores
    scores = pd.DataFrame.from_dict(d).fillna(-1)

    # Combine all markers for each cell type
    aggr = pd.concat(dfs)

    # Labels
    labels = {key:scores.loc[[key],:].T.sort_values(key,ascending=False).index[0] for key in list(scores.index)}

    # Add labeling
    aggr['label'] = aggr['cluster'].apply(lambda x: labels[x])

    return scores, aggr, labels

def sub_cluster_and_rename(adata, group, group_vars, markers, res=1, method='t-test_overestim_var', tweaks={}):
    """
    Subclusters & Rename
    ---------------------
    Takes in an adata object and information about grouping/variables, does subclustering, and
    reassigns variables.

    Inputs:
        - adata: scanpy anndata object
        - group: grouping to do sub-clustering on (ex. louvain)
        - group_vars: sub-groups within grouping to do sub-clustering on (ex. ['1','2'] for louvain clusters 1 & 2)
        - markers: marker list
        - res: resolution of louvain clustering
        - method: DE method for marker assignment
        - tweaks: dictionary of specific edits to make to labeling {'tcell,1': 'tcell_cd4', 'tcell,2':'tcell_cd8'}

    Outputs:
        - genes: marker genes used to label clusters
    """
    new_name = group+'_R'

    # sub-cluster
    sc.tl.louvain(adata, resolution=res, restrict_to = (group, group_vars), key_added=new_name)

    # compute umap + rank-genes
    sc.pl.umap(adata, color=[new_name,])
    sc.tl.rank_genes_groups(adata, new_name, method=method)

    # find markers
    scores, genes, labels = build_marker_set(adata, markers)

    # custom update if needed
    labels.update(tweaks)

    # reassign
    adata.obs[new_name] = adata.obs[new_name].apply(lambda x: labeler(labels,x))
    sc.pl.umap(adata, color=[new_name])

    return genes
