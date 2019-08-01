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

def build_marker_set(adata, markers_df, groupby='louvain', metric='sum', thresh=1e-5, key_added='cell_type', cell_type_idx='Cell-Type', gene_idx='Gene', tweaks=None, **kwargs):
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
        - key_added: key to add to adata object
        - cell_type_idx: markers dataframe cell_type ID
        - gene_id: markers dataframe gene ID
        - tweaks: specific changes to label dict
        - **kwargs: for sc.tl.rank_genes_groups call - any specific paramters for selecting marker genes

    Outputs:
        - scores: dataframe of scores that have the sum z-scores for each cell-type
            for a given cluster; clusters with no marker genes are assigned -1
        - aggr: aggregate list of marker genes by cell type; includes pvals, z-scores,
            initial cluster assignment, gene name, what cell-type the marker comes from,
            and what cell type it was ultimately labeled
        - labels: dictionary mapping from cluster to labeled cell-type
    """
    if tweaks is None:
        tweaks = {}

    assert groupby in list(adata.obs), 'Please ensure {} is a field in AnnData object.'.format(groupby)

    sc.tl.rank_genes_groups(adata, groupby=groupby, **kwargs)

    markers = aggr_markers(adata)
    markers = markers[markers['pvals_adj']<thresh]

    d = {}
    dfs = list()
    for cell_type in set(markers_df[cell_type_idx]):
        filt = markers[markers.names.isin(markers_df[markers_df[cell_type_idx]==cell_type][gene_idx].values)]
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
    labels.update(tweaks)

    # Add labeling
    aggr['label'] = aggr['cluster'].apply(lambda x: labels[x])
    adata.obs[key_added] = adata.obs[groupby].apply(lambda x: labeler(labels, x))

    return scores, aggr, labels

def sub_cluster_and_rename(adata, group, group_vars, markers, new_name=None, res=1, method='t-test_overestim_var', plot_maps=True, **kwargs):
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
        - new_name: new name for sub-clustered labels
        - res: resolution of louvain clustering
        - method: DE method for marker assignment
        - kwargs:
            - tweaks: dictionary of specific edits to make to labeling {'tcell,1': 'tcell_cd4', 'tcell,2':'tcell_cd8'}
            (useful if assignment is incorrect for rare populations)


    Outputs:
        - genes: marker genes used to label clusters
    """
    if new_name is None:
        new_name = group+'_R'

    # sub-cluster
    sc.tl.louvain(adata, resolution=res, restrict_to = (group, group_vars), key_added=new_name)

    # plot UMAP
    if plot_maps: sc.pl.umap(adata, color=[new_name])

    # Build marker set and reassign
    scores, genes, labels = build_marker_set(adata, markers, groupby=new_name, key_added=new_name, **kwargs)

    # Plot new maps
    if plot_maps: sc.pl.umap(adata, color=[new_name])

    return genes
