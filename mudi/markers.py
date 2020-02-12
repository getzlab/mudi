import gc
import pandas as pd
from tqdm import tqdm
import os
import scanpy as sc
import scanpy.external as sce
from anndata import AnnData
import sys
from collections import defaultdict
import numpy as np
from typing import Union

from .utils import aggr_markers

# ---------------------------------
# Cell-type Calling
# ---------------------------------
def labeler(d,x):
    try:
        return d[x]
    except:
        return 'n/a'

def convert_marker_dict_to_df(marker_dict: dict):
    """
    Converts marker dictionary format to a marker dataframe.
    --------------------------------------
    Example marker_dict:

    {'tcell_activated': ['CD69', 'IL2RA'],
     'tcell_effector': ['CD3D', 'B3GAT1', 'PDCD1', 'FAS', 'CCR7'],
     'tcell_regulatory': ['CD4',
      'IL2RA',
      'FOXP3',
      'SMAD3',
      'STAT5A',
      'STAT5B',
      'IL10'],
     'tcell_exhausted': ['CD3D',
      'PDCD1',
      'FAS',
      'CCR7', .......

    Example output:
    --------------------------------------
    Cell-Type	Gene
    0	tcell_activated	CD69
    1	tcell_activated	IL2RA
    2	tcell_effector	CD3D
    3	tcell_effector	B3GAT1
    4	tcell_effector	PDCD1
    5	tcell_effector	FAS
    6	tcell_effector	CCR7
    7	tcell_regulatory	CD4
    8	tcell_regulatory	IL2RA

    """
    pairs = list()
    for key in marker_dict:
        for gene in marker_dict[key]:
            pairs.append((key,gene))
    return pd.DataFrame(pairs).rename(columns={0:'Cell-Type',1:'Gene'})

def build_marker_set(
    adata: AnnData, \
    markers_df: pd.DataFrame, \
    groupby: str = 'louvain', \
    metric: str = 'sum', \
    pval_thresh: float = 1e-5, \
    lfc_thresh: float = 1.0, \
    key_added: str = 'cell_type', \
    cell_type_idx: str = 'Cell-Type', \
    gene_idx: str = 'Gene', \
    tweaks: Union[None, dict] = None, \
    **kwargs \
    ):
    """
    Build Marker set.
    ---------------------
    This takes in an Anndata object, a dictionary mapping from
    cell-types to markers of interest.

    Inputs:
        - adata: scanpy anndata object
        - markers_df: markers dataframe or dictionary
        - metric: metric for assigning cell-type
        - pval_thresh: threshold for markers (DE adj. pval)
        - lfc_thresh: threshold for markers (DE LFC)
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
    markers = markers[(markers['pvals_adj']<pval_thresh) & (markers['logfoldchanges']>lfc_thresh)]

    if isinstance(markers_df, dict):
        markers_df = convert_marker_dict_to_df(markers_df)

    d = {}
    dfs = list()
    for cell_type in np.unique(markers_df[cell_type_idx]):
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

def sub_cluster_and_rename(
    adata: AnnData, \
    group: str, \
    group_vars: list, \
    markers: pd.DataFrame, \
    new_name: Union[None, str] = None, \
    res: float = 1.0, \
    method: str = 't-test_overestim_var', \
    plot_maps: bool = True, \
    **kwargs
    ):
    """
    Subcluster & Rename
    ---------------------
    Takes in an adata object and information about grouping/variables, does subclustering, and
    reassigns variables.

    Inputs:
        - adata: scanpy anndata object
        - group: grouping to do sub-clustering on (ex. louvain)
        - group_vars: sub-groups within grouping to do sub-clustering on (ex. ['1','2'] for louvain clusters 1 & 2)
        - markers: pd.DataFrame of markers
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

def compile_df_de(adata: pd.DataFrame, markers: pd.DataFrame):
    """
    Compiles mapped markers to a de-gene dictionary and adds
    this information into the resulting dataframe.
    ---------------------
    Args:
        * adata: AnnData object with rank_genes_groups computed for a grouping
        * markers: markers dataframe or dictionary for inputs

    Returns:
        * dataframe with marker genes including a "subtypes" column with all
            overlapping subtypes
    """
    if isinstance(markers, dict):
        markers = convert_marker_dict_to_df(markers)

    df_de = aggr_markers(adata)
    df_de['cluster'] = df_de.cluster.astype('category')

    # Reverse Marker Mapping
    # Map genes --> [putative_cell_types...,]
    r_marker_mapping = defaultdict(list)
    for idx,row in markers.iterrows():
        r_marker_mapping[row['Gene']].append(row['Cell-Type'])

    df_de['subtypes'] = df_de['names'].apply(lambda x: ','.join(r_marker_mapping[x]))
    df_de['cluster'] = df_de['cluster'].apply(lambda x: x.replace('/',''))
    df_de['cluster'] = df_de['cluster'].apply(lambda x: x.replace('?',''))
    return df_de

# ---------------------------------
# Write to Excel file
# ---------------------------------
def to_xlsx(filename, df):
    """
    Write markers dataframe file to excel file.
    ---------------------
    Args:
        * filename: name of output excel file
        * df: marker dataframe
            ** computed either by aggr_markers or compile_de_df
    """
    if not filename.endswith('.xlsx'):
        filename = filename + '.xlsx'

    with pd.ExcelWriter(filename) as writer:
        for c in sorted(df.cluster.cat.categories, key=str):
            df[df.cluster == c].set_index('no.').drop(columns='cluster').to_excel(writer, sheet_name=f'{c}')
