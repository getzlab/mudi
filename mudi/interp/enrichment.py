from gprofiler import gprofiler
import pandas as pd

from ..markers import aggr_markers

def gprof(adata,
          key_added='gsea',
          custom_bg='raw',
          ordered_query=True,
          src_filter=None,
          **kwargs):
    """
    G-profiler wrapper for gene set enrichment analysis.
    See here: https://biit.cs.ut.ee/gprofiler/page/apis

    Gprofiler query: list of genes/proteins; can be mixed types of gene IDs, SNP IDs,
        chromosomal intervals, or term IDs

    Requires differentially expressed genes to be computed first via
        sc.tl.rank_genes_groups(adata, groupby=<CLUSTER>)

    Inputs:
        * key_added: key to add in adata.uns
            if None   --> then does not modify AnnData input object, just returns the
                results
            if (str)  --> adds key to adata.uns with enrichment data
        * custom_bg (list): custom background set;
            if 'raw'  --> adata.raw.varnames <default>
            if None   --> use default background selection for gProfiler
            if (list) --> use provided genes as background
        * ordered_query (bool): g:Profiler gene lists may be interpreted as ordered
            lists where elements (genes, proteins, probesets) are in the order
            of decreasing importance. The ordered query option is useful when
            the genes can be placed in some biologically meaningful order.
            For instance, according to differential expression in a given
            RNA-seq experiment. g:Profiler then performs incremental enrichment analysis
            with increasingly larger numbers of genes starting from the top of
            the list. This procedure identifies specific functional terms that
            associate to most dramatic changes in given experimental setup, as
            well as broader terms that characterise the gene set as a whole;
            True <default>
        * src_filter (list): list of the following database sources:
            Source ID:
                GO:MF	molecular function
                GO:CC	cellular component
                GO:BP	biological process
                KEGG	Kyoto Encyclopedia of Genes and Genomes
                REAC	Reactome
                WP	    WikiPathways
                TF	    Transfac
                MIRNA	miRTarBase
                HPA	    Human Protein Atlas
                CORUM	CORUM protein complexes
                HP	    Human Phenotype Ontology

    kwargs:
        * organism (str): 'hsapiens' <default>
        * signifcant (bool): True <default>
        * exclude_iea (bool): include electronic GO annotations: annotations assigned to
            genes using in silico curation methods adn will ahve an IEA evidence
            code to back this up; ecluding IEA annotations may help reduce bias towards
            abundant and ubiquitous housekeeping genes like ribosomal genes;
            False <default>
        * max_p_value (float): 1.0 <default>
        * max_set_size (int): 0 <default>
        * correction_method (str): 'analytical': gcustom g:SCS for reducing signifance -
            traditional multiple tests corrections are for tests independent of one another
            which is not correct for GO analysis; 'fdr': benjamini-hochberg, 'bc': bonferroni correction; 'analytical' <default>


    Outputs:
        if key_added is provided:
            Adds the enrichments for each rank-gene cluster to adata.uns. This may
            be retrieved by calling "get_uns(adata, <KEY_ADDED>)".
        if key_added is None:
            Returns a dataframe for enrichment results for each possible cluster

    """
    assert adata.uns['rank_genes_groups'], \
        "Compute differentially expressed genes first."

    if custom_bg is 'raw':
        custom_bg = list(adata.raw.var_names)
    elif custom_bg is None:
        custom_bg = list()

    de_markers = aggr_markers(adata)
    enrichments = list()

    for clust in set(de_markers['cluster']):
        df = de_markers[de_markers['cluster']==clust]

        enrichment = gprofiler(
            df.names,
            custom_bg=custom_bg,
            ordered_query=ordered_query,
            src_filter=src_filter,
            **kwargs
        )

        enrichment['cluster'] = clust
        enrichment = enrichment.sort_values('p.value').loc[:,['cluster','p.value','term.size','overlap.size','term.name','domain','intersection']]
        enrichments.append(enrichment)

    enrichments = pd.concat(enrichments)

    if key_added is None:
        return enrichments
    else:
        _enrichment = {}
        _enrichment['cols'] = list(enrichments)
        _enrichment['rows'] = enrichments.index
        _enrichment['values'] = enrichments.values
        _enrichment['src'] = src_filter
        adata.uns[key_added] = _enrichment



















    #
