# IMMUNE_MARKERS : less granularity
# IMMUNE_BROAD_MARKERS : more granularity

# ------------------------------------------
# Markers
# ------------------------------------------
bcell = {'bcell_activated': ['CD19','IL2RA','CD30'],
         'bcell_plasma': ['CD27','CD38','SDC1','SLAMF7','IL6'],
         'bcell_memory': ['MS4A1','CD27','CD40','CD80','PDCD1LG2', 'CXCR3','CXCR4','CXCR5','CXCR6'],
         'bcell_marginal': ['CD1A','CR2','CD37','NOTCH2'],
         'bcell_follicular': ['CR2','CD22','FCER2'],
         'bcell_regulatory': ['CD1A','CD5','CR2','CD24','TLR4','IL10','TGFB1']}

tcell = {#'tcell_general': ['CD3D','CD3E','CD3G','PTPRC'],
         'tcell_activated': ['CD69','IL2RA'],
         'tcell_effector': ['CD3D','B3GAT1','PDCD1','FAS','CCR7'],
         'tcell_regulatory': ['CD4','IL2RA','FOXP3','SMAD3','STAT5A','STAT5B','IL10'],
         'tcell_exhausted': ['CD3D','PDCD1','FAS','CCR7','SELL','LAG3','HAVCR2','TIGIT','ENTPD1'], #CD38?
         'tcell_helper_1': ['CD4','CXCR3','TBX21','STAT1','STAT6','IFNG'],
         'tcell_helper_2': ['CD4','CCR4','PTGDR2','GATA3','STAT5','STAT6','IL4','IL5','IL13'],
         'tcell_naive': ['IL7R','S100A4','CD3D','SELL','CD27'],
         'tcell_cytotoxic': ['CD8A', 'CD8B'],
         'tcell_memory': ['IL7R','CD3D','CCR7','SELL'],
         'tcell_helper_17': ['CD4','CCR6','RORC','RORA','STAT3','IL17A','IL17F']}

monocyte = {'monocyte_inflammatory': ['CCR2'],
             'monocyte_resident': ['CXCR1'],
             'monocyte_CD14+': ['CD14', 'FCN1', 'LYZ'],
             'monocyte_FCGR3A+': ['CD14','FCGR3A','MS4A7'],
             'macrophages': ['CD68','CCR5','TFRC','ITGAM','FCGR1A','CSF1R','MRC1','CD163']}

nk = {'nk_cell': ['NKG7','GNLY','NCR1','KLRD1','NCAM1']}

dc = {'DC': ['FCER1A', 'ITGAX', 'CD83', 'THBD','CD209','CD1C', 'LYZ']}

pdc =  {'pDC': ['IL3RA','CLEC4C','NRP1']}

# PPBP, CD41b, CD42a, CD42b, CD61
megakaryocyte = {'megakaryocyte': ['PPBP', 'ITGA2B','GP9', 'GP1BA', 'ITGB3']}

# CD235a - GYPA
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5957548/
erythrocyte = {'erythrocyte': ['GYPA', 'BLVRB', 'HBB', 'HBA1']}

neutrophil = {'neutrophil': ['FUT4','ITGAM','FCGR3A','ITGB2','FCGR2A','CD44','CD55']}

eosinophil = {'eosionophil': ['IL5RA','CCR3','EMR1','ADGRE1','SIGLEC8']}

basophil = {'basophil': ['IL3RA', 'ENPP3']}

mast = {'mast': ['KIT','CD33','TPSAB1', 'CD9']}

hsc = {'hsc': ['CD34','THY1', 'CDK6']}

# Full Markers
markers = {'tcell': tcell,
           'bcell': bcell,
           'monocytes': monocyte,
           'nk':nk,
           'DC': dc,
           'pDC': pdc,
           'megakaryocyte': megakaryocyte,
           'erythrocyte': erythrocyte,
           'neutrophil': neutrophil,
           'eosinophil': eosinophil,
           'basophil': basophil,
           'mast': mast,
           'hsc': hsc
          }

# Specific Markers
IMMUNE_MARKERS = {
    **tcell,
    **bcell,
    **monocyte,
    **nk,
    **dc,
    **pdc,
    **megakaryocyte,
    **erythrocyte,
    **neutrophil,
    **eosinophil,
    **basophil,
    **mast,
    **hsc
}

# General Markers
IMMUNE_BROAD_MARKERS = {key:{item for row in markers[key].values() for item in row} for key in markers}
