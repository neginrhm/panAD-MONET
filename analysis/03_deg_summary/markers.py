"""
markers.py
==========
Define canonical and subtype-specific marker gene sets for major brain cell
types in the 13-AD dataset, and generate a dot plot visualization.

Inputs
------
- A preprocessed AnnData/Scanpy object (variable `all1`) with Leiden clustering.

Outputs
-------
- Dot plot of marker expression across clusters (displayed interactively).

Usage
-----
    # Run within an interactive session (e.g., Jupyter) after loading `all1`:
    exec(open("markers.py").read())
"""

import scanpy as sc

# --- Marker gene definitions -------------------------------------------------

marker_dict = {

    'Vasculature': ['PDGFRB', 'FLT1'],
    'Excitatory Neuron': ['RORB', 'THEMIS'],
    'Inhibitory Neuron': ['GAD2', 'SST'],

    'Endothelial cells': ['PECAM1', 'VWF', 'VEGFA'],

    'Microglia': ['P2RY12', 'CX3CR1', 'CD4', 'CD68'],
    'Micro-PVM_1': ['MRC1', 'RGL1', 'AGFG1'],
    'Monocyte': ['FCN1', 'CD300E', 'LILRA5'],
    'Micro-PVM_2': ['FAM117B', 'DLEU7', 'ITM2B'],
    'Micro-PVM_2_1': ['SPC25', 'CIT', 'NSD2'],
    'Lymphocyte': ['SKAP1', 'SLC38A1', 'PARP8'],
    'Micro-PVM_2_3': ['CSMD1', 'KCNIP4', 'SH3RF3'],
    'Micro-PVM_3': [
        'LINC02698', 'PTPRG', 'TRHDE', 'AC131956.2', 'CADM1', 'TPRG1',
        'WIPF3', 'PARVG', 'RALGAPA2', 'KCNJ5', 'MGLL', 'FPR3', 'ZNF804A',
        'SLC47A1', 'AL138756.1', 'DSCAM', 'SAMD4A', 'LIPA', 'NHSL1',
        'RASGRP3', 'TBC1D16', 'TANC2', 'ATG7', 'DOP1B', 'SNX29', 'GLDN',
        'FRMD4B', 'RNF130', 'SLC1A3', 'PLXDC2'
    ],
    'Micro-PVM_4': ['RAPGEF1', 'PDE8A', 'MAP4K3'],

    'Astrocyte': ['GFAP', 'SLC1A2'],
    'Astro_1': ['NIN', 'MGAT4C', 'IQCA1'],
    'Astro_2': [
        'WIF1', 'GRM3', 'SGCD', 'RORA', 'DNAH7', 'SH3RF1', 'CADM2',
        'HPSE2', 'VAV3', 'PDE3B', 'SYNE1', 'ADGRB3', 'GPM6A', 'FLRT2',
        'LRP1B', 'MERTK', 'COL5A3', 'ARHGAP24', 'SLC1A2', 'ZNF98',
        'CABLES1', 'LSAMP', 'SOX5', 'CARMIL1', 'NRXN1', 'ADGRL3', 'PLCB1',
        'GPC5', 'DNM3', 'CADM1'
    ],
    'Astro_3': ['NFASC', 'LINC00609', 'KAZN'],
    'Astro_4': ['RASGRF1', 'LDB2', 'SH3GL2'],
    'Astro_5': ['ADAMTS18', 'SHISA6', 'ADGRV1'],
    'Astro_6': ['CHI3L1', 'F3', 'NRCAM'],

    'Oligodendrocyte': ['MOBP', 'MOG', 'MBP', 'PLP1'],
    'Oligo_1': ['AFF3', 'RASGRF1', 'RBFOX1'],
    'Oligo_2': [
        'CLASP2', 'TJP1', 'CDC42BPA', 'PIP4K2A', 'SHTN1', 'PDE1C',
        'SPOCK1', 'ZBTB20', 'ZNF536', 'DIP2B', 'ANK2', 'LIMCH1',
        'IL1RAPL1', 'KLHL32', 'FRMD5', 'DYSF', 'PTPRD', 'DIP2C',
        'SLC44A1', 'CTNNA3', 'ARAP2', 'PDE4B', 'COBL', 'FNBP1', 'FRMD4B',
        'LINC01608', 'PLEKHH1', 'PLXDC2', 'AUTS2', 'PALM2-AKAP2'
    ],
    'Oligo_3': ['CNTN1', 'GRIK2', 'KCNIP4'],
    'Oligo_4': [
        'RPL18', 'RPL10', 'RPS28', 'CBR1', 'PTGDS', 'RPLP1', 'LRRTM3',
        'ENO4', 'HEPN1', 'GCNT1', 'RPS27', 'B2M', 'EIF1', 'PCA3', 'FTL',
        'CRYAB', 'FTH1', 'H3F3B', 'PDE5A', 'SELENOP', 'RNASE1', 'DBNDD2',
        'ANKUB1', 'TUBA1A', 'MALAT1', 'PTMA', 'CD81', 'SPP1', 'CNP',
        'PLCG2'
    ],

    'OPC': ['PDGFRA', 'CSPG4'],
    'OPC_1': ['FYN', 'BCAS1', 'FRMD4A'],
    'OPC_2': [
        'S100B', 'SUMF1', 'LINC01572', 'LINC00402', 'MGST1', 'FMO3',
        'OMG', 'TSBP1', 'AF064866.1', 'RGCC', 'IL7', 'MROH9', 'MFSD4B',
        'OLIG1', 'TTC6', 'GLIPR1L1', 'PLCG2', 'HEPN1', 'LINC00486',
        'BCL2L14', 'MAK', 'BCHE', 'P2RY14', 'ITM2C', 'PDE5A', 'LINC00960',
        'CPB2', 'SHOC1', 'PRSS48', 'SPARCL1'
    ],
    'OPC_2_2': ['ARHGEF3', 'CAMK2D', 'SORBS1']
}

# --- Plot --------------------------------------------------------------------

sc.set_figure_params(dpi=100)
sc.pl.dotplot(
    all1, marker_dict, 'leiden',
    dendrogram=True, standard_scale='var', swap_axes=True, cmap='RdPu'
)
