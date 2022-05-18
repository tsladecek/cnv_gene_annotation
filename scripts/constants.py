#!/usr/bin/env python3
# -*- coding: utf-8 -*-

IMPUTE_VALUE = "mean"
SCALING = "standardize"
SD_DEFAULT = 0
CNV_SPECIFIC_ATT_NUM = 23


DBNSFP_GENE_ATTRIBUTES_MODELLING = [
    # 'Gene_name', 'Ensembl_gene', 'chr', 'Gene_old_names',
    #    'Gene_other_names', 'Uniprot_acc(HGNC/Uniprot)',
    #    'Uniprot_id(HGNC/Uniprot)', 'Entrez_gene_id', 'CCDS_id', 'Refseq_id',
    #    'ucsc_id', 'MIM_id', 'OMIM_id', 'Gene_full_name', 'Pathway(Uniprot)',
    #    'Pathway(BioCarta)_short', 'Pathway(BioCarta)_full',
    
    # Binarize
    # -------------------------------------------------------------
    'Pathway(ConsensusPathDB)',
    #  'Pathway(KEGG)_id', 'Pathway(KEGG)_full',
       'Function_description', 'Disease_description', 
    #  'MIM_phenotype_id',
       'MIM_disease', 
    #  'Orphanet_disorder_id',
       'Orphanet_disorder',
    #  'Orphanet_association_type', 
       'Trait_association(GWAS)',
    # 'GO_biological_process', 'GO_cellular_component',
    # 'GO_molecular_function',
    # 'Tissue_specificity(Uniprot)',
       'Expression(egenetics)', 
    # 'Expression(GNF/Atlas)',
    #  'Interactions(IntAct)', 'Interactions(BioGRID)',
    # Count
    # -------------------------------------------------------------
       'Interactions(ConsensusPathDB)',

    # NUMERICAL
    # -------------------------------------------------------------
       'P(HI)', 'HIPred_score',
    # CATEGORICAL
    # -------------------------------------------------------------
       'HIPred',  # 2 classes = ['Y', 'N']
    # -------------------------------------------------------------
       'GHIS', 'P(rec)',

    # CATEGORICAL
    # -------------------------------------------------------------
       'Known_rec_info',  # 2 classes = ['lof-tolerant', 'recessive']
    # -------------------------------------------------------------
    # NUMERICAL
    # -------------------------------------------------------------
       'RVIS_EVS', 'RVIS_percentile_EVS',
       'LoF-FDR_ExAC', 'RVIS_ExAC', 'RVIS_percentile_ExAC', 'ExAC_pLI',
       'ExAC_pRec', 'ExAC_pNull', 'ExAC_nonTCGA_pLI', 'ExAC_nonTCGA_pRec',
       'ExAC_nonTCGA_pNull', 'ExAC_nonpsych_pLI', 'ExAC_nonpsych_pRec',
       'ExAC_nonpsych_pNull', 'gnomAD_pLI', 'gnomAD_pRec', 'gnomAD_pNull',
       'ExAC_del.score', 'ExAC_dup.score', 'ExAC_cnv.score',

    # CATEGORICAL
    # -------------------------------------------------------------
       'ExAC_cnv_flag',  # 2 classes = ['Y', 'N']
    # -------------------------------------------------------------
    # NUMERICAL
    # -------------------------------------------------------------
       'GDI', 'GDI-Phred',
    # CATEGORICAL
    # -------------------------------------------------------------
    # 3 classes = ['low', 'medium', 'high']
       'Gene damage prediction (all disease-causing genes)',
       'Gene damage prediction (all Mendelian disease-causing genes)',
       'Gene damage prediction (Mendelian AD disease-causing genes)',
       'Gene damage prediction (Mendelian AR disease-causing genes)',
       'Gene damage prediction (all PID disease-causing genes)',
       'Gene damage prediction (PID AD disease-causing genes)',
       'Gene damage prediction (PID AR disease-causing genes)',
       'Gene damage prediction (all cancer disease-causing genes)',
       'Gene damage prediction (cancer recessive disease-causing genes)',
       'Gene damage prediction (cancer dominant disease-causing genes)',
    # -------------------------------------------------------------
    # NUMERICAL
    # -------------------------------------------------------------
       'LoFtool_score', 'SORVA_LOF_MAF0.005_HetOrHom',
       'SORVA_LOF_MAF0.005_HomOrCompoundHet', 'SORVA_LOF_MAF0.001_HetOrHom',
       'SORVA_LOF_MAF0.001_HomOrCompoundHet',
       'SORVA_LOForMissense_MAF0.005_HetOrHom',
       'SORVA_LOForMissense_MAF0.005_HomOrCompoundHet',
       'SORVA_LOForMissense_MAF0.001_HetOrHom',
       'SORVA_LOForMissense_MAF0.001_HomOrCompoundHet',

    # CATEGORICAL
    # -------------------------------------------------------------
       'Essential_gene',
       'Essential_gene_CRISPR', 'Essential_gene_CRISPR2',
       'Essential_gene_gene-trap', 
    # NUMERICAL
    # -------------------------------------------------------------       
       'Gene_indispensability_score',
    # CATEGORICAL
    # -------------------------------------------------------------      
       'Gene_indispensability_pred',
       # 'MGI_mouse_gene', 'MGI_mouse_phenotype',
       # 'ZFIN_zebrafish_gene', 'ZFIN_zebrafish_structure',
       # 'ZFIN_zebrafish_phenotype_quality', 'ZFIN_zebrafish_phenotype_tag'
    ]

DBNSFP_NUMERICAL = [
    'P(HI)', 'HIPred_score',
    'GHIS', 'P(rec)',
    'RVIS_EVS', 'RVIS_percentile_EVS',
    'LoF-FDR_ExAC', 'RVIS_ExAC', 'RVIS_percentile_ExAC', 'ExAC_pLI',
    'ExAC_pRec', 'ExAC_pNull', 'ExAC_nonTCGA_pLI', 'ExAC_nonTCGA_pRec',
    'ExAC_nonTCGA_pNull', 'ExAC_nonpsych_pLI', 'ExAC_nonpsych_pRec',
    'ExAC_nonpsych_pNull', 'gnomAD_pLI', 'gnomAD_pRec', 'gnomAD_pNull',
    'ExAC_del.score', 'ExAC_dup.score', 'ExAC_cnv.score',
    'GDI', 'GDI-Phred',
    'LoFtool_score', 'SORVA_LOF_MAF0.005_HetOrHom',
    'SORVA_LOF_MAF0.005_HomOrCompoundHet', 'SORVA_LOF_MAF0.001_HetOrHom',
    'SORVA_LOF_MAF0.001_HomOrCompoundHet',
    'SORVA_LOForMissense_MAF0.005_HetOrHom',
    'SORVA_LOForMissense_MAF0.005_HomOrCompoundHet',
    'SORVA_LOForMissense_MAF0.001_HetOrHom',
    'SORVA_LOForMissense_MAF0.001_HomOrCompoundHet',
    'Gene_indispensability_score'
    ]

BINARIZE = [
         'Pathway(ConsensusPathDB)',
         'Function_description', 
         'Disease_description',
         'MIM_disease',
         'Orphanet_disorder',
         'Trait_association(GWAS)',
         'Expression(egenetics)', 
        
    ]

COUNT = ['Interactions(ConsensusPathDB)']

DUMMY = [
        'HIPred',
        'Known_rec_info',
        'ExAC_cnv_flag',
        'Gene damage prediction (all disease-causing genes)',
        'Gene damage prediction (all Mendelian disease-causing genes)',
        'Gene damage prediction (Mendelian AD disease-causing genes)',
        'Gene damage prediction (Mendelian AR disease-causing genes)',
        'Gene damage prediction (all PID disease-causing genes)',
        'Gene damage prediction (PID AD disease-causing genes)',
        'Gene damage prediction (PID AR disease-causing genes)',
        'Gene damage prediction (all cancer disease-causing genes)',
        'Gene damage prediction (cancer recessive disease-causing genes)',
        'Gene damage prediction (cancer dominant disease-causing genes)',
        'Essential_gene',
        'Essential_gene_CRISPR',
        'Essential_gene_CRISPR2',
        'Essential_gene_gene-trap', 
        'Gene_indispensability_pred'
    ]
