#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Preprocessing of the genes dataset
"""

# %% imports
import pandas as pd
from constants import DBNSFP_GENE_ATTRIBUTES_MODELLING
from data_preprocessing import preprocess

# %% load datasets
genes_raw = pd.read_csv(snakemake.input.protein_coding, sep='\t')# .drop("Unnamed: 0", axis=1)
# genes_raw = pd.read_csv('../data/generated_data/gencode_dbnsfp_protein_coding.tsv', sep='\t')# .drop("Unnamed: 0", axis=1)
# %% DATA PREPROCESSING
columns = ['Gene_name', 'chr', 'start', 'end']
columns.extend(DBNSFP_GENE_ATTRIBUTES_MODELLING)

genes0 = genes_raw.loc[:, columns]

# %%
# replace all missing numerical values with mean values
genes = preprocess(genes0,
                   impute_numerical=True,
                   scale_numerical=False)

genes.to_csv(snakemake.output.genes_preprocessed, sep='\t')

# %% scale
genes = preprocess(genes0,
                   impute_numerical=True,
                   scale_numerical=True)

genes.to_csv(snakemake.output.genes_preprocessed_scaled, sep='\t')

# Unimputed unscaled
genes = preprocess(genes0,
                   impute_numerical=False,
                   scale_numerical=False)

genes.to_csv(snakemake.output.genes_preprocessed_unimputed, sep='\t')
