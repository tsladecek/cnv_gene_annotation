#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 14:42:32 2020

@author: tomas
"""

# %% imports
import pandas as pd
import numpy as np

# %%
dbnsfp = pd.read_csv(snakemake.input.dbnsfp, compression='gzip', sep='\t', low_memory=False)
dbnsfp.replace('.', np.nan, inplace=True)

gencode_genes = pd.read_csv(snakemake.input.gencode_genes, sep='\t')
gencode_genes.rename(columns = {'gene_name':'Gene_name'}, inplace=True)
# %% add "chr" before chromosome name
dbnsfp.chr = [f'chr{i}' for i in dbnsfp.chr]

# %% Merge dataframes on Gene name 

gencode_dbnsfp = gencode_genes.merge(dbnsfp, how='inner', on=['Gene_name', 'chr'])

# %% Only protein coding genes
protein_coding = gencode_dbnsfp.query("gene_type == 'protein_coding'")

protein_coding = protein_coding.reset_index(drop=True)

#protein_coding.to_csv('../data/generated_data/gencode_dbnsfp_protein_coding.tsv', sep='\t')
protein_coding.to_csv(snakemake.output[0], sep='\t', index=False)
