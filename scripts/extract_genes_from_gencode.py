#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract genes info from gencode table
"""

# %%
import pandas as pd
import gzip

# %%
gencode_gene_list = []

with gzip.open(snakemake.input.gencode, 'r') as f:

    # Header
    for i in range(7):
        f.readline()

    # Data
    for line in f.readlines():
        line = line.decode('UTF-8')

        line_list = line.split('\t')

        if len(line_list) > 2 and line_list[2] == 'gene':
            chr_id, start, end = line_list[0], line_list[3], line_list[4]
            info = line_list[8].split(';')
            
            for i in info:
                if i.startswith('gene_name'):
                    gene_name = i.split('=')[1]

                if i.startswith('gene_type'):
                    gene_type = i.split('=')[1]

            gencode_gene_list.append([gene_name, chr_id, start, end, gene_type])

gencode_genes = pd.DataFrame(gencode_gene_list,
                             columns=['gene_name', 'chr', 'start', 'end', 'gene_type'])

# %%
gencode_genes.to_csv(snakemake.output.gencode_genes, sep='\t', index=False)
