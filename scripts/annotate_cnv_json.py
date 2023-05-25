#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# %%
import pandas as pd
import json
import numpy as np
import gzip


class NumpyEncoder(json.JSONEncoder):
    """ Custom encoder for numpy data types """
    
    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                            np.int16, np.int32, np.int64, np.uint8,
                            np.uint16, np.uint32, np.uint64)):

            return int(obj)

        elif isinstance(obj, (np.float_, np.float16, np.float32, np.float64)):
            return float(obj)

        elif isinstance(obj, (np.complex_, np.complex64, np.complex128)):
            return {'real': obj.real, 'imag': obj.imag}

        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()

        elif isinstance(obj, (np.bool_)):
            return bool(obj)

        elif isinstance(obj, (np.void)):
            return None

        return json.JSONEncoder.default(self, obj)
# %%
# cnv_chrom, cnv_start, cnv_end, cnv_type = 18, 79805149, 79942630, 'gain'
cnv_chrom = snakemake.wildcards.chrom
cnv_start = int(snakemake.wildcards.start)
cnv_end = int(snakemake.wildcards.end)
cnv_type = snakemake.wildcards.cnv_type

minimal_overlap = 1

# genes = pd.read_csv("../data/generated_data/genes_preprocessed.tsv", sep='\t', index_col=0)
# genes_scaled = pd.read_csv("../data/generated_data/genes_preprocessed_scaled.tsv", sep='\t', index_col=0)
# genes_unimputed = pd.read_csv("../data/generated_data/genes_preprocessed_unimputed.tsv", sep='\t', index_col=0)

genes = pd.read_csv(snakemake.input.genes_preprocessed, sep='\t', index_col=0)
genes_scaled = pd.read_csv(snakemake.input.genes_preprocessed_scaled, sep='\t', index_col=0)
genes_unimputed = pd.read_csv(snakemake.input.genes_preprocessed_unimputed, sep='\t', index_col=0)
genes_raw = pd.read_csv(snakemake.input.genes_raw, sep='\t', index_col=0)

# %%

candidates = genes.query(f"chr == 'chr{cnv_chrom}'")

# %%
# Overlapping genes

# from left side
left = candidates.query(f"start < {cnv_start} & end > {cnv_start + minimal_overlap}")

# inside
inside = candidates.query(f"start >= {cnv_start} & end <= {cnv_end}")

# from right side
right = candidates.query(f"start < {cnv_end - minimal_overlap} & end > {cnv_end}")

genes_in_cnv = pd.concat([left, inside, right])    
genes_in_cnv = genes_in_cnv.drop_duplicates()

genes_in_cnv_scaled = genes_scaled.iloc[genes_in_cnv.index]
genes_in_cnv_unimputed = genes_unimputed.iloc[genes_in_cnv.index]
genes_in_cnv_raw = genes_raw.iloc[genes_in_cnv.index]

# %% Make json file
genes_json = []

for i in range(genes_in_cnv.shape[0]):
    g = {}
    for j in range(4, genes_in_cnv.shape[1]):
        attribute = genes_in_cnv.columns[j]

        try:
            raw_value = genes_in_cnv_raw.iloc[i].loc[attribute]
        except KeyError:
            # the attribute is dummied (value is appended to attribute name
            # as {attribute}_{value})
            att = '_'.join(attribute.split('_')[:-1])
            raw_value = genes_in_cnv_raw.iloc[i].loc[att]

        g[attribute] = {
            "imputed_unscaled": genes_in_cnv.iloc[i].loc[attribute],
            "imputed_scaled": genes_in_cnv_scaled.iloc[i].loc[attribute],
            "unimputed_unscaled": genes_in_cnv_unimputed.iloc[i].loc[attribute],
            "raw": raw_value
        }
    genes_json.append(g)

result = {
    "chrom": cnv_chrom,
    "start": cnv_start,
    "end": cnv_end,
    "cnv_type": cnv_type,
    "genes": genes_json
}
    
# %%
with gzip.open(snakemake.output.cnv, 'wt') as f:
    json.dump(result, f, indent=2, cls=NumpyEncoder)
