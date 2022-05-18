# GENE ANNOTATION

Get genes inside of a CNV and annotate them with dbNSFP4

1. Create Conda Environment

```
conda env create --file environment.yml 
conda activate gene_annot
```

2. Run pipeline

```
snakemake --cores <number of cores> 
```
