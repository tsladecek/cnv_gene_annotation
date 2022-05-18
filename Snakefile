include: "scripts/input_function.py"


rule all:
    input:
        all_cnvs


rule download_gencode_data:
    output:
        "data/gencode.v35.annotation.gff3.gz"
    run:
        shell("wget -O data/gencode.v35.annotation.gff3.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gff3.gz")

rule download_dbnsfp_data:
    output:
        "data/dbNSFP4.0_gene.complete.gz"
    shell:
        "wget -O data/dbNSFP4.0_gene.complete.gz https://s3.eu-central-1.amazonaws.com/jbrowse.genovisio.com/data_hg38/tracks/genes/dbNSFP4.0_gene.complete.gz"

rule gencode_genes:
    input:
        gencode = "data/gencode.v35.annotation.gff3.gz"
    output:
        gencode_genes = "data/generated_data/gencode_genes.tsv"
    script:
        "scripts/extract_genes_from_gencode.py"

rule annotate_gencode_with_dbnsfp:
    input:
        dbnsfp = "data/dbNSFP4.0_gene.complete.gz",
        gencode_genes = "data/generated_data/gencode_genes.tsv"
    output:        
        protein_coding = "data/generated_data/gencode_dbnsfp_protein_coding.tsv"
    script:
        "scripts/gencode_dbnsfp.py"

rule preprocess_genes_data:
    input:
        protein_coding = "data/generated_data/gencode_dbnsfp_protein_coding.tsv"
    output:
        genes_preprocessed = "data/generated_data/genes_preprocessed.tsv",
        genes_preprocessed_scaled = "data/generated_data/genes_preprocessed_scaled.tsv",
        genes_preprocessed_unimputed = "data/generated_data/genes_preprocessed_unimputed.tsv"
    script:
        "scripts/genes_preprocessing.py"


rule annotate_cnv:
    input:
        genes_preprocessed = "data/generated_data/genes_preprocessed.tsv",
        genes_preprocessed_scaled = "data/generated_data/genes_preprocessed_scaled.tsv",
        genes_preprocessed_unimputed = "data/generated_data/genes_preprocessed_unimputed.tsv",
    output:
        cnv = "data/annotated_cnvs_json/{chrom}_{start}_{end}_{cnv_type}.json.gz"
    script:
        "scripts/annotate_cnv_json.py"
