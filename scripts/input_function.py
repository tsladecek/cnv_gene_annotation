import os


cnv_types = { 'dup': 'gain', 'gain': 'gain', 'del': 'loss', 'loss': 'loss' }

def all_cnvs(wildcards):
    """This should be populated with all CNVs desired to be annotated
    
    see example below with list of 2 CNVs
    """
    input_filename = 'input_example.bed' if 'input.bed' not in os.listdir() else 'input.bed'
    with open(input_filename) as f:
        
        outputs = []
        for cnv in f.readlines():
            chrom, start, end, cnv_type = cnv.strip().split('\t')
            chrom = chrom.strip('chr')
            try:
                cnv_type = cnv_types[cnv_type.lower()]
            except KeyError:
                exit(1)
            else:
                outputs.append(os.path.join('data', 'annotated_cnvs_json', f'{chrom}_{start}_{end}_{cnv_type}.json.gz'))

    return outputs

