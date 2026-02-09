import pandas as pd
import os

# Load configuration
configfile: "config.yml"

# Load samples
samples_df = pd.read_csv(config["input_tsv"], sep="\t").set_index("sample", drop=False)

rule all:
    input:
        expand(os.path.join(config["output_base"], "{sample}.{type}", "out_pileup/cellSNP.base.vcf.gz"), 
               zip, 
               sample=samples_df['sample'], 
               type=samples_df['type'])

rule run_cellsnp:
    input:
        bam = lambda wc: os.path.join(samples_df.loc[wc.sample, 'path'], "possorted_bam.bam"),
        barcodes = lambda wc: os.path.join(samples_df.loc[wc.sample, 'path'], "filtered_peak_bc_matrix/barcodes.tsv")
    output:
        vcf = os.path.join(config["output_base"], "{sample}.{type}", "out_pileup/cellSNP.base.vcf.gz")
    params:
        outdir = lambda wc: os.path.join(config["output_base"], f"{wc.sample}.{wc.type}", "out_pileup"),
        vcfdb = config["vcfdb"],
        cellsnp = config["cellsnp"],
        minMAF = config["minMAF"],
        minCOUNT = config["minCOUNT"]
    threads: config["threads"]
    resources:
        # Snakemake 9 prefers 'mem_mb'. '8G' from your config is converted to 8000
        mem_mb = 8000,
        runtime = "12h" # You can now specify human-readable time
    shell:
        """
        {params.cellsnp} \
            -s {input.bam} \
            -b {input.barcodes} \
            -O {params.outdir} \
            -R {params.vcfdb} \
            -p {threads} \
            --minMAF {params.minMAF} \
            --minCOUNT {params.minCOUNT}
        """
