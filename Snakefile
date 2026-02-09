import pandas as pd
import os

# Load configuration
configfile: "config.yml"

# Read headerless TSV: Column 0 is sample ID, Column 1 is the base path
# We assign names manually since there is no title row
samples_df = pd.read_csv(
    config["input_tsv"], 
    sep="\s+", 
    header=None, 
    names=["sample", "path"]
).set_index("sample", drop=False)

rule all:
    input:
        # Assuming we still want to keep the {sample}.{type} structure in output
        # If 'type' isn't in your TSV anymore, you can hardcode it or add it to the TSV
        expand(os.path.join(config["output_base"], "{sample}", "out_pileup/cellSNP.base.vcf.gz"), 
               sample=samples_df['sample'])

rule run_cellsnp:
    input:
        # Constructed as requested: $path/$sample/outs/...
        bam = lambda wc: os.path.join(samples_df.loc[wc.sample, 'path'], wc.sample, "outs/possorted_bam.bam"),
        barcodes = lambda wc: os.path.join(samples_df.loc[wc.sample, 'path'], wc.sample, "outs/filtered_peak_bc_matrix/barcodes.tsv")
    output:
        vcf = os.path.join(config["output_base"], "{sample}", "out_pileup/cellSNP.base.vcf.gz")
    params:
        outdir = lambda wc: os.path.join(config["output_base"], wc.sample, "out_pileup"),
        vcfdb = config["vcfdb"],
        cellsnp = config["cellsnp"],
        maf = config["minMAF"],
        count = config["minCOUNT"]
    threads: config["threads"]
    resources:
        # Snakemake 9 uses these to automatically communicate with bsub
        mem_mb = 8000,
        runtime = "24h",
        slurm_partition = "large_mem" # Or lsf_queue if using LSF
    shell:
        """
        {params.cellsnp} \
            -s {input.bam} \
            -b {input.barcodes} \
            -O {params.outdir} \
            -R {params.vcfdb} \
            -p {threads} \
            --minMAF {params.maf} \
            --minCOUNT {params.count}
        """