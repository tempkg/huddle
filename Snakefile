# Snakefile - Snakemake v9
configfile: "config.yml"

# Read sample file (no header)
with open(config["sample_file"]) as f:
    SAMPLES = [line.strip().split() for line in f if line.strip()]

# Convert to dict: {sample: path}
SAMPLE_PATHS = {sample: path for sample, path in SAMPLES}

# Main rule
rule all:
    input:
        expand("{output_base}/{sample}/cellSNP.base.vcf.gz", 
               output_base=config["output_base"],
               sample=SAMPLE_PATHS.keys())

# Processing rule
rule cellsnp:
    input:
        bam=lambda w: f"{SAMPLE_PATHS[w.sample]}/{w.sample}/outs/possorted_bam.bam",
        barcodes=lambda w: f"{SAMPLE_PATHS[w.sample]}/{w.sample}/outs/filtered_peak_bc_matrix/barcodes.tsv"
    output:
        directory("{output_base}/{sample}")
    resources:
        cpus=config["threads"],
        mem_mb=int(config["memory"].rstrip("G")) * 1024
    shell:
        """
        {config[cellsnp]} \
            -s {input.bam} \
            -b {input.barcodes} \
            -O {output} \
            -R {config[vcfdb]} \
            -p {resources.cpus} \
            --minMAF {config[minMAF]} \
            --minCOUNT {config[minCOUNT]}
        """