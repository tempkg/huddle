# Snakefile
import pandas as pd

# Load configuration
configfile: "config.yml"

# Read sample file (no header)
samples = []
with open(config["sample_file"]) as f:
    for line in f:
        if line.strip():  # Skip empty lines
            parts = line.strip().split()
            if len(parts) >= 2:
                samples.append({"sample": parts[0], "path": parts[1]})

# Convert to pandas DataFrame for easier access
samples_df = pd.DataFrame(samples)
SAMPLES = list(samples_df['sample'])

# Main rule
rule all:
    input:
        expand("{output_base}/{sample}/cellSNP.base.vcf.gz", 
               output_base=config["output_base"],
               sample=SAMPLES)

# Processing rule
rule cellsnp:
    input:
        bam=lambda wildcards: f"{samples_df.loc[samples_df['sample']==wildcards.sample, 'path'].values[0]}/{wildcards.sample}/outs/possorted_bam.bam",
        barcodes=lambda wildcards: f"{samples_df.loc[samples_df['sample']==wildcards.sample, 'path'].values[0]}/{wildcards.sample}/outs/filtered_peak_bc_matrix/barcodes.tsv"
    output:
        directory("{output_base}/{sample}")
    threads: config["threads"]
    resources:
        mem_mb=int(str(config["memory"]).replace("G", "")) * 1024
    shell:
        """
        {config[cellsnp]} \
            -s {input.bam} \
            -b {input.barcodes} \
            -O {output} \
            -R {config[vcfdb]} \
            -p {threads} \
            --minMAF {config[minMAF]} \
            --minCOUNT {config[minCOUNT]}
        """