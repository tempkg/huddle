#!/usr/bin/env nextflow

// Simple cellsnp pipeline
// Usage: nextflow run cellsnp.nf

// Read samples from TSV
samples = Channel.fromPath("samples.tsv")
    .splitCsv(sep: '\t', header: false)
    .map { row -> [row[0], row[1]] }

// Simple process
process cellsnp {
    tag "${sample}"
    
    input:
    tuple val(sample), val(path)
    
    output:
    path("${sample}.pileup/**")
    
    script:
    """
    # Create output directory
    outdir="${params.output_base}/${sample}.pileup"
    mkdir -p \$outdir
    
    # Run cellsnp-lite
    ${params.cellsnp} \
        -s ${path}/${sample}/outs/possorted_bam.bam \
        -b ${path}/${sample}/outs/filtered_peak_bc_matrix/barcodes.tsv \
        -O \$outdir \
        -R ${params.vcfdb} \
        -p ${params.threads} \
        --minMAF ${params.minMAF} \
        --minCOUNT ${params.minCOUNT}
    """
}

workflow {
    cellsnp(samples)
}