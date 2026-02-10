#!/usr/bin/env nextflow

params.input_tsv = "./samples.tsv"
params.output_base = "/research/groups/mackgrp/home/common/Analysis/Project_PrimaryRelapse/snatac/out_cellsnps"
params.config = "./config.yml"

// Read config from YAML file
config = file(params.config).withReader { 
    new groovy.yaml.YamlSlurper().parse(it) 
}

// Read samples
samples = file(params.input_tsv).readLines()
    .collect { it.split('\t') }
    .collectEntries { [it[0], it[1]] }

workflow {
    // Process each sample
    samples.each { sample, path ->
        processCellSNP(sample, path)
    }
}

process processCellSNP {
    tag "$sample"
    
    input:
    val sample
    val path
    
    output:
    path("${sample}.pileup/**")
    
    script:
    """
    mkdir -p ${params.output_base}/${sample}.pileup
    
    ${config.cellsnp} \
        -s $path/$sample/outs/possorted_bam.bam \
        -b $path/$sample/outs/filtered_peak_bc_matrix/barcodes.tsv \
        -O ${params.output_base}/${sample}.pileup \
        -R ${config.vcfdb} \
        -p ${config.threads} \
        --minMAF ${config.minMAF} \
        --minCOUNT ${config.minCOUNT}
    """
}