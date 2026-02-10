nextflow.enable.dsl = 2

params.input_tsv   = params.input_tsv
params.output_base = params.output_base
params.cellsnp     = params.cellsnp
params.vcfdb       = params.vcfdb
params.threads     = params.threads ?: 1
params.memory      = params.memory ?: '4 GB'
params.minMAF      = params.minMAF ?: 0
params.minCOUNT    = params.minCOUNT ?: 1

/*
 * Read TSV (no header)
 * col1 = sample
 * col2 = path
 */
Channel
    .fromPath(params.input_tsv)
    .splitCsv(sep: '\t', header: false)
    .map { row -> 
        tuple(row[0], row[1]) 
    }
    .set { samples_ch }


process CELLSNP {

    tag { sample }

    cpus params.threads
    memory params.memory

    publishDir "${params.output_base}", mode: 'copy'

    input:
    tuple val(sample), val(path)

    output:
    path "${sample}.cellsnp/**"

    script:
    """
    ${params.cellsnp} \
        -s ${path}/${sample}/outs/possorted_bam.bam \
        -b ${path}/${sample}/outs/filtered_peak_bc_matrix/barcodes.tsv \
        -O ${sample}.cellsnp/out_pileup \
        -R ${params.vcfdb} \
        -p ${task.cpus} \
        --minMAF ${params.minMAF} \
        --minCOUNT ${params.minCOUNT}
    """
}


workflow {
    samples_ch | CELLSNP
}
