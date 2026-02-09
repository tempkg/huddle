nextflow.enable.dsl=2

params.config = './config.yml'

/*
 * Load YAML config
 */
def cfg = new groovy.yaml.YamlSlurper().parse(file(params.config))

Channel
    .fromPath(cfg.input_tsv)
    .splitText()
    .map { line ->
        def (sample, path) = line.tokenize()
        tuple(sample, path)
    }
    .set { samples_ch }


process CELLSNP {

    tag { "${sample}" }

    publishDir "${cfg.output_base}/${sample}.cellsnp",
               mode: 'copy'

    cpus cfg.threads
    memory cfg.memory

    input:
    tuple val(sample), val(path)

    output:
    path "out_pileup", emit: pileup

    script:
    """
    ${cfg.cellsnp} \
        -s ${path}/${sample}/outs/possorted_bam.bam \
        -b ${path}/${sample}/outs/filtered_peak_bc_matrix/barcodes.tsv \
        -O out_pileup \
        -R ${cfg.vcfdb} \
        -p ${cfg.threads} \
        --minMAF ${cfg.minMAF} \
        --minCOUNT ${cfg.minCOUNT}
    """
}


workflow {
    CELLSNP(samples_ch)
}
