nextflow.enable.dsl=2



// Pick the experiment based on seq_type
def seq_params = params.seq_type.find { key, val -> params.seq in val.aliases }?.value

if (!seq_params) {
    error "Invalid seq_type: ${params.seq}. Must be one of: ${params.seq_type.collect { it.value.aliases }.flatten()}"
}

def celltag = seq_params.celltag
def umitag  = seq_params.umitag
def minMaf = seq_params.minMaf
def minCount = seq_params.minCount

println "Sequencing type: ${params.seq}"
println "Cell tag: ${celltag}"
println "UMI tag: ${umitag}"
println "Min MAF: ${minMaf}"
println "Min COUNT: ${minCount}"


/*
 * Load samples from TSV (no header)
 * column 1: sample name
 * column 2: base path to CellRanger outputs
 */
Channel
    .fromPath(params.input_tsv, checkIfExists: true)
    .splitCsv(sep: '\t', header: false)
    .map { row -> tuple(row[0], row[1]) }
    .set { samples_ch }

/*
 * cellsnp-lite process
 */
process CELLSNP {

    tag { sample }

    cpus { params.threads }

    input:
    tuple val(sample), val(base_path)

    output:
    path("${params.output_base}/${sample}.pileup")

    script:
    """
    set -euo pipefail

    # ---- Input validation ----
    bam="${base_path}/${sample}/outs/possorted_bam.bam"
    barcodes="${base_path}/${sample}/outs/filtered_peak_bc_matrix/barcodes.tsv"

    if [[ ! -f "\$bam" ]]; then
        echo "ERROR: BAM not found: \$bam" >&2
        exit 1
    fi

    if [[ ! -f "\$barcodes" ]]; then
        echo "ERROR: barcodes not found: \$barcodes" >&2
        exit 1
    fi

    # ---- Output directory ----
    outdir="${params.output_base}/${sample}.pileup"
    mkdir -p "\$outdir"


    # ---- Run cellsnp-lite ----
    ${params.cellsnp} \
        -s "\$bam" \
        -b "\$barcodes" \
        -O "\$outdir" \
        -p ${params.threads} \
        --cellTAG ${celltag} \
        --UMItag ${umitag} \
        --minMAF ${minMaf} \
        --minCOUNT ${minCount}
    """
}

/*
 * Workflow
 */
workflow {
    CELLSNP(samples_ch)
}
