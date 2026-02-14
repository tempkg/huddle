nextflow.enable.dsl=2


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
process CALLHLA {

    tag { sample }

    input:
    tuple val(sample), val(base_path)

    script:
    """
    set -euo pipefail

    # ---- Input validation ----
    bam="${base_path}/${sample}/outs/possorted_bam.bam"

    if [[ ! -f "\$bam" ]]; then
        echo "ERROR: BAM not found: \$bam" >&2
        exit 1
    fi

    # ---- Output directory ----
    outdir="${params.output_base}/${sample}.hla"
    mkdir -p "\$outdir"

    # ---- Run cellsnp-lite ----
    conda run -n hla arcasHLA extract "\$bam" -o "\$outdir" -t 8 -v
    """
}

/*
 * Workflow
 */
workflow {
    new File(params.output_base).mkdirs()
    CALLHLA(samples_ch)
}


