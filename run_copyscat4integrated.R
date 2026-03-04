
library(CopyscAT)
library(scvar)
library(data.table)
library(dplyr)



RunCopyscAT4IntegrationData(
    ref_dir = 'hg38_references',
    py_pfrag = 'software/script/process_fragment_file.py',
    fragment = 'merged_atac_fragments.tsv.gz',
    fbc = 'cell_info.txt',
    fctrl = 'immune_cells.txt',
    outdir = 'out_copyscat'
)






