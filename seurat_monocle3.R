library(GetoptLong)
library(Seurat)
library(pseudotime)
library(dplyr)
library(slingshot)
library(data.table)
library(monocle3)

rds <- 'sc_celltype_anno.rds'
outdir <- 'out_monocle3'

#GetoptLong(
#    "rds=s",       "matrix file",
#    "outdir=s",    "output path"
#)


dir.create(outdir)

seu <- readRDS(rds)

RunMonocle3Pipe(seu=seu, root_cells='M_EPN_IUE0_AATTTCCTCCCGTTAC-1', outdir=outdir)






