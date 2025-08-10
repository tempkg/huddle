library(GetoptLong)
library(Seurat)
library(pseudotime)
library(dplyr)
library(slingshot)
library(data.table)

rds <- ''
fmarker <- 'findmarkers.xls'
outdir <- 'out_slingshot'

#GetoptLong(
#    "rds=s",       "matrix file",
#    "outdir=s",    "output path"
#)


dir.create(outdir)

seu <- readRDS(rds)
d_marker <- fread(fmarker, data.table=F)

# reduction - check 'seu' for rna umap name.
RunSlingshotPipe_Seurat(seu=seu, d_marker=d_marker, 
    assay='SCT', reduction='UMAP.RNA', 
    root='5', clus='seurat_clusters', 
    nknots = 10, lineage = 1, km=4, 
    w=4.5, h=3, outdir=outdir)






