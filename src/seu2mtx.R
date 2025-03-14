# 2024-11-15 v2.5

library(Seurat)
library(dplyr)
library(yaml)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)
print(args)

dyml <- yaml.load_file(args[1])
rds <- dyml$rds
assay <- dyml$assay
ctl <- dyml$ctl
form <- dyml$form
outdir <- dym$outdir

dir.create(outdir)


seurat_obj <- readRDS(rds)
metadata <- seurat_obj@meta.data

cell_anno <- metadata[,c('seurat_clusters', 'cell_type2', 'orig.ident')]
cell_anno$cell <- rownames(metadata)

ref_group <- strsplit(ctl, ',')[[1]]
print(ref_group)

# make cell anno data
if (form == 'cluster'){
    cell_anno$group <- cell_anno$seurat_clusters
    cell_anno$group <- paste0('C', cell_anno$group)
}

if (form == 'sample'){
    cell_anno$group <- cell_anno$orig.ident
}

if (form == 'celltype'){
    cell_anno$group <- cell_anno$cell_type2
}

cell_anno$group <- as.character(cell_anno$group)
cell_anno$group[cell_anno$cell_type2 %in% ref_group] <- 'Control'

# output cell anno data
cell_anno <- cell_anno[, c('cell', 'group')]
fwrite(cell_anno, file=paste0(outdir, '/cell.anno'), sep='\t', row.names=F, col.names=F)


# output exp count data (fast)
assay_count <- seurat_obj@assays[[assay]]@counts[,cell_anno$cell]
fwrite(as.matrix(assay_count), file=paste0(outdir, '/count.mtx'), sep = "\t", quote=F, row.names=F, col.names=F)
writeLines(row.names(assay_count), paste0(outdir, '/genes.txt'))
writeLines(colnames(assay_count), paste0(outdir, '/barcodes.txt'))


