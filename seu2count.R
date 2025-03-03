
library(Seurat)
library(data.table)


rds <- 'multiome_integrated.plus.rds'
fmeta <- 'out_celltype/metaData.cellType.xls'
outdir <- 'out_exp_mtx'
dir.create(outdir)

seu <- readRDS(rds)
meta <- fread(fmeta, data.table=F)
seu$cell_type2 <- meta$cell_type2[match(rownames(seu@meta.data), meta$V1)]

for (query in unique(meta$cell_type2)){
    query_cell <- rownames(meta[meta[[ident]]==query, ])
    sub_obj <- subset(seu, cells = query_cell)
    sub_count <- as.matrix(sub_obj@assays[['RNA']]@counts)
    write.table(sub_count, paste0(outdir, '/', query, '.exp.txt'), sep='\t', quote=F, col.names=NA)
}


