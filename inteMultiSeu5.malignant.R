
library(Seurat)
library(seumultiome)
library(dplyr)
library(ggplot2)


rds <- 'multiome_integrated_plus.regout_mt.rds'
finfo <- 'sample.group.txt'
outdir <- 'out_malignant_cells'

dim_rna <- 30
dim_atac <- 30
seed <- 47


dir.create(outdir)


seu <- readRDS(rds)

info <- read.table(finfo, sep='\t', header=F)
colnames(info) <- c('sample', 'group')
seu$group <- info$group[match(seu$orig.ident, info$sample)]
Idents(seu) <- 'group'

# 1
malignant_cells <- readLines(fcell)
# 2
seu$UMAP_1 <- seu@reductions[['wnn.umap']]@cell.embeddings[,1]
seu$UMAP_2 <- seu@reductions[['wnn.umap']]@cell.embeddings[,2]
# filter noise by umap
filtered_noise <- renames(seu@meta.data[seu$UMAP_1 <7 & seu$UMAP_2 <7 & seu@UMAP_2 > -7,])
# final malignant cells
final_cells <- intersect(malignant_cells, filtered_noise)

# final obj
malignant_obj <- subset(seu, cells=final_cells)


# extract subgroup
# 'reduc_rna' depend on original integration multiome pipeline setting
seu <- RunWNNPlus(obj=malignant_obj, reduc_rna='pca', reduc_atac='integrated_lsi', dim_max_rna=dim_rna, dim_max_atac=dim_atac, res_clus=0.4, seed=seed)
saveRDS(seu, file=paste0(outdir,'/2.multiome_integrated_plus.regout_mt.malignantCellOnly.WNN.', dim_rna, '_', dim_atac, '.rds'))

p <- DimPlot(seu, reduction = 'wnn.umap', group.by = 'group', pt.size = 0.1)
ggsave(paste0(outdir, '/2.multiome_integrated_plus.regout_mt.malignantCellOnly.WNN.', dim_rna, '_', dim_atac, '.pdf'), width = 6.5, height = 5)









