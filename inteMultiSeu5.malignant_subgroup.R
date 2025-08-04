
library(Seurat)
library(seumultiome)
library(dplyr)
library(ggplot2)


rds <- 'multiome_integrated.plus.rds'
finfo <- 'sample.group.txt'
outdir <- 'out_subgroup.method1'

dim_rna <- 30
dim_atac <- 20
seed <- 47


dir.create(outdir)


seu <- readRDS(rds)




# read malignant cells
malignant_cells <- readLines(fcell)

info <- read.table(finfo, sep='\t', header=F)
colnames(info) <- c('sample', 'group')

seu <- readRDS(rds)
seu$group <- info$group[match(seu$orig.ident, info$sample)]
Idents(seu) <- 'group'

# 
seu$UMAP_1 <- seu@reductions[['wnn.umap']]@cell.embeddings[,1]
seu$UMAP_2 <- seu@reductions[['wnn.umap']]@cell.embeddings[,2]
# filter noise by umap
filtered_noise <- rownames(seu@meta.data[seu$UMAP_1 < 7 & seu$UMAP_2 < 7 & seu@UMAP_2 > -7,])
# final malignant cells
final_cells <- intersect(malignant_cells, filtered_noise)


# final obj
malignant_obj <- subset(seu, cells=final_cells)


# extract subgroup
#for (g in info$group){
for (g in c('ST', 'ZR', 'PF')){
    if (g == 'ST'){  
        subseu <-  subset(malignant_obj, subset=group %in% c('EWP', 'ZR'))
    } else {
        subseu <-  subset(malignant_obj, subset=group == g)
    }

    subseu <-  subset(malignant_obj, subset=group %in% g)
    subseu <- RunWNNPlus(obj=subseu, reduc_rna='pca', reduc_atac='integrated_lsi', dim_max_rna=dim_rna, dim_max_atac=dim_atac, res_clus=0.4, seed=seed)

    saveRDS(subseu, file=paste0(outdir,'/', g, '.multiome_integrate.', dim_rna, '_', dim_atac, '.rds'))

    p <- DimPlot(subseu, reduction = 'wnn.umap', group.by = 'Sample', pt.size = 0.1)
    ggsave(paste0(outdir, '/', g, '.multiome_integrate.wnnUMAP.', dim_rna, '_', dim_atac, '.pdf'), width = 6.5, height = 5)

}








