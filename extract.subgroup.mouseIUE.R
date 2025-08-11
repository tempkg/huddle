# Mouse IUE

library(Seurat)
library(seumultiome)
library(dplyr)
library(ggplot2)


rds <- 'integrated.seurat5.1_v6.2/multiome_integrated.plus.rds'
outdir <- 'out_subgroup_tumor_celltype'

dir.create(outdir)


seu <- readRDS(rds)
meta <- seu@meta.data

meta$group <- 'GBM'
meta$group[meta$Sample %in% c('Y38','Y55','Y56')] <- 'YAP1'
meta$group[meta$Sample %in% c('ZR0','ZR1','ZR2')] <- 'ZR'
seu$group <- meta$group

yap1_tumor_cells <- rownames(meta %>% filter(seurat_clusters %in% c(3,14)))
gbm_tumor_cells <- rownames(meta %>% filter(seurat_clusters %in% c(1,5,12)))
zr_tumor_cells <- rownames(meta %>% filter(group == 'ZR') %>% filter(seurat_clusters %in% c(2,6,19,23)))


# all tumor cell data
tumor_obj <- RunWNNPlus(obj=seu, reduc_rna='pca', reduc_atac='integrated.atac', dim_max_rna=30, dim_max_atac=30, res_clus=0.4, seed=42)
saveRDS(tumor_obj, paste0(outdir, '/zr_yap1_gbm_tumor_cells.rds'))

p <- DimPlot(tumor_obj, reduction = 'wnn.umap', group.by = 'group', pt.size = 0.1)
ggsave(paste0(outdir, '/zr_yap1_gbm_tumor_cells.group.pdf'), width = 6.5, height = 5)



# EPN-ZR-Fus
zr_tumor_obj <- subset(seu, cells = zr_tumor_cells)
zr_tumor_obj <- RunWNNPlus(obj=zr_tumor_obj, reduc_rna='pca', reduc_atac='integrated.atac', dim_max_rna=30, dim_max_atac=30, res_clus=0.4, seed=42)
saveRDS(zr_tumor_obj, paste0(outdir, '/epn_zr_tumor_cells.rds'))




