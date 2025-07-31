
library(Seurat)
library(seumultiome)
library(dplyr)
library(ggplot2)


rds <- 'multiome_integrated.plus.rds'
finfo <- 'sample.group.txt'
outdir <- 'out_subgroup'
group <- 'ZR'

dir.create(outdir)


seu <- readRDS(rds)

info <- read.table(finfo, sep='\t', header=F)
colnames(info) <- c('sample', 'group')
seu$group <- info$group[match(seu$orig.ident, info$sample)]


# Seurat v5 (if RNA didn't process JoinLayers)
#seu[['RNA']] <- JoinLayers(seu[['RNA']])

#samples <- c('H_09_2459','H_09_2461','H_15_4174','H_16_05327','H_93_0050','H_99_0398')

# tumor cell
Idents(seu) <- 'group'
tumor_obj <- subset(seu, subset=group == 'ZR')

# re-umap
# default: dim_max_rna=30, dim_max_atac=30, res_clus=0.4, seed=42
tumor_obj <- RunWNNPlus(obj=tumor_obj, reduc_rna='pca', reduc_atac='integrated_lsi')

# output
p <- DimPlot(tumor_obj, group.by='orig.ident', reduction = 'wnn.umap')
ggsave(paste0(outdir, '/', group, '.umap.sample.pdf'), w=8, h=6)

#saveRDS(tumor_obj, paste0(outdir, '/epn_zr_tumor_cells.rds'))



