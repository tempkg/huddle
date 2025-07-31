
library(Seurat)
library(seumultiome)
library(dplyr)
library(ggplot2)


rds <- 'multiome_integrated.plus.rds'
finfo <- 'sample.group.txt'
outdir <- 'out_subgroup'

dim_rna <- 30
dim_atac <- 30
seed <- 47


dir.create(outdir)


seu <- readRDS(rds)

info <- read.table(finfo, sep='\t', header=F)
colnames(info) <- c('sample', 'group')
seu$group <- info$group[match(seu$orig.ident, info$sample)]
Idents(seu) <- 'group'

# extract subgroup
for (g in unique(info$group)){
    subseu <-  subset(seu, subset=group == g)
    subseu <- RunWNNPlus(obj=subseu, reduc_rna='pca', reduc_atac='integrated_lsi', dim_max_rna=dim_rna, dim_max_atac=dim_atac, res_clus=0.4, seed=seed)

    saveRDS(subseu, file=paste0(outdir,'/', g, '.multiome_integrate.', dim_rna, '_', dim_atac, '.rds'))

    p <- DimPlot(subseu, reduction = 'wnn.umap', group.by = 'Sample', pt.size = 0.1)
    ggsave(paste0(outdir, '/', g, '.multiome_integrate.wnnUMAP.', dim_rna, '_', dim_atac, '.pdf'), width = 6.5, height = 5)

}








