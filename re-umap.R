# Hua Sun
# 2025-08-04

library(Seurat)
library(Signac)
library(ggplot2)
library(seuproc)
library(GetoptLong)



rds <- 'scrna_sct2_integrated.harmony.rds'
outdir <- 'out_reumap'

dir.create(outdir)
seed <- 47

seu <- readRDS(rds)


#method='umap-learn' not run in general macbook
subseu <- RunUMAP_Plus(obj=seu, 
    method='uwot',
    #method='umap-learn',
    reduc='harmony_rna',
    reduc_name='umap',
    reduc_key='UMAP_',
    min_dim=1,
    max_dim=30,
    algorithm_clus=1,
    res_clus=0.4,
    seed=seed)


p <- DimPlot(seu, reduction = 'rna.umap', pt.size = 0.1, label=T)
ggsave(paste0(outdir, '/scRNA_integrated.re-umap.pdf'), width = 6.5, height = 5)


#saveRDS(seu, file=paste0(outdir,'/scRNA_integrated.re-umap.rds'))












