# Hua Sun


library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(seumultiome)
library(GetoptLong)
library(sctransform)


set.seed(47)
options(future.globals.maxSize=2000000000000)


rds <- 'multiome_integrated.plus.v2.malignantTumorOnly.rds'
finfo <- 'sample.group.txt'
outdir <- 'out_subgroup.re-integration'


dim_rna <- 30
dim_atac <- 20


GetoptLong(
    "rds=s",       "rds",
    "finfo=s",      "info group",
    "outdir=s",    "outdir"
)


dir.create(outdir)
seed <- 47

seu <- readRDS(rds)

info <- read.table(finfo, sep='\t', header=F)
colnames(info) <- c('sample', 'group')
seu$group <- info$group[match(seu$orig.ident, info$sample)]



ReRunFromNormalizeToIntegration <- function(obj, group, dim_rna, dim_atac, seed, outdir){
     # snRNA
    obj <- SCTransform(obj, vst.flavor='v2', variable.features.n=3000)
    obj <- RunPCA(obj, seed.use=seed, verbose = F)

    reduc_rna <- 'harmony_rna'
    obj <- harmony::RunHarmony(object=obj, group.by.vars = 'orig.ident', reduction = 'pca', assay.use = 'SCT', reduction.save=reduc_rna)


    # snATAC
    obj <- ATACIntegrationHarmonyViaObject(obj=obj, assay_use='peaks', new_reduc='harmony_atac', max_dim=dim_atac, seed=seed)

    # wnn
    obj <- RunWNNPlus(obj=obj, reduc_rna=reduc_rna, reduc_atac='harmony_atac', dim_max_rna=dim_rna, dim_max_atac=dim_atac, weight_name='SCT.weight', res_clus=0.4, seed=seed)

    saveRDS(obj, file=paste0(outdir,'/', group, '.multiome_integrate.', dim_rna, '_', dim_atac, '.rds'))

    p <- DimPlot(obj, reduction = 'wnn.umap', group.by = 'orig.ident', pt.size = 0.1)
    ggsave(paste0(outdir, '/', group, '.multiome_integrate.wnnUMAP.', dim_rna, '_', dim_atac, '.pdf'), width = 6.5, height = 5)
}





# extract subgroup
for (g in unique(info$group)){
    tarcells <- rownames(seu@meta.data)[seu$group == g]
    subseu <- subset(seu, cells=tarcells)
    
    ReRunFromNormalizeToIntegration(obj=subseu, group=g, dim_rna=dim_rna, dim_atac=dim_atac, seed=seed, outdir=outdir)

}












