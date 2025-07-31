# Hua Sun


library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(seumultiome)
library(GetoptLong)


set.seed(47)
options(future.globals.maxSize=2000000000000)


rds <- 'multiome_integrated.plus.v2.malignantTumorOnly.rds'
finfo <- 'sample.group.txt'
outdir <- 'out_re-integration'


dim_rna <- 30
dim_atac <- 30


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
seu$group <- info$group[match(seu$Sample, info$sample)]



# snRNA
seu <- SCTransform(seu, vst.flavor='v2', variable.features.n=3000)
seu <- RunPCA(seu, seed.use=seed, verbose = F)

reduc_rna <- 'harmony_rna'
seu <- harmony::RunHarmony(object=seu, group.by.vars = 'orig.ident', reduction = 'pca', assay.use = 'SCT', reduction.save=reduc_rna)


# snATAC
seu <- ATACIntegrationHarmonyViaObject(obj=seu, assay_use='peaks', new_reduc='harmony_atac', max_dim=dim_atac, seed=seed)

# wnn
multiome_wnn <- RunWNNPlus(obj=seu, reduc_rna=reduc_rna, reduc_atac='harmony_atac', dim_max_rna=dim_rna, dim_max_atac=dim_atac, weight_name='SCT.weight', res_clus=0.4, seed=seed)


# save integrated data
#saveRDS(multiome_wnn, file=paste0(outdir,'/multiome_integrate.', dim_rna, '_', dim_atac, '.rds'))

p <- DimPlot(multiome_wnn, reduction = 'wnn.umap', group.by = 'group', pt.size = 0.1)
ggsave(paste0(outdir, '/multiome_integrate.wnnUMAP.', dim_rna, '_', dim_atac, '.pdf'), width = 6.5, height = 5)
#pdf(paste0(outdir, '/multiome_integrate.wnnUMAP.', dim_rna, '_', dim_atac, '.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
#print(p)
#dev.off()





# extract subgroup
for (g in unique(info$group)){
    tarcells <- rownames(seu@meta.data)[seu$group == g]
    subseu <- subset(seu, cells=tarcells)
    subseu <- RunWNNPlus(obj=subseu, reduc_rna=reduc_rna, reduc_atac='harmony_atac', dim_max_rna=dim_rna, dim_max_atac=dim_atac, weight_name='SCT.weight', res_clus=0.4, seed=seed)

    saveRDS(subseu, file=paste0(outdir,'/', g, '.multiome_integrate.', dim_rna, '_', dim_atac, '.rds'))

    p <- DimPlot(subseu, reduction = 'wnn.umap', group.by = 'Sample', pt.size = 0.1)
    ggsave(paste0(outdir, '/', g, '.multiome_integrate.wnnUMAP.', dim_rna, '_', dim_atac, '.pdf'), width = 6.5, height = 5)

}














