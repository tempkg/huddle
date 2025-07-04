# Hua Sun


library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(seumultiome)
library(GetoptLong)


set.seed(47)
options(future.globals.maxSize=2000000000000)


rds <- 'multiome_integrated_plus.seurat5.rds'
approach <- 'dim30'
outdir <- 'out_subgroup'


GetoptLong(
    "rds=s",       "rds",
    "approach=s",  "approach",
    "outdir=s",    "outdir"
)



dir.create(outdir)

seu <- readRDS(rds)


seed <- 47

# snRNA
seu <- RunPCAToUMAP(obj=seu, assay_use='SCT', max_dim=30, seed=seed)
reduc_rna <- 'pca'

# approach=2
if (approach == 2){
    seu <- Seurat5Integration(obj=seu, integ_method='harmony', norm_method='SCT', new_reduc='harmony_rna')
    reduc_rna <- 'harmony_rna'
}


# snATAC
seu <- ATACIntegrationHarmonyViaObject(obj=seu, assay_use='peaks', new_reduc='harmony_atac', max_dim=30, seed=seed)

# wnn
multiome_wnn <- RunWNNPlus(obj=seu, reduc_rna=reduc_rna, reduc_atac='harmony_atac', dim_max_rna=30, dim_max_atac=30, weight_name='SCT.weight', res_clus=0.4, seed=seed)


# save integrated data
saveRDS(multiome_wnn, file=paste0(outdir,'/multiome_integrate.', approach, '.rds'))

p <- DimPlot(multiome_wnn, reduction = 'wnn.umap', group.by = 'Sample', pt.size = 0.1)
pdf(paste0(outdir, '/multiome_integrate.wnnUMAP.', approach, '.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()




