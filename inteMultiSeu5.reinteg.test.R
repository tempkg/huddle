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
outdir <- 'out_subgroup'


dim_rna <- 30
dim_atac <- 30


GetoptLong(
    "rds=s",       "rds",
    "outdir=s",    "outdir"
)



dir.create(outdir)

seu <- readRDS(rds)


seed <- 47

# snRNA
seu <- RunPCAToUMAP(obj=seu, assay_use='SCT', max_dim=dim_rna, seed=seed)
reduc_rna <- 'pca'


#seu <- Seurat5Integration(obj=seu, integ_method='harmony', norm_method='SCT', new_reduc='harmony_rna')
#reduc_rna <- 'harmony_rna'



# snATAC
seu <- ATACIntegrationHarmonyViaObject(obj=seu, assay_use='peaks', new_reduc='harmony_atac', max_dim=dim_atac, seed=seed)

# wnn
multiome_wnn <- RunWNNPlus(obj=seu, reduc_rna=reduc_rna, reduc_atac='harmony_atac', dim_max_rna=dim_rna, dim_max_atac=dim_atac, weight_name='SCT.weight', res_clus=0.4, seed=seed)


# save integrated data
saveRDS(multiome_wnn, file=paste0(outdir,'/multiome_integrate.', dim_rna, '_', dim_atac, '.rds'))

p <- DimPlot(multiome_wnn, reduction = 'wnn.umap', group.by = 'Sample', pt.size = 0.1)
pdf(paste0(outdir, '/multiome_integrate.wnnUMAP.', dim_rna, '_', dim_atac, '.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()




