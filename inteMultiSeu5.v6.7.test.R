# Hua Sun
# Seurat v5.1
# v6.7  7/29/25


library(Seurat)
library(Signac)
library(dplyr)
library(seumultiome)


set.seed(47)
options(future.globals.maxSize=2000000000000)


outdir <- 'out_multiome_integrated'

dir.create(outdir)



multiome_filtered_list <- readRDS('filtered.remDoublet.rds')
integ_atac <- readRDS('snatac_integ.withumap.rds')



snrna <- MergeObjects(multiome_filtered_list)
DefaultAssay(snrna) <- 'RNA'
snrna <- Seurat::SCTransform(snrna, vst.flavor = "v2", vars.to.regress = regress, variable.features.n = 2000)
snrna <- Seurat::RunPCA(snrna, seed.use = 47)
saveRDS(snrna, file=paste0(outdir,'/snrna_normalized_sct2.rds'))

multiome_filtered_list <- NULL


integ_rna <- Seurat5Integration(obj=snrna, integ_method='harmony', norm_method='SCT', new_reduc='harmony_rna')
integ_rna <- RunUMAP_Plus(obj=integ_rna, reduc='harmony_rna', reduc_name='rna.umap', reduc_key='rnaUMAP_', min_dim=1, max_dim=30, seed=47, res_clus=0.4)
saveRDS(integ_rna, file=paste0(outdir,'/snrna_integ.withumap.rds'))

p <- DimPlot(integ_rna, reduction = 'rna.umap', group.by = 'Sample', pt.size = 0.1)
pdf(paste0(outdir, '/snrna_integrated.wnnUMAP.sample.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()

p <- DimPlot(integ_rna, reduction = 'rna.umap', group.by = 'seurat_clusters', pt.size = 0.1)
pdf(paste0(outdir, '/snrna_integrated.wnnUMAP.cluster.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()




integrated_multiome <- RunWNN_Customized(integ_rna=integ_rna, integ_atac=integ_atac, 
                            reduc_rna='harmony_rna', reduc_atac='integrated_lsi', dim_max_rna=30, dim_max_atac=30, 
                            weight_name = "RNA.weight", res_clus=0.4, seed=47)



p <- DimPlot(integrated_multiome, reduction = 'wnn.umap', group.by = 'Sample', pt.size = 0.1)
pdf(paste0(outdir, '/multiome_integrated.wnnUMAP.sample.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()

p <- DimPlot(integrated_multiome, reduction = 'wnn.umap', group.by = 'seurat_clusters', pt.size = 0.1)
pdf(paste0(outdir, '/multiome_integrated.wnnUMAP.cluster.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()







