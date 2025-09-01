# Hua Sun
# Seurat v5.1
# v6.6  7/18/25


library(Seurat)
library(Signac)
library(JASPAR2022)
library(dplyr)
library(stringr)
library(stringi)
library(ggplot2)
library(seumultiome)
library(scqc)
library(GetoptLong)

options(future.globals.maxSize=2000000000000)

ref <- 'mm10'
regress <- NULL
outdir <- 'out_multiome_integrated'


# set
annotation <- Annotations(ref)
genome <- UCSCBSGenome(ref)



multiome_filtered_list <- readRDS(paste0(outdir,'/filtered.remDoublet_emptyDrops.rds'))
integ_atac <- readRDS(paste0(outdir,'/snatac_integ.withumap.rds'))



integrated_multiome <- MultiomeIntegration_v1(obj_list=multiome_filtered_list, regress=regress, snatac_obj=integ_atac)

# save integrated data
#saveRDS(integrated_multiome, file=paste0(outdir,'/multiome_integrated.sct_lsi.rds'))
#write.table(integrated_multiome@meta.data, file = paste0(outdir, '/multiome_integrated.sct_lsi.metadata.xls'), sep = "\t", quote=F, col.names = NA)

p <- DimPlot(integrated_multiome, reduction = 'wnn.umap', group.by = 'Sample', pt.size = 0.1)
pdf(paste0(outdir, '/multiome_integrated.sct_lsi.wnnUMAP.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()




#pfm <- readRDS(paste0(outdir, "/pfm.jaspar2022_vertebrates_core.rds"))E
#integrated_multiome <- AddGeneActivityAndMotif(obj=integrated_multiome, genome=genome, pfm=pfm, assay_use='peaks')
#saveRDS(integrated_multiome, file=paste0(outdir,'/multiome_integrated.sct_lsi.plus.rds'))


# Seurat v5
#integrated_multiome[['RNA']] <- JoinLayers(integrated_multiome[['RNA']])
#saveRDS(integrated_multiome, file=paste0(outdir,'/multiome_integrated.sct_lsi.plus2.rds'))








