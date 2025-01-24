# Hua Sun
# v5  1/22/25


library(Seurat)
library(Signac)
library(JASPAR2022)
library(dplyr)
library(stringr)
library(stringi)
library(GetoptLong)
library(SeuMultiome)

set.seed(47)
options(future.globals.maxSize=2000000000000)


rdir <- 'out_cellranger_arc'
ref <- 'mm10'
yml <- 'filter.yml'
min_cells <- 3
max_cells <- 1000000
outdir <- 'out_multiome_integrated'

GetoptLong(
    "rdir=s",      "dir",
    "ref=s",       "ref",
    "yml=s",       "yml",
    "min_cells=i", "min_cells",
    "max_cells=i", "set max cell counts",
    "outdir=s",    "outdir"
)


dir.create(outdir)


# set
macs2 <- '/research/groups/mackgrp/home/common/Software/miniconda3/envs/macs2/bin/macs2'
annotation <- Annotations(ref)
genome <- UCSCBSGenome(ref)

# make object list
print('[INFO] make object ...')
multiome_list <- DirCreateMultiomeObject(rdir, ref, annotation, min_cells, outdir)
print(multiome_list)
#saveRDS(multiome_list, file=paste0(outdir,'/rawdata.list.rds'))


# filter
print('[INFO] filter ...')
multiome_filtered_list <- lapply(X = multiome_list, FUN = function(x) { x <- FilterCellsMultiomeYML(x, yml) })
#saveRDS(multiome_filtered_list, file=paste0(outdir,'/filtered.list.rds'))


######################
##      snRNA 
######################

multiome_filtered_list <- lapply(X = multiome_filtered_list, FUN = function(x) { x <- FilterDoubletByDoubletFinder(x, outdir) })
saveRDS(multiome_filtered_list, file=paste0(outdir,'/snrna.remDoublet.rds'))


# remove processed
multiome_list <- NULL


# extract top x feature cells
if (max_cells < 1000000){
    print('[INFO] extract cells ...')
    multiome_filtered_list <- base::lapply(X = multiome_filtered_list, FUN = function(x) { x <- ExtractTopNCells(x, ncells=max_cells) })
}


# Integration
integ_rna <- Seurat5SCTIntegration(obj_list=multiome_filtered_list, integ_method='rpca', res_clus=0.5)
saveRDS(integ_rna, file=paste0(outdir,'/snrna_sct_integ.withumap.rds'))
write.table(integ_rna@meta.data, file = paste0(outdir, '/snrna_sct_integ.metadata.xls'), sep = "\t", quote=F, col.names = NA)


# plot umap
p <- DimPlot(integ_rna, group.by = 'Sample', reduction = "rna.umap", pt.size = 0.1)
pdf(paste0(outdir, '/snrna_sct_integ.umap.pdf'), width = 5.5, height = 4, useDingbats=FALSE)
print(p)
dev.off()



# remove processed
snrna_list <- NULL
combined_snrna <- NULL




######################
##      snATAC 
######################

# snATAC: 1.Call peaks using MACS2
print('[INFO] macs2 ...')
snatac_macs2_list <- lapply(X = multiome_filtered_list, FUN = function(x) { x <- CallPeaksUsingMACS2(x, macs2, annotation) })
saveRDS(snatac_macs2_list, file=paste0(outdir, '/snatac_macs2.list.rds'))


# snATAC: 2.Combine peaks & recall fragment
print('[INFO] combine peaks ...')
snatac_peak_list <- ObjectListConvertBedToGRanges(snatac_macs2_list)
saveRDS(snatac_peak_list, file=paste0(outdir, '/snatac_peaks.list.rds'))


# snATAC: 3.Add gene activity 
print('[INFO] call pfm ...')
pfm <- CallPFMFromJASPAR(JASPAR2022)
saveRDS(pfm, file = paste0(outdir, "/pfm.jaspar2022.", ref, ".rds"))

print('snatac_peakExtra_list')
snatac_peakExtra_list <- lapply(X = snatac_peak_list, FUN = function(x) { x <- AddGeneActivityAndMotifPlus(obj=x, genome=genome, pfm=pfm) })
saveRDS(snatac_peakExtra_list, file=paste0(outdir, '/snatac_peaksExtra.list.rds'))

# remove processed
snatac_macs2_list <- NULL
snatac_peak_list <- NULL


# snATAC: 4.Integration 
integ_atac <- RunIntegrationATAC(snatac_list=snatac_peakExtra_list)
saveRDS(integ_rna, file=paste0(outdir,'/snatac_integ.withumap.rds'))
write.table(integ_atac@meta.data, file = paste0(outdir, '/snatac_integ.metadata.xls'), sep = "\t", quote=F, col.names = NA)

p <- DimPlot(integ_atac, group.by = 'Sample', reduction = "atac.umap", pt.size = 0.1)
pdf(paste0(outdir, '/snatac_integ.umap.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()


# remove processed
snatac_peakExtra_list <- NULL
combined_snatac <- NULL




# Combind snRNA and snATAC 
print('[INFO] integrate multiome ...')
integrated_multiome <- CustomizedCombinationForIntegratedData(integ_rna=integ_rna, integ_atac=integ_atac)

# save integrated data
saveRDS(integrated_multiome, file=paste0(outdir,'/multiome_integrated.rds'))
write.table(integrated_multiome@meta.data, file = paste0(outdir, '/multiome_integrated.metadata.xls'), sep = "\t", quote=F, col.names = NA)

p <- DimPlot(integrated_multiome, reduction = 'wnn.umap', group.by = 'Sample', pt.size = 0.1)
pdf(paste0(outdir, '/multiome_integrated.wnnUMAP.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()







