# Hua Sun
# v6.2  1/24/25


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
saveRDS(multiome_list, file=paste0(outdir,'/rawdata.list.rds'))


# filter
print('[INFO] filter ...')
multiome_filtered_list <- lapply(X = multiome_list, FUN = function(x) { x <- FilterCellsMultiomeYML(x, yml) })
#saveRDS(multiome_filtered_list, file=paste0(outdir,'/filtered.list.rds'))


multiome_filtered_list <- lapply(X = multiome_filtered_list, FUN = function(x) { x <- FilterDoubletByDoubletFinder(x, outdir) })
saveRDS(multiome_filtered_list, file=paste0(outdir,'/snrna.remDoublet.rds'))


# extract top x feature cells
if (max_cells < 1000000){
    print('[INFO] extract cells ...')
    multiome_filtered_list <- base::lapply(X = multiome_filtered_list, FUN = function(x) { x <- ExtractNCells(x, ncells=max_cells) })
}


# remove processed
multiome_list <- NULL



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

# remove processed
snatac_macs2_list <- NULL


# snATAC: 3.Integration 
snatac_peak_list <- lapply(X = snatac_peak_list, FUN = function(x) { x <- ComputeLSI(x, 'peaks', 10) })
integ_atac <- RunIntegrationATAC(snatac_list=snatac_peak_list)
saveRDS(integ_atac, file=paste0(outdir,'/snatac_integ.withumap.rds'))
write.table(integ_atac@meta.data, file = paste0(outdir, '/snatac_integ.metadata.xls'), sep = "\t", quote=F, col.names = NA)

p <- DimPlot(integ_atac, group.by = 'Sample', reduction = "atac.umap", pt.size = 0.1)
pdf(paste0(outdir, '/snatac_integ.umap.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()



# Combind snATAC and snRNA
print('[INFO] integrate multiome ...')
integrated_multiome <- MultiModalNeighborsMultiomePipe(obj_list=multiome_filtered_list, snatac_obj=integ_atac)


# save integrated data
saveRDS(integrated_multiome, file=paste0(outdir,'/multiome_integrated.rds'))
write.table(integrated_multiome@meta.data, file = paste0(outdir, '/multiome_integrated.metadata.xls'), sep = "\t", quote=F, col.names = NA)

p <- DimPlot(integrated_multiome, reduction = 'wnn.umap', group.by = 'Sample', pt.size = 0.1)
pdf(paste0(outdir, '/multiome_integrated.wnnUMAP.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()




######################
##  add gene activity and motifs 
######################

#snatac_peak_list <- readRDS(paste0(outdir, '/snatac_peaks.list.rds'))
#integrated_multiome <- readRDS(paste0(outdir,'/multiome_integrated.rds'))

# snATAC: Add gene activity 
print('[INFO] call pfm ...')
pfm <- CallPFMFromJASPAR(JASPAR2022)
saveRDS(pfm, file = paste0(outdir, "/pfm.jaspar2022.rds"))

print('snatac_peakExtra_list')
snatac_peakExtra_list <- lapply(X = snatac_peak_list, FUN = function(x) { x <- AddGeneActivityAndMotif(obj=x, genome=genome, pfm=pfm) })
saveRDS(snatac_peakExtra_list, file=paste0(outdir, '/snatac_peaksExtra.list.rds'))

# remove processed
snatac_peak_list <- NULL

combined_snatac <- MergeObjects(snatac_peakExtra_list)

# remove processed
snatac_peakExtra_list <- NULL

integrated_multiome[['geneact']] <- combined_snatac[['geneact']]
integrated_multiome[['chromvar']] <- combined_snatac[['chromvar']]

# remove processed
combined_snatac <- NULL

saveRDS(integrated_multiome, file=paste0(outdir,'/multiome_integrated.plus.rds'))



