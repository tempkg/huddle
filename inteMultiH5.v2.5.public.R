# Hua Sun
# v2.5  12/27/24


library(Seurat)
library(Signac)
library(JASPAR2022)
library(JASPAR2024)
library(dplyr)
library(stringr)
library(stringi)
library(harmony)
library(GetoptLong)
library(SeuMultiome)

set.seed(47)
options(future.globals.maxSize=2000000000000)


rdir <- 'out_cellranger_arc'
ref <- 'hg38'
filter <- 'v2'
min_cells <- 3
outdir <- 'out_multiome_integrated'

GetoptLong(
    "rdir=s",      "dir",
    "ref=s",       "ref",
    "filter=s",    "filter version",
    "min_cells=i", "min_cells",
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
multiome_filtered_list <- c()
if (filter == 'v2'){
    multiome_filtered_list <- lapply(X = multiome_list, FUN = function(x) { x <- FilterCellsMultiome_v2(x) })
}
if (filter == 'v3'){
    multiome_filtered_list <- lapply(X = multiome_list, FUN = function(x) { x <- FilterCellsMultiome_v3(x) })
}
saveRDS(multiome_filtered_list, file=paste0(outdir,'/filtered.list.rds'))


multiome_filtered_list <- lapply(X = multiome_filtered_list, FUN = function(x) { x <- FilterDoubletByDoubletFinder(x, outdir) })
saveRDS(multiome_filtered_list, file=paste0(outdir,'/filtered.doublet.list.rds'))


# remove processed
multiome_list <- NULL



print('[INFO] normalize ...')
snrna_list <- base::lapply(X = multiome_filtered_list, FUN = function(x) { x <- SCTransform(x, vst.flavor = "v2", variable.features.n = 3000) })

# sct-integration & pca
features <- SelectIntegrationFeatures(object.list = snrna_list, nfeatures = 3000)
write.table(features, file = paste0(outdir, '/snrna_selectedIntegrationFeatures.xls'), sep = "\t", quote=F, col.names = NA)


# Merge normalized samples for harmony
print('[INFO] combine snRNAs ...')
combined_snrna <- ObjectListMergeObjects(snrna_list)
# save merged
saveRDS(combined_snrna, file=paste0(outdir,'/snrna_merged.rds'))




print('[INFO] harmony snRNAs ...')
harmony_rna <- CustomizedHarmonyRNA1(combined_snrna, features)

# save harmony data
saveRDS(harmony_rna, file=paste0(outdir,'/snrna_harmony.rds'))
write.table(harmony_rna@meta.data, file = paste0(outdir, '/snrna_harmony.metadata.xls'), sep = "\t", quote=F, col.names = NA)

# plot umap
p <- DimPlot(harmony_rna, group.by = 'Sample', reduction = "rna.umap", pt.size = 0.1)
pdf(paste0(outdir, '/snrna_harmony.umap.pdf'), width = 5.5, height = 4, useDingbats=FALSE)
print(p)
dev.off()


# remove processed
snrna_list <- NULL
combined_snrna <- NULL



print('[INFO] macs2 ...')
snatac_macs2_list <- lapply(X = multiome_filtered_list, FUN = function(x) { x <- CallPeaksUsingMACS2(x, macs2, annotation) })
saveRDS(snatac_macs2_list, file=paste0(outdir, '/snatac_macs2.list.rds'))



print('[INFO] combine peaks ...')
snatac_peak_list <- ObjectListConvertBedToGRanges(snatac_macs2_list)

saveRDS(snatac_peak_list, file=paste0(outdir, '/snatac_peaks.list.rds'))



print('[INFO] call pfm ...')
pfm <- CallPFMFromJASPAR(JASPAR2024)
saveRDS(pfm, file = paste0(outdir, "/pfm.jaspar2024.", ref, ".rds"))

print('snatac_peakExtra_list')
snatac_peakExtra_list <- lapply(X = snatac_peak_list, FUN = function(x) { x <- AddGeneActivityAndMotif(x, genome, pfm) })
saveRDS(snatac_peakExtra_list, file=paste0(outdir, '/snatac_peaksExtra.list.rds'))



print('[INFO] merge snATAC ...')
combined_snatac <- ObjectListMergeObjects(snatac_peakExtra_list)
saveRDS(combined_snatac, file=paste0(outdir, '/snatac_merged.rds'))

# remove processed
snatac_macs2_list <- NULL
snatac_peak_list <- NULL
snatac_peakExtra_list <- NULL



print('[INFO] harmony snATAC ...')
harmony_atac <- CustomizedHarmonyATAC1(combined_snatac)

# save harmony data
saveRDS(harmony_atac, file=paste0(outdir,'/snatac_harmony.rds'))
write.table(harmony_atac@meta.data, file = paste0(outdir, '/snatac_harmony.metadata.xls'), sep = "\t", quote=F, col.names = NA)

p <- DimPlot(harmony_atac, group.by = 'Sample', reduction = "atac_merged.umap", pt.size = 0.1)
pdf(paste0(outdir, '/snatac_merged.umap.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()

p <- DimPlot(harmony_atac, group.by = 'Sample', reduction = "atac.umap", pt.size = 0.1)
pdf(paste0(outdir, '/snatac_harmony.umap.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()



# remove processed
combined_snatac <- NULL



print('[INFO] integrate multiome ...')
integrated_multiome <- CustomizedCombinationForHarmonizedData1(harmony_atac, harmony_rna)

saveRDS(integrated_multiome, file=paste0(outdir,'/multiome_integrated.rds'))
write.table(integrated_multiome@meta.data, file = paste0(outdir, '/multiome_integrated.metadata.xls'), sep = "\t", quote=F, col.names = NA)

p <- DimPlot(integrated_multiome, reduction = 'wnn.umap', group.by = 'Sample', pt.size = 0.1)
pdf(paste0(outdir, '/multiome_integrated.wnnUMAP.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()







