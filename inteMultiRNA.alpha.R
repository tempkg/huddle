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
ref <- 'mm10'
filter <- 'v3b'
min_cells <- 3
outdir <- 'out_snRNAFilterAlpha'

GetoptLong(
    "rdir=s",      "dir",
    "ref=s",       "ref",
    "filter=s",    "filter version",
    "min_cells=i", "min_cells",
    "outdir=s",    "outdir"
)


dir.create(outdir)


# set
annotation <- Annotations(ref)


# make object list
print('[INFO] make object ...')
#multiome_list <- DirCreateMultiomeObject(rdir, ref, annotation, min_cells, outdir)
#print(multiome_list)
#saveRDS(multiome_list, file=paste0(outdir,'/rawdata.list.rds'))

multiome_list <- readRDS(paste0(outdir,'/rawdata.list.rds'))

# filter
print('[INFO] filter ...')
multiome_filtered_list <- c()
if (filter == 'v2'){
    multiome_filtered_list <- lapply(X = multiome_list, FUN = function(x) { x <- FilterCellsMultiome_v2(x) })
}
if (filter == 'v2a'){
    multiome_filtered_list <- lapply(X = multiome_list, FUN = function(x) { x <- FilterCellsMultiome_v2a(x) })
}
if (filter == 'v2d'){
    multiome_filtered_list <- lapply(X = multiome_list, FUN = function(x) { x <- FilterCellsMultiome_v2d(x) })
}
if (filter == 'v2e'){
    cutoff <- CalculateMultiomeMaxCutoffN(multiome_list)
    multiome_filtered_list <- lapply(X = multiome_list, FUN = function(x) { x <- FilterCellsMultiome_v2e(x, cutoff[1], cutoff[2], cutoff[3]) })
}
if (filter == 'v2f'){
    multiome_filtered_list <- lapply(X = multiome_list, FUN = function(x) { x <- FilterCellsMultiome_v2f(x) })
}
if (filter == 'v3'){
    multiome_filtered_list <- lapply(X = multiome_list, FUN = function(x) { x <- FilterCellsMultiome_v3(x) })
}
if (filter == 'v3b'){
    cutoff <- CalMedianMaxCutoffN(multiome_list)
    print(cutoff)
    multiome_filtered_list <- lapply(X = multiome_list, FUN = function(x) { x <- FilterCellsMultiome_v3b(x, cutoff[1], cutoff[2], cutoff[3]) })
}
#saveRDS(multiome_filtered_list, file=paste0(outdir,'/filtered.list.rds'))


multiome_filtered_list <- lapply(X = multiome_filtered_list, FUN = function(x) { x <- FilterDoubletByDoubletFinder(x, outdir) })
#saveRDS(multiome_filtered_list, file=paste0(outdir,'/filtered.doublet.list.rds'))


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
#saveRDS(combined_snrna, file=paste0(outdir,'/snrna_merged.rds'))

# just for test
#merged_snrna <- combined_snrna %>%
#                RunPCA(assay = 'SCT', features=features) %>%
#                RunUMAP(assay = 'SCT', dims = 1:30) %>%
#                FindNeighbors(dims = 1:30) %>%
#                FindClusters(resolution=0.3)
#saveRDS(merged_snrna, file=paste0(outdir,'/snrna_merged.withumap.rds'))

# remove processed
snrna_list <- NULL
#merged_snrna <- NULL


print('[INFO] harmony snRNAs ...')
harmony_rna <- CustomizedHarmonyRNA1(combined_snrna, features)

# save harmony data
saveRDS(harmony_rna, file=paste0(outdir,'/snrna_harmony.rds'))
#write.table(harmony_rna@meta.data, file = paste0(outdir, '/snrna_harmony.metadata.xls'), sep = "\t", quote=F, col.names = NA)

# plot umap
p <- DimPlot(harmony_rna, group.by = 'Sample', reduction = "umap.rna", pt.size = 0.1)
pdf(paste0(outdir, '/snrna_harmony.umap.pdf'), width = 5.5, height = 4, useDingbats=FALSE)
print(p)
dev.off()







