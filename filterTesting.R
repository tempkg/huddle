# Hua Sun
# v2.5  12/27/24


library(Seurat)
library(Signac)
library(stringr)
library(stringi)
library(GetoptLong)
library(SeuMultiome)

options(future.globals.maxSize=2000000000000)


rds <- 'out_rawobject/rawdata.list.rds'
ref <- 'mm10'
yml <- 'filtercutoff.yml'
outdir <- 'out_snRNAFilterAlpha'

GetoptLong(
    "fobj=s",      "fobj",
    "ref=s",       "ref",
    "yml=s",       "yml",
    "outdir=s",    "outdir"
)


dir.create(outdir)


# set
annotation <- Annotations(ref)

multiome_list <- readRDS(rds)

multiome_filtered_list <- lapply(X = multiome_list, FUN = function(x) { x <- FilterCellsMultiomeYML(x, yml) })
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
p <- DimPlot(harmony_rna, group.by = 'Sample', reduction = "rna.umap", pt.size = 0.1)
pdf(paste0(outdir, '/snrna_harmony.umap.pdf'), width = 5.5, height = 4, useDingbats=FALSE)
print(p)
dev.off()







