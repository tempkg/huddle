# Hua Sun
# v2.6  1/16/25


library(Seurat)
library(Signac)
library(stringr)
library(stringi)
library(harmony)
library(GetoptLong)
library(SeuMultiome)

options(future.globals.maxSize=2000000000000)


rds <- 'out_rawobject/rawdata.list.rds'
ref <- 'mm10'
yml <- 'filter.yml'
outdir <- 'out_filterMerge'

GetoptLong(
    "rds=s",       "rds",
    "ref=s",       "ref",
    "yml=s",       "yml",
    "outdir=s",    "outdir"
)


dir.create(outdir)


# set
annotation <- Annotations(ref)

multiome_list <- readRDS(rds)

multiome_filtered_list <- lapply(X = multiome_list, FUN = function(x) { x <- FilterCellsMultiomeYML(x, yml) })
multiome_filtered_list <- lapply(X = multiome_filtered_list, FUN = function(x) { x <- FilterDoubletByDoubletFinder(x, outdir) })
saveRDS(multiome_filtered_list, file=paste0(outdir,'/snrna.remDoublet.rds'))

# remove processed
multiome_list <- NULL


print('[INFO] normalize ...')
snrna_list <- base::lapply(X = multiome_filtered_list, FUN = function(x) { x <- SCTransform(x, vst.flavor = "v2", variable.features.n = 3000) })

features <- SelectIntegrationFeatures(object.list = snrna_list, nfeatures = 3000)

# Merge normalized samples for harmony
print('[INFO] combine snRNAs ...')
combined_snrna <- ObjectListMergeObjects(snrna_list)

# just for test
merged_snrna <- combined_snrna %>%
                RunPCA(assay = 'SCT', features=features) %>%
                RunUMAP(assay = 'SCT', dims = 1:30, reduction.name = "umap.mergedrna") %>%
                FindNeighbors(dims = 1:30) %>%
                FindClusters(resolution=0.4)
saveRDS(merged_snrna, file=paste0(outdir,'/snrna_sct_merged.withumap.rds'))




#------- regress out mt
print('[INFO] normalize ...')
snrna_list <- base::lapply(X = multiome_filtered_list, FUN = function(x) { x <- SCTransform(x, vst.flavor = "v2", variable.features.n = 3000, vars.to.regress = "percent.mt") })

features <- SelectIntegrationFeatures(object.list = snrna_list, nfeatures = 3000)

# Merge normalized samples for harmony
print('[INFO] combine snRNAs ...')
combined_snrna <- ObjectListMergeObjects(snrna_list)

# just for test
merged_snrna <- combined_snrna %>%
                RunPCA(assay = 'SCT', features=features) %>%
                RunUMAP(assay = 'SCT', dims = 1:30, reduction.name = "umap.mergedrna") %>%
                FindNeighbors(dims = 1:30) %>%
                FindClusters(resolution=0.4)
saveRDS(merged_snrna, file=paste0(outdir,'/snrna_sct_merged.regressout_mt.withumap.rds'))





#------- normalize log
print('[INFO] normalize ...')
snrna_list <- base::lapply(X = multiome_filtered_list, FUN = function(x) { 
    x <- NormalizeData(x) 
    x <- FindVariableFeatures(x)
    x <- ScaleData(x)
})

features <- SelectIntegrationFeatures(object.list = snrna_list, nfeatures = 2000)

# Merge normalized samples for harmony
print('[INFO] combine snRNAs ...')
combined_snrna <- ObjectListMergeObjects(snrna_list)

# just for test
merged_snrna <- combined_snrna %>%
                RunPCA(features=features) %>%
                RunUMAP(dims = 1:30, reduction.name = "umap.mergedrna") %>%
                FindNeighbors(dims = 1:30) %>%
                FindClusters(resolution=0.4)
saveRDS(merged_snrna, file=paste0(outdir,'/snrna_log_merged.withumap.rds'))





#------- normalize log, regress out mt
print('[INFO] normalize ...')
snrna_list <- base::lapply(X = multiome_filtered_list, FUN = function(x) { 
    x <- NormalizeData(x) 
    x <- FindVariableFeatures(x)
    x <- ScaleData(x, vars.to.regress = "percent.mt")
})

features <- SelectIntegrationFeatures(object.list = snrna_list, nfeatures = 2000)

# Merge normalized samples for harmony
print('[INFO] combine snRNAs ...')
combined_snrna <- ObjectListMergeObjects(snrna_list)

# just for test
merged_snrna <- combined_snrna %>%
                RunPCA(features=features) %>%
                RunUMAP(dims = 1:30, reduction.name = "umap.mergedrna") %>%
                FindNeighbors(dims = 1:30) %>%
                FindClusters(resolution=0.4)
saveRDS(merged_snrna, file=paste0(outdir,'/snrna_log_merged.regressout_mt.withumap.rds'))













