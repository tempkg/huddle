# Hua Sun
# v2.6  1/16/25
# seurat v4 and v5


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
VariableFeatures(combined_snrna) <- features

# just for test
merged_snrna <- RunPCAToUMAP(obj=combined_snrna, assay_use='SCT', reduc_name='umap.mergedrna', res_clus=0.4)
saveRDS(merged_snrna, file=paste0(outdir,'/snrna_sct_merged.withumap.rds'))



snrna_cca <- RunCCAIntegrationForSCT(combined_snrna)
saveRDS(snrna_cca, file=paste0(outdir,'/snrna_sct_integ_cca.withumap.rds'))





