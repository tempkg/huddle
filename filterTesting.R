# Hua Sun
# v2.8  1/22/25
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
outdir <- 'out_filterTesting'

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
snrna_list <- Seurat4RNANormalize(multiome_filtered_list)
features <- SelectIntegrationFeatures(object.list = snrna_list, nfeatures = 3000)

# Merge normalized samples
print('[INFO] combine snRNAs ...')
combined_snrna <- MergeObjects(snrna_list)
VariableFeatures(combined_snrna) <- features
# just for test
merged_snrna <- RunPCAToUMAP(obj=combined_snrna, assay_use='SCT', reduc_name='umap.mergedrna', res_clus=0.4)
saveRDS(merged_snrna, file=paste0(outdir,'/snrna_sct_merged.withumap.v4.rds'))



# v4 (not work in Seurat v5)
snrna_rpca <- Seurat4SCTIntegration(obj_list=snrna_list)
saveRDS(snrna_rpca, file=paste0(outdir,'/snrna_sct_integ_rpca.withumap.v4.rds'))


# v5
snrna_integ <- Seurat5SCTIntegration(obj_list=multiome_filtered_list, integ_method='cca')
saveRDS(snrna_integ, file=paste0(outdir,'/snrna_sct_integ_cca.withumap.v5.rds'))

snrna_integ <- Seurat5SCTIntegration(obj_list=multiome_filtered_list, integ_method='rpca')
saveRDS(snrna_integ, file=paste0(outdir,'/snrna_sct_integ_rpca.withumap.v5.rds'))






