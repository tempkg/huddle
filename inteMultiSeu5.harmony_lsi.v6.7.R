# Hua Sun
# Seurat v5.1
# v6.7  7/29/25


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


rdir <- 'out_cellranger_arc'
seed <- 42
ref <- 'hg38'
yml <- 'filter.yml'
min_cells <- 3
max_cells <- 20000
regress <- NULL
outdir <- 'out_multiome_integrated.seurat5.harmony_lsi'

rm_cells <- ''
fpfm <- ''

GetoptLong(
    "rdir=s",      "dir",
    "seed=i",      "seed",
    "ref=s",       "ref",
    "yml=s",       "yml",
    "min_cells=i", "min_cells",
    "max_cells=i", "set max cell counts",
    "regress=s",   "regress out",
    "rm_cells=s",  "remove target cells",
    "fpfm=s",       "pfm data",
    "outdir=s",    "outdir"
)


dir.create(outdir)

set.seed(seed)


# set
macs2 <- '/research/groups/mackgrp/home/common/Software/miniconda3/envs/macs2/bin/macs2'
annotation <- Annotations(ref)
genome <- UCSCBSGenome(ref)


# save results for singl sample
outdir_process <- paste0(outdir, '/processed')
dir.create(outdir_process)

# make object list
print('[INFO] make object ...')
multiome_list <- DirCreateMultiomeObject(rdir, ref, annotation, min_cells, outdir_process)
print(multiome_list)
saveRDS(multiome_list, file=paste0(outdir,'/rawdata.list.rds'))


# filter
print('[INFO] filter ...')
multiome_filtered_list <- lapply(X = multiome_list, FUN = function(x) { x <- FilterCellsMultiomeYML(x, yml) })
#saveRDS(multiome_filtered_list, file=paste0(outdir,'/filtered.list.rds'))

# filter doublet
multiome_filtered_list <- lapply(X = multiome_filtered_list, FUN = function(x) { x <- FilterDoubletByDoubletFinder(x, outdir_process) })
saveRDS(multiome_filtered_list, file=paste0(outdir,'/filtered.remDoublet.rds'))

# filter empty drop
multiome_filtered_list <- lapply(X = multiome_filtered_list, FUN = function(x) { 
    name <- unique(x$orig.ident)
    x <- RunFilterEmptyDrops(obj=x, dir=paste0(rdir, '/', name, '/outs'), form='cellranger-arc', outdir=paste0(outdir_process, '/emptydrop.', name)) 
})
saveRDS(multiome_filtered_list, file=paste0(outdir,'/filtered.remDoublet_emptyDrops.rds'))


# remove processed
multiome_list <- NULL




# filter - multiplet by atac (optional)
if (nchar(rm_cells) > 0){
    poorcells <- read.table(rm_cells, sep='\t', header=F)
    colnames(poorcells) <- c('cell', 'sample')
    multiome_filtered_list <- base::lapply(X = multiome_filtered_list, FUN = function(x) { x <- FilterCellsByTable(obj=x, df=poorcells) })
    saveRDS(multiome_filtered_list, file=paste0(outdir,'/filtered.remDoublet_emptyDrops_multiplet.rds'))
}



# limit maximum clean cell (optional)
if (max_cells < 20000){
    print('[INFO] extract cells ...')
    multiome_filtered_list <- base::lapply(X = multiome_filtered_list, FUN = function(x) { x <- ExtractNCells(x, ncells=max_cells) })
    saveRDS(multiome_filtered_list, file=paste0(outdir,'/filtered.remDoublet_emptyDrops.', max_cells, '.rds'))
}



# sct normalize
snrna <- MergeObjects(multiome_filtered_list)
DefaultAssay(snrna) <- 'RNA'
snrna <- Seurat::SCTransform(snrna, vst.flavor = "v2", vars.to.regress = regress, variable.features.n = 2000)
snrna <- Seurat::RunPCA(snrna, seed.use = seed)
saveRDS(snrna, file=paste0(outdir,'/snrna_normalized_sct2.rds'))

# integ snrna
integ_rna <- Seurat5Integration(obj=snrna, integ_method='harmony', norm_method='SCT', new_reduc='harmony_rna')
integ_rna <- RunUMAP_Plus(obj=integ_rna, reduc='harmony_rna', reduc_name='rna.umap', reduc_key='rnaUMAP_', min_dim=1, max_dim=30, seed=seed, res_clus=0.4)
saveRDS(integ_rna, file=paste0(outdir,'/snrna_integ.withumap.rds'))


p <- DimPlot(integ_rna, reduction = 'rna.umap', group.by = 'Sample', pt.size = 0.1)
pdf(paste0(outdir, '/snrna_integ.umap.sample.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()

p <- DimPlot(integ_rna, reduction = 'rna.umap', group.by = 'seurat_clusters', pt.size = 0.1)
pdf(paste0(outdir, '/snrna_integ.umap.cluster.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()



######################
##      snATAC 
######################

# snATAC: 1.Call peaks using MACS2
print('[INFO] macs2 ...')
snatac_macs2_list <- lapply(X = multiome_filtered_list, FUN = function(x) { x <- CallPeaksUsingMACS2(x, macs2, annotation) })
saveRDS(snatac_macs2_list, file=paste0(outdir, '/snatac_macs2.list.rds'))

# remove processed
multiome_filtered_list <- NULL


# snATAC: 2.Combine peaks & recall fragment
print('[INFO] combine peaks ...')
snatac_peak_list <- ObjectListConvertBedToGRanges(snatac_macs2_list)
snatac_peak_list <- lapply(X = snatac_peak_list, FUN = function(x) { x <- ComputeLSI(x, 'peaks', 10) })
saveRDS(snatac_peak_list, file=paste0(outdir, '/snatac_peaks.list.rds'))

# remove processed
snatac_macs2_list <- NULL


# snATAC: 3.Integration 
integ_atac <- ATACIntegration(snatac_list=snatac_peak_list, max_dim=30, seed=seed)
saveRDS(integ_atac, file=paste0(outdir,'/snatac_integ.withumap.rds'))
write.table(integ_atac@meta.data, file = paste0(outdir, '/snatac_integ.metadata.xls'), sep = "\t", quote=F, col.names = NA)

p <- DimPlot(integ_atac, group.by = 'Sample', reduction = "atac.umap", pt.size = 0.1)
pdf(paste0(outdir, '/snatac_integ.umap.sample.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()

p <- DimPlot(integ_atac, group.by = 'seurat_clusters', reduction = "atac.umap", pt.size = 0.1)
pdf(paste0(outdir, '/snatac_integ.umap.cluster.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()


# remove processed
snatac_peak_list <- NULL


# Combind snATAC and snRNA
print('[INFO] integrate multiome ...')
integrated_multiome <- RunWNN_Customized(integ_rna=integ_rna, integ_atac=integ_atac, 
                            reduc_rna='harmony_rna', reduc_atac='integrated_lsi', dim_max_rna=30, dim_max_atac=30, 
                            weight_name = "RNA.weight", res_clus=0.4, seed=seed)


# save integrated data
saveRDS(integrated_multiome, file=paste0(outdir,'/multiome_integrated.harmony_lsi.rds'))
write.table(integrated_multiome@meta.data, file = paste0(outdir, '/multiome_integrated.harmony_lsi.metadata.xls'), sep = "\t", quote=F, col.names = NA)

p <- DimPlot(integrated_multiome, reduction = 'wnn.umap', group.by = 'Sample', pt.size = 0.1)
pdf(paste0(outdir, '/multiome_integrated.harmony_lsi.wnnUMAP.sample.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()

p <- DimPlot(integrated_multiome, reduction = 'wnn.umap', group.by = 'seurat_clusters', pt.size = 0.1)
pdf(paste0(outdir, '/multiome_integrated.harmony_lsi.wnnUMAP.cluster.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()



######################
##  add gene activity and motifs 
######################

print('[INFO] call pfm ...')
if (nchar(fpfm) == 0){
    pfm <- CallPFMFromJASPAR(JASPAR2022)
    saveRDS(pfm, file = paste0(outdir, "/pfm.jaspar2022_vertebrates_core.rds"))
} else {
    pfm <- readRDS(fpfm)
}

integrated_multiome <- AddGeneActivityAndMotif(obj=integrated_multiome, genome=genome, pfm=pfm, assay_use='peaks')
#saveRDS(integrated_multiome, file=paste0(outdir,'/multiome_integrated.harmony_lsi.plus.rds'))


# Seurat v5
integrated_multiome[['RNA']] <- JoinLayers(integrated_multiome[['RNA']])
saveRDS(integrated_multiome, file=paste0(outdir,'/multiome_integrated.harmony_lsi.plus2.rds'))








