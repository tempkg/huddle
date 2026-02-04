# Hua Sun
# 2026-1-20 v2.5


library(Seurat)
library(Signac)
library(data.table)
library(stringr)
library(stringi)
library(seuproc)
library(dplyr)
library(TFBSTools)
library(GetoptLong)
library(ggplot2)

options(future.globals.maxSize=100000000000)


ids <- 'sample_list.txt'
ref <- 'hg38'  # mm10,hg38

sdir <- 'out_signac.v2'
frds <- '2.snatac.filtered.rmDoublet.rds'


method <- 'harmony'  # 'rlsi'
max_dim <- 30
seed <- 42

fpfm <- 'pfm.jaspar2024.human_core.all.rds'

outdir <- 'out_snatac_harmony'

GetoptLong(
    "ids=s",      "id list",
    "ref=s",       "ref",
    "sdir=s",     "seurat result dir.",
    "frds=s",     "read file",
    "method=s",    "integration method",
    "seed=i",      "seed",
    "max_dim=i",   "max_dim",
    "extra",       "extra",
    "fpfm=s",       "fpfm",
    "outdir=s",    "outdir"
)

set.seed(seed)

dir.create(outdir)

macs2 <- '/research/groups/mackgrp/home/common/Software/miniconda3/envs/macs2/bin/macs2'
annotation <- Annotations(ref)

# read
sample_list <- readLines(ids)


seurat_list <- NULL

for (name in sample_list){
        print(name)

        fseu <- paste0(sdir, '/', name, '/', frds)
        temp_seu <- readRDS(fseu)

        # make data list
        seurat_list <- append(seurat_list, temp_seu)
        
        temp_seu <- NULL
}

# save removed doublet list
saveRDS(seurat_list, file=paste0(outdir,'/snatac.remDoublet.list.rds'))



# snATAC: 1.Call peaks using MACS2
print('[INFO] macs2 ...')
snatac_macs2_list <- lapply(X = seurat_list, FUN = function(x) { x <- CallPeaksUsingMACS2(obj=x, macs2=macs2, new_assay='peaks', annotation=annotation) })
saveRDS(snatac_macs2_list, file=paste0(outdir, '/snatac_macs2.list.rds'))


# snATAC: 2.Combine peaks & recall fragment
print('[INFO] combine peaks ...')
snatac_peak_list <- ObjectListConvertBedToGRanges(snatac_macs2_list, annotation)
snatac_peak_list <- lapply(X = snatac_peak_list, FUN = function(x) { x <- ComputeLSI(x, 'peaks', 10) })
saveRDS(snatac_peak_list, file=paste0(outdir, '/snatac_peaks.list.rds'))

# remove processed
snatac_macs2_list <- NULL


# snATAC: 3.Integration 
if (method == 'harmony'){
    integ_atac <- ATACIntegrationHarmony (snatac_list=snatac_peak_list, group='Sample', new_reduc='harmony_atac', max_dim=max_dim, res_clus=0.4, seed=seed)
} else {
    integ_atac <- ATACIntegration(snatac_list=snatac_peak_list, group='Sample', new_reduc='integrated_lsi', max_dim=max_dim, res_clus=0.4, seed=seed)
}
saveRDS(integ_atac, file=paste0(outdir,'/snatac_integrated_',method,'.withumap.rds'))
write.table(integ_atac@meta.data, file = paste0(outdir,'/snatac_integrated_',method,'.metadata.xls'), sep = "\t", quote=F, col.names = NA)

p <- DimPlot(integ_atac, group.by = 'Sample', reduction = "atac.umap", pt.size = 0.1)
pdf(paste0(outdir,'/snatac_integrated_',method,'.umap.pdf'), width = 6.5, height = 5, useDingbats=FALSE)
print(p)
dev.off()


# remove processed
snatac_peak_list <- NULL




##---------------- Extra ----------------##
if (extra){
    pfm <- readRDS(fpfm)
    genome <- UCSCBSGenome(ref)
    integ_atac <- AddMotifActivity(obj=integ_atac, genome=genome, pfm=pfm, assay='peaks')
    saveRDS(integ_atac, file=paste0(outdir,'/snatac_integrated_',method,'.with_motif.rds'))

    DefaultAssay(integ_atac) <- 'peaks'
    gene.activity <- Signac::GeneActivity(integ_atac)
    saveRDS(gene.activity, file=paste0(outdir,'/snatac_integrated.gene_activity.rds'))

    integ_atac <- AddGeneActivity(obj=integ_atac, gene_activity=gene.activity, assay='peaks', new_assay='RNA')
    saveRDS(integ_atac, file=paste0(outdir,'/snatac_integrated_',method,'.with_motif_geneact.rds'))
}





