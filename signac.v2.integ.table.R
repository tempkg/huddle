# Hua Sun
# Seurat >=5.1
# 2025-11-08 v2


library(Seurat)
library(Signac)
library(harmony)
library(data.table)
library(seuproc)
library(scqc)
library(dplyr)
library(JASPAR2024)
library(GetoptLong)
library(ggplot2)


options(future.globals.maxSize=100000000000)



table <- 'data.table'
ref <- 'hs38'  # mm10,hg38
min_cells <- 10
macs2 <- 'macs2'
method <- 'harmony'  # 'rlsi'
max_dim <- 30
seed <- 42

outdir <- 'out_snatac_harmony'

GetoptLong(
    "table=s",     "table",
    "ref=s",       "ref",
    "min_cells=i", "min_cells",
    "method=s",    "integration method",
    "macs2=s",     "MACS2 path",
    "seed=i",      "seed",
    "max_dim=i",   "max_dim",
    "extra",       "extra",
    "outdir=s",    "outdir"
)

set.seed(seed)
dir.create(outdir)

macs2 <- '/research/groups/mackgrp/home/common/Software/miniconda3/envs/macs2/bin/macs2'
annotation <- Annotations(ref)
genome <- UCSCBSGenome(ref)



outdir_process <- paste0(outdir, '/processed')
dir.create(outdir_process)

seurat_list <- NULL

dt <- fread(table, header=FALSE, data.table=F)



for (name in dt$V1){
        print(name)
        
        sample_dir <- dt$V2[dt$V1==name]
        h5 <- paste0(sample_dir, '/filtered_peak_bc_matrix.h5')
        fragment <- paste0(sample_dir, '/fragments.tsv.gz')
        fsc <- paste0(sample_dir, '/singlecell.csv')

        # Create object
        signac_obj <- CreateSignacObjectFromH5(h5=h5, fsc=fsc, min_cells=min_cells, assay='ATAC', fragment=fragment, ref=ref, annotation=annotation)
        print('[INFO] Created object')
        print(dim(signac_obj))

        signac_obj$Sample <- name
        signac_obj$ref <- ref


        # QC plot
        p <- QCPlotATAC(signac_obj)
        ggsave(paste0(outdir_process, '/', name, '.rawdata.qc.pdf'), w=6, h=1.5, useDingbats=T)


        # Filter
        # cell-level filtering
        max_ncount_atac <- Max_nCountX(signac_obj@meta.data$nCount_ATAC, 0.9)
        signac_obj <- subset(
            x = signac_obj,
            subset = nCount_ATAC > 3000 &
                nCount_ATAC < max_ncount_atac &
                nucleosome_signal < 4 &
                TSS.enrichment > 4 &
                TSS.enrichment < 15 &
                blacklist_ratio < 0.01 &
                pct_reads_in_peaks > 40
        )
        print(dim(signac_obj))

        # out filtered
        #saveRDS(signac_obj, file=paste0(outdir, '/1.snatac.filtered.rds'))


        #------// doublet filter
        blacklist <- blacklist_hg38_unified
        if (ref == 'mm10'){ blacklist <- blacklist_mm10 }

        df <- Run_scDblFinder_ATAC(obj=signac_obj, fragfile=fragment, repeats=blacklist)
        write.table(df, paste0(outdir, '/', name, '.out_scDblFinder_pval.tsv'), sep='\t', quote=F, col.names=NA)

        cells_singl <- rownames(df[!df$is_doublet,])
        #writeLines(cells_singl, 'filtered_doublet.cell')
        signac_obj <- subset(signac_obj, cells=cells_singl)

        saveRDS(signac_obj, file=paste0(outdir_process, '/', name, '.snatac.filtered.rmDoublet.rds'))
        #------\\


        # make data list
        seurat_list <- append(seurat_list, signac_obj)
}


# snATAC: 1.Call peaks using MACS2
print('[INFO] macs2 ...')
snatac_macs2_list <- lapply(X = multiome_filtered_list, FUN = function(x) { x <- CallPeaksUsingMACS2(obj=x, macs2=macs2, new_assay='peaks', annotation=annotation) })
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
    integ_atac <- ATACIntegrationHarmony (snatac_list=snatac_peak_list, group='Sample',  new_reduc='harmony_atac', max_dim=max_dim, res_clus=0.4, seed=seed)
} else {
    integ_atac <- ATACIntegration(snatac_list=snatac_peak_list, group='Sample',  new_reduc='integrated_lsi', max_dim=max_dim, res_clus=0.4, seed=seed)
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
    print('[INFO] call pfm ...')
    if (nchar(fpfm) == 0){
        if (ref == 'hs38'){
            # human
            pfm <- getMatrixSet(x = JASPAR2024, opts = list(species = 9606, all_versions = FALSE))
        } else {
            #pfm <- CallPFMFromJASPAR(JASPAR2024)
            pfm <- getMatrixSet(x = JASPAR2024, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
        }
        saveRDS(pfm, file = paste0(outdir, "/pfm.jaspar2024_vertebrates_core.rds"))
    } else {
        pfm <- readRDS(fpfm)
    }

    integ_atac <- AddMotifActivity(obj=integ_atac, genome=genome, pfm=pfm, assay='peaks')
    integ_atac <- AddGeneActivity(obj=integ_atac, assay='peaks', new_assay='RNA')
    
    saveRDS(integ_atac, file=paste0(outdir,'/snatac_integrated_lsi.plus.rds'))
}





