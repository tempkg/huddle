# Hua Sun
# 2025-11-08 v2


library(Seurat)
library(Signac)
library(seuproc)
library(scqc)
library(GenomicRanges)
library(future)

library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(limma)
library(biovizBase)

library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SeuratDisk)

library(chromVAR)
library(JASPAR2024)
library(TFBSTools)
library(motifmatchr)

library(stringr)
library(stringi)

library(patchwork)
library(data.table)


library(GetoptLong)


dir <- 'sample/outs'
name <- 'sample'
ref <- 'hs38'  # mm10,hg38
min_cells <- 10
macs2 <- 'macs2'
seed <- 42
outdir <- 'out_signac'

GetoptLong(
    "dir=s",       "dir outs",
    "name=s",      "sample name",
    "ref=s",       "ref",
    "min_cells=i", "min_cells",
    "macs2=s",     "MACS2 path",
    "seed=i",      "seed",
    "extra",       "Extra process",
    "outdir=s",    "outdir"
)

set.seed(seed)

dir.create(outdir)

h5 <- paste0(dir, '/filtered_peak_bc_matrix.h5')
fragment <- paste0(dir, '/fragments.tsv.gz')
fsc <- paste0(dir, '/singlecell.csv')

# set
macs2 <- '/research/groups/mackgrp/home/common/Software/miniconda3/envs/macs2/bin/macs2'
annotation <- Annotations(ref)
genome <- UCSCBSGenome(ref)

print(name)
print(ref)


# Create object
signac_obj <- CreateSignacObjectFromH5(h5=h5, fsc=fsc, min_cells=min_cells, assay='ATAC', fragment=fragment, ref=ref, annotation=annotation)
print('[INFO] Created object')
print(dim(signac_obj))

signac_obj$Sample <- name
signac_obj$ref <- ref


# QC plot
p <- QCPlotATAC(signac_obj)
ggsave(paste0(outdir, '/', name, '.rawdata.qc.pdf'), w=6, h=1.5, useDingbats=T)



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
saveRDS(signac_obj, file=paste0(outdir, '/1.snatac.filtered.rds'))


#------// doublet filter
blacklist <- blacklist_hg38_unified
if (ref == 'mm10'){ blacklist <- blacklist_mm10 }

df <- Run_scDblFinder_ATAC(obj=signac_obj, fragfile=fragment, repeats=blacklist)
write.table(df, paste0(outdir, '/out_scDblFinder_pval.tsv'), sep='\t', quote=F, col.names=NA)

cells_singl <- rownames(df[!df$is_doublet,])
#writeLines(cells_singl, 'filtered_doublet.cell')
signac_obj <- subset(signac_obj, cells=cells_singl)

saveRDS(signac_obj, file=paste0(outdir, '/2.snatac.filtered.rmDoublet.rds'))
#------\\


# Re-Call peaks using MACS2 
# ref in seurat object
signac_obj <- CallPeaksUsingMACS2(obj=signac_obj, macs2=macs2, new_assay='peaks', annotation = annotation)
saveRDS(signac_obj, file=paste0(outdir, '/3.snatac.filtered.macs2.rds'))


# Normalization and UMAP
signac_obj <- ComputeLSI(obj=signac_obj, assay='peaks', min_cutoff=10)
signac_obj <- RunUMAP_Plus(obj=signac_obj, reduc='lsi', min_dim=2, max_dim=30, algorithm_clus=3, res_clus=0.4, seed=seed)

# out filtered
saveRDS(signac_obj, file=paste0(outdir, '/4.snatac.filtered.normalized.rds'))
write.table(signac_obj@meta.data, file = paste0(outdir, "/snatac.filtered.normalized.metaData.xls"), sep="\t", quote=F, col.names = NA)


# out umap
pdf(paste0(outdir, '/umap.snatac.filtered.normalized.pdf'), w=6, h=6)
p <- DimPlot(signac_obj, reduction='umap', pt.size=0.1)
print(p)
dev.off()





##---------------- Extra ----------------##
if (extra){
    # Add gene activity
    signac_obj <- AddGeneActivity(obj=signac_obj, assay='peaks', new_assay='RNA')


    # Add motif activity
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

    signac_obj <- AddMotifActivity(obj=signac_obj, genome=genome, pfm=pfm, assay='peaks')

    # output
    saveRDS(signac_obj, file=paste0(outdir, '/5.snatac.filtered.normalized.with_act_motif.rds'))
    #write.table(signac_obj@meta.data, file = paste0(outdir, "/snatac.filtered.normalized.with_act_motif.metaData.xls"), sep="\t", quote=F, col.names = NA)
}










