# Hua Sun
# 2026-02-02 v3


library(Seurat)
library(Signac)
library(seuproc)
library(scqc)
library(dplyr)
library(stringr)
library(stringi)
library(ggplot2)
library(yaml)
library(GetoptLong)


dir <- 'sample/outs'
name <- 'sample'
ref <- 'hg38'  # mm10,hg38
f <- 'filter.yml'
outdir <- 'out_signac'

GetoptLong(
    "dir=s",     "dir outs",
    "name=s",    "sample name",
    "ref=s",     "ref",
    "f=s",         "yaml",
    "outdir=s",  "outdir"
)

dir.create(outdir)


h5 <- paste0(dir, '/filtered_peak_bc_matrix.h5')
fragment <- paste0(dir, '/fragments.tsv.gz')
fsc <- paste0(dir, '/singlecell.csv')

# set
dyaml <- yaml.load_file(f)

min_cells <- dyaml$min_cells
ncount_atac_min <- dyaml$ncount_atac_min
nucle_signal <- dyaml$nucleosome_signal
tss_enrich_min <- dyaml$tss_enrich_min
tss_enrich_max <- dyaml$tss_enrich_max
bl_ratio <- dyaml$blacklist_ratio
pct_in_peaks <- dyaml$pct_reads_in_peaks
qt <- dyaml$qt
seed <- dyaml$seed
max_dim <- dyaml$max_dim
clus_resolution <- dyaml$clus_resol 

macs2 <- dyaml$macs2
annotation <- Annotations(ref)



set.seed(seed)

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
ggsave(paste0(outdir, '/', name, '.rawdata.qc.pdf'), w=6, h=3, useDingbats=T)



# Filter
signac_obj <- High_quality_filter_ATAC(
    obj = signac_obj,
    min_ncount=ncount_atac_min,
    max_nucleo=nucle_signal,
    min_tss=tss_enrich_min,
    max_tss=tss_enrich_max,
    bl_ratio=bl_ratio,
    pct_in_peaks=pct_in_peaks,
    qt_cutoff=qt
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
# FindTopFeatures(signac_obj, min.cutoff = 10)
signac_obj <- ComputeLSI(obj=signac_obj, assay='peaks', min_cutoff=10)
signac_obj <- RunUMAP_Plus(obj=signac_obj, reduc='lsi', min_dim=2, max_dim=max_dim, algorithm_clus=3, res_clus=clus_resolution, seed=seed)

# out filtered
saveRDS(signac_obj, file=paste0(outdir, '/4.snatac.filtered.normalized.rds'))
write.table(signac_obj@meta.data, file = paste0(outdir, "/snatac.filtered.normalized.metaData.xls"), sep="\t", quote=F, col.names = NA)


# out umap
pdf(paste0(outdir, '/umap.snatac.filtered.normalized.pdf'), w=6, h=6)
p <- DimPlot(signac_obj, reduction='umap', pt.size=0.1)
print(p)
dev.off()



# add gene activity
DefaultAssay(signac_obj) <- 'peaks'
gene.activity <- Signac::GeneActivity(signac_obj)
saveRDS(gene.activity, file=paste0(outdir,'/gene_activity.rds'))

signac_obj <- AddGeneActivity(obj=signac_obj, gene_activity=gene.activity, assay='peaks', new_assay='RNA')
saveRDS(signac_obj, file=paste0(outdir,'/5.snatac.filtered.normalized.with_geneact.rds'))













