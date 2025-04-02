# Hua Sun
# 4/2/25 v0.2
# seurat v5.1.0

library(Seurat)
library(data.table)
library(seuproc)
library(dplyr)
library(GetoptLong)
library(ggplot2)

set.seed(42)
options(future.globals.maxSize=10000000000)



h5 <- 'filtered_feature_bc_matrix.h5'
name <- 'sample'
ref <- 'mm10'
min_cells <- 3
outdir <- 'out_seurat5.single'


GetoptLong(
    "name=s",      "sample name",
    "h5=s",        "h5 file",
    "ref=s",       "ref",
    "group=s",     "group",
    "min_cells=i", "min_cells",
    "barcode_metrics=s", "barcode_metrics",
    "outdir=s",    "outdir"
)


dir.create(outdir)

print(name)


seu <- CreateSeuratObjectFromH5(h5=h5, min_cells=min_cells, name=name, ref=ref)
print(dim(seu))


# QC plot
p <- QCPlotRNA(seu)
ggsave(paste0(outdir, '/', name, '.rawdata.qc.pdf'), w=5, h=2.5, useDingbats=T)


# filter
cutoff <- CalMaxCutoffRNA(obj=seu, maxratioCount=0.99, maxratioFeature=0.9, revise_val=TRUE)

print(cutoff)
seu <- FilterRNACells(obj=seu, ncount_rna_min=2000, nfeature_rna_min=500, ncount_rna_max=cutoff[1], nfeature_rna_max=cutoff[2], percent_mt=10)

seu <- RunDoubletFinder(seu, outdir)
saveRDS(seu, file=paste0(outdir, '/filtered.remDoublet.rds'))


# normalize
# vst.flavor = 'v2'
seu <- Seurat5SCTNormalize(obj=seu, regress="percent.mt", max_dim=30, res_clus=0.6)
saveRDS(seu, file=paste0(outdir, '/scrna.remDuplet.sct2.rds'))


