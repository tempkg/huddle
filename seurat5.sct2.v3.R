# Hua Sun
# 2026-1-26 v3
# seurat v5.4.0

library(Seurat)
library(data.table)
library(seuproc)
library(scqc)
library(dplyr)
library(yaml)
library(ggplot2)
library(GetoptLong)

options(future.globals.maxSize=100000000000)


dir <- 'outs'
name <- 'sample'
ref <- 'human' # mouse
f <- 'filter.yml'
outdir <- 'out_seurat5'


GetoptLong(
    "name=s",     "sample name",
    "dir=s",      "outs cellranger",
    "ref=s",      "ref",
    "f=s",        "yaml",
    "outdir=s",   "outdir"
)

dir.create(outdir)
print(name)

dyaml <- yaml.load_file(f)

min_cells <- dyaml$min_cells
percent_mt <- dyaml$percent_mt
ncount_rna_min <- dyaml$ncount_rna_min
nfeature_rna_min <- dyaml$nfeature_rna_min
max_nfeature <- dyaml$nfeature_rna_max
qt <- dyaml$qt

doublet_rate <- dyaml$doublet_rate

seed <- dyaml$seed
max_dim <- dyaml$max_dim

clus_resolution <- dyaml$clus_resol 

set.seed(seed)

h5 <- paste0(dir, '/filtered_feature_bc_matrix.h5')
seu <- CreateSeuratObjectFromH5(h5=fh5, min_cells=min_cells, name=name, ref=ref)
print(dim(seu))


# QC plot
p <- QCPlotRNA(seu)
ggsave(paste0(outdir, '/', name, '.rawdata.qc.pdf'), w=5, h=2.5, useDingbats=T)


# filter v2026-01-27
seu <- High_quality_filter(obj=seu, 
        min_ncount=ncount_rna_min, 
        min_nfeature=nfeature_rna_min, 
        percent_mt=percent_mt, 
        max_nfeature=max_nfeature,
        qt_cutoff=qt)
print(nrow(seu@meta.data))

seu <- RunDoubletFinder2(obj=seu, rate=doublet_rate, max_dim=10, singlet=TRUE, outdir=outdir)
print(nrow(seu@meta.data))
seu <- RunFilterEmptyDrops(dir=dir, obj=seu, outdir=outdir)
print(nrow(seu@meta.data))
saveRDS(seu, file=paste0(outdir, '/filtered.remDoublet_emptyDrops.rds'))


# normalize
# vst.flavor = 'v2'
seu <- Seurat5SCTNormalize(obj=seu, regress="percent.mt", variable=3000)
seu <- RunUMAP(seu, dims = 1:max_dim, reduction = "pca", seed.use=seed)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:max_dim)
seu <- FindClusters(seu, resolution = clus_resolution)
saveRDS(seu, file=paste0(outdir, '/scrna.remDoublet_emptyDrops.sct2.rds'))

# out metadata
write.table(seu@meta.data, file=paste0(outdir, '/scrna.remDoublet_emptyDrops.sct2.metadata.xls'), sep='\t', quote=F, col.names=NA)


# umap
p <- DimPlot(seu)
ggsave(paste0(outdir, '/scrna.remDoublet_emptyDrops.sct2.umap.pdf'), w=6, h=5)




