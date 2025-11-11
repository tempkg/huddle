# Hua Sun
# 2025-11-08 v2.3
# seurat v5.1.0

library(Seurat)
library(data.table)
library(seuproc)
library(scqc)
library(dplyr)
library(GetoptLong)
library(ggplot2)


options(future.globals.maxSize=100000000000)


dir <- 'outs'
name <- 'sample'
ref <- 'mouse' # human
min_cells <- 3
percent_mt <- 10
seed <- 42
outdir <- 'out_seurat5'


GetoptLong(
    "name=s",      "sample name",
    "dir=s",       "outs cellranger",
    "ref=s",       "ref",
    "min_cells=i", "min_cells",
    "percent_mt=i", "percent_mt",
    "seed=i",      "seed",
    "outdir=s",    "outdir"
)

set.seed(seed)
dir.create(outdir)

print(name)


h5 <- paste0(dir, '/filtered_feature_bc_matrix.h5')
seu <- CreateSeuratObjectFromH5(h5=h5, min_cells=min_cells, name=name, ref=ref)
print(dim(seu))


# QC plot
p <- QCPlotRNA(seu)
ggsave(paste0(outdir, '/', name, '.rawdata.qc.pdf'), w=5, h=2.5, useDingbats=T)


# filter
cutoff <- CalMaxCutoffRNA(obj=seu, maxratioCount=0.99, maxratioFeature=0.9, revise_val=TRUE)

print(cutoff)
seu <- subset(x = seu,
            subset = nCount_RNA > 2000 &
                nFeature_RNA > 500 &
                nCount_RNA < cutoff[1] &
                nFeature_RNA < cutoff[2] &
                percent.mt < percent_mt
        )
print(nrow(seu@meta.data))

seu <- RunDoubletFinder(obj=seu, outdir=outdir, singlet=TRUE)
print(nrow(seu@meta.data))
seu <- RunFilterEmptyDrops(dir=dir, obj=seu, outdir=outdir)
print(nrow(seu@meta.data))
saveRDS(seu, file=paste0(outdir, '/filtered.remDoublet_emptyDrops.rds'))


# normalize
# vst.flavor = 'v2'
seu <- Seurat5SCTNormalize(obj=seu, regress="percent.mt", variable=3000)
seu <- RunUMAP(seu, dims = 1:30, reduction = "pca", seed.use=seed)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30)
seu <- FindClusters(seu, resolution = 0.4)
saveRDS(seu, file=paste0(outdir, '/scrna.remDoublet_emptyDrops.sct2.rds'))

# out metadata
write.table(seu@meta.data, file=paste0(outdir, '/scrna.remDoublet_emptyDrops.sct2.metadata.xls'), sep='\t', quote=F, col.names=NA)


# umap
p <- DimPlot(seu)
ggsave(paste0(outdir, '/scrna.remDoublet_emptyDrops.sct2.umap.pdf'), w=6, h=5)










