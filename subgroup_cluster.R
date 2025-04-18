
library(Seurat)
library(ggplot2)

rds <- '' 
clu <- c(8,4,1,0,10,6,7,5,11,2)
outdir <- 'out_subgroup'


dir.create(outdir)


seu <- readRDS(rds)
print(dim(seu))
subseu <- subset(x = seu, idents = seurat_clusters %in% clu)
print(dim(subseu))
print(unique(subseu$seurat_clusters))

subseu <- RunUMAP(subseu, dims = 1:30, reduction = 'integrated.rna', reduction.name = 'rna.umap')
subseu <- FindNeighbors(subseu, reduction = 'integrated.rna', dims = 1:30)
subseu <- FindClusters(subseu, resolution = 0.5)

# plot umap
p <- DimPlot(subseu, group.by = 'Sample', reduction = "rna.umap", pt.size = 0.1)
ggsave(paste0(outdir, '/subgroup.umap.pdf'), width = 5.5, height = 4)




