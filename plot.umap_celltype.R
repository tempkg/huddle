
library(GetoptLong)
library(Seurat)
library(seuextra)
library(ggplot2)
library(data.table)

#source('~/Documents/_Toolkit/_gitpkg/doc.colorset/color_celltype/color.celltype.R')



rds <- 'seurat5.1_v6.2/multiome_integrated.plus.rds'
umap <- 'wnn.umap'
fmeta <- 'out_celltype/cluster_cellType.xls'
outdir <- 'out_celltype'

GetoptLong(
    "rds=s",    ".rds seurat file",
    "umap=s",   "umap",
    "fmeta=s",  "fmeta",
    "outdir=s",  "outdir"
)


seu <- readRDS(rds)
meta <- fread(fmeta, data.table=F)
seu$cell_type <- meta$cell_type[match(seu$seurat_clusters, meta$seurat_clusters)]
seu$cell_type2 <- meta$cell_type2[match(seu$seurat_clusters, meta$seurat_clusters)]

# seurat plot
#tumor_cluster <- seu$seurat_clusters[!seu$seurat_clusters %in% c(28,29,22,26,2)]
#seu$cell_type2[seu$seurat_clusters %in% tumor_cluster] <- 'Tumor'


# scCustom plot	
#p <- scCustom_DimPlot(obj=seu, reduction=umap, group_by='cell_type2', colors=colorS2, label=T)
p <- scCustom_DimPlot(obj=seu, reduction=umap, group_by='cell_type2', label=T)
ggsave(paste0(outdir, '/umap.celltype.corrected.pdf', w=8, h=6)




