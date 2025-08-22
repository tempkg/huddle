library(seuextra)
library(ggplot2)
library(Seurat)
library(data.table)

rds <- 'scrna_sct2_harmony.Lineage-1.rds'
outpdf <- 'iue_zr.exp.pdf'

seu <- readRDS(rds)
info <- fread('out_celltype/cluster_cellType.xls', data.table=F)
seu$cell_type2 <- info$cell_type2[match(seu$seurat_clusters, info$seurat_clusters)]

DefaultAssay(seu) <- 'SCT'
Idents(seu) <- 'cell_type2'

levels(seu) <- c('CycProg', 'Astrocyte', 'Neuron', 'Unclassified')

features <- c('Htr2a','Htr2b','Gphn','Homer1','Mpp2','Shank3')


p <- Seurat_DotPlot(obj=seu, features=features, assay='SCT', col_set="exp")
ggsave(outpdf, w=3, h=3)






