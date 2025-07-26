
library(Seurat)
library(scCustomize)

rds <- 'seurat.rds' 
finfo <- 'merged.clean_lineage_cells.tsv'

seu <- readRDS(rds)

info <- fread(finfo, data.table=F)
seu$lineage_group <- info$group[match(rownames(seu@meta.data), info$cell)]

# https://samuel-marsh.github.io/scCustomize/articles/Gene_Expression_Plotting.html
Idents(seu) <- 'lineage_group'
Cluster_Highlight_Plot(seurat_object = seu, cluster_name = c("Non-Lineage", "Lineage-1"), highlight_color = c("navy", "red"))



