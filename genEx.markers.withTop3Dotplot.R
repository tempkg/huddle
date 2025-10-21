
library(Seurat)
library(seuextra)
library(stringr)
library(dplyr)
library(GetoptLong)
library(RColorBrewer)
library(ggplot2)
library(data.table)


# known genes
genes <- c('GFAP','AQP4','APOE', 'MKI67','TOP2A','CDK1','DNAH12','CFAP157','FOXJ1','COL3A1','COL1A1','COL1A2','PID1','SLITRK5','GRM8','VIM','NES','HES1')


rds <- 'ZR.multiome_integrate.30_30.rds'
fct <- 'cluster_cellType.corrected.xls'
fmarker <- 'meta.data'

GetoptLong(
    "rds=s",       ".rds seurat file",
    "fct=s",       "cell type",
    "fmarker=s",   "marker data"
)

dir.create(outdir)


seu <- readRDS(rds)
DefaultAssay(seu) <- 'SCT'

d_ct <- fread(fct, data.table=F)
seu$cell_type2 <- d_ct$cell_type2[match(seu$seurat_clusters, d_ct$seurat_clusters)]

d_markers <- fread(fmarker, data.table=F)
d_markers <- d_markers %>% filter(gene %in% genes)

dx_log2fc <- d_markers %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = T)
p <- Seurat_DotPlotMarkers(obj=seu, group='cell_type2', d_markers=dx_log2fc)
ggsave('top3_markers.logfc.pdf', w=5, h=3)








