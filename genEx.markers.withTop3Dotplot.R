
library(Seurat)
library(seuextra)
library(stringr)
library(dplyr)
library(GetoptLong)
library(RColorBrewer)
library(ggplot2)
library(data.table)


rds <- 'ZR.multiome_integrate.30_30.rds'
assay <- 'SCT'
idents <- 'cell_type2'
finfo <- ''
outdir <- 'out_expMarkers.zr'

GetoptLong(
    "rds=s",       ".rds seurat file",
    "assay=s",     "assay",
    "idents=s",    "idents",
    "finfo=s",     "info",
    "outdir=s",    "output path"
)

dir.create(outdir)


seu <- readRDS(rds)
DefaultAssay(seu) <- assay


d_markers <- fread(finfo, data.table=F)
# known genes
genes <- c('GFAP','AQP4','APOE', 'MKI67','TOP2A','CDK1','DNAH12','CFAP157','FOXJ1','COL3A1','COL1A1','COL1A2','PID1','SLITRK5','GRM8','VIM','NES','HES1')

#
d_markers <- d_markers %>% filter(gene %in% genes)

dx_log2fc <- d_markers %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = T)
p <- Seurat_DotPlotMarkers(obj=seu, group=idents, d_markers=dx_log2fc)
ggsave(paste0(outdir, '/top3_markers.logfc.pdf'), w=5, h=3)








