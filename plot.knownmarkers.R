# 2/1/25

library(Seurat)
library(seuextra)
library(ggplot2)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(GetoptLong)
library(data.table)
library(tidyr)
library(dplyr)


rds <- 'ZR.multiome_integrate.30_20.rds'
fdb <- 'out_expMarkers.cell_type2/cell_type2.findallmarker.geneExp.xls'
fmeta <- 'final.cluster_cellType.corrected.xls'

markers <- c('GFAP','AQP4','APOE',...) 

width <- 2.5
height <- 2.1
outdir <- 'out_checkKnownMarkers'

GetoptLong(
    "rds=s",       ".rds seurat file",
    "fdb=s",      "cluster markers",
    "fmeta=s", "known marker gene",
    "width=f",     "width",
    "height=f",    "height",
    "outdir=s",    "output path"
)


dir.create(outdir)

# read object
seu <- readRDS(rds)
meta <- fread(fmeta, data.table=F)
seu$cell_type2 <- meta$cell_type2[match(seu$seurat_clusters, meta$seurat_clusters)]

diff <- fread(fdb, data.table=F)

diff <- diff %>% 
        filter(avg_log2FC > 0.585 & diff_pct > 0.05 & pct.1 > 0.2) %>%
        filter(gene %in% markers)

write.table(diff, paste0(outdir, "/knownmarkers.cell_type2.xls"), sep='\t', quote=F, row.names=F)


# re-sort - check both of them & choice the best one
dx_log2fc <- diff %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = T)
p <- Seurat_DotPlotMarkers(obj=seu, group='cell_type2', d_markers=dx_log2fc)
ggsave(paste0(outdir, '/topX_markers.resorted_log2fc.pdf'), w=width, h=height)




