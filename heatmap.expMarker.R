# Hua Sun



library(Seurat)
library(data.table)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(scapekit)
library(RColorBrewer)
library(randomcoloR)
library(circlize)
library(GetoptLong)



tf <- '~/OneDrive/mydoc/database/cistarget/tf_lists/allTFs_hg38.txt'
group_by <- 'group'
width <- 3.5
height <- 1.5

rds <- '2.multiome_integrated_plus.regout_mt.malignantCellOnly.rds'
finfo <- 'sample.group.txt'
fmarker <- '5.out_genExpGroupMarkers/group.findallmarker.geneExp.xls'
outdir <- '5.out_genExpGroupMarkers'


dir.create(outdir)


# 1.read object
print('[INFO] reading .rds ...')
seurat_obj <- readRDS(rds)
info <- read.table(finfo, sep='\t', header=F)
colnames(info) <- c('sample', 'group')
seurat_obj$group <- info$group[match(seurat_obj$Sample, info$sample)]


d_exp_markers <- read.table(fmarker, sep='\t', header=T)
# "p_val"  "avg_log2FC" "pct.1"  "pct.2"  "p_val_adj"  "cluster"  "gene"  "diff_pct"  

# read tf list
tf_list <- readLines(tf)
# add TF to marker
d_exp_markers$TF <- tf_list[match(d_exp_markers$gene, tf_list)]


d_exp_marker_filter <- d_exp_markers %>% 
            filter(p_val_adj < 0.01) %>%
            filter(avg_log2FC > 1) %>% 
            filter(pct.1 > 0.35) %>% 
            filter(pct.2 < 0.15)
#write.table(d_exp_marker_filter, paste0(outdir, '/expMarker.filter.xls'), sep='\t', quote=F, row.names=F)


Idents(seurat_obj) <- group_by
seurat_obj_sub <- subset(x=seurat_obj, downsample=500)
metadata <- seurat_obj_sub@meta.data

d_matrix <- as.matrix(seurat_obj_sub$SCT@data)

show_gene <- intersect(gene_info$gene, tf_list)

# best for z-score
outfile <- paste0(outdir,'/heatmap.expMarker.filtered.1.pdf')
ht <- Heatmap_DiffMarkers_Vertical(
    data=d_matrix, 
    meta=metadata, 
    diff_marker=d_exp_marker_filter,
    group='group', 
    sort_group=c('ZR', 'EWP', 'PF'),
    col_group=c("ZR" = "#E21F26", "EWP" = "#F57F20", "PF" = "#7C50A0"),
    show_gene=show_gene
)

pdf(paste0(outdir,'/heatmap.expMarker.filtered.1.pdf'), width=2.2, height=3, useDingbats=FALSE)
draw(ht)
dev.off()


