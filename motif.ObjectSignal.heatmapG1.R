library(scapekit)
library(ComplexHeatmap)
library(data.table)




# extract markers
d_motif_markers <- fread('motifMarkers.celltypes.remDup.top100.xls', data.table=F)
d_motif_markers <- d_motif_markers[, c('cluster', 'gene', 'TF')]
colnames(d_motif_markers) <- c('cell_type2', 'motif', 'gene')
rownames(d_motif_markers) <- d_motif_markers$motif
d_motif_markers$motif <- NULL


# extract motif data from Seurat object
seu <- readRDS('../seurat5.1_v6.2/multiome_integrated.plus.rds')
fmeta <- '../out_celltype/metaData.cellType.xls'


Idents(seu) <- 'cell_type2'
sub_seu <- subset(seu, downsample=500)
meta <- sub_seu@meta.data

d_matrix <- as.matrix(sub_seu@assays$chromvar@data)



colorG1 <- c(
        "RGC" = "#BC243C",
        "CycProg" = "#EFC050",
        "Astrocyte" = "#EF563C",
        "Epend" = "#4eaf49",
        "Neuron" = "#0F4C81"
        )



# best for z-score
outfile <- paste0(outdir,'/heatmap.motifMarker.pdf')
ht <- Heatmap_DiffMarkers_Vert(
    data=d_matrix, 
    meta=metadata, 
    diff_marker=d_motif_markers,
    group='group', 
    sort_group=c('RGC', 'CycProg', 'Astrocyte', 'Epend', 'Neuron'),
    col_group=colorG1,
    show_gene=c('PLAGL2(MA1548.1)', 'Plagl1(MA1615.1)', 'PLAG1(MA0163.1)')
)

pdf(paste0(outdir,'/heatmap.expMarker.filtered.pdf'), width=2.2, height=3, useDingbats=FALSE)
draw(ht)
dev.off()






