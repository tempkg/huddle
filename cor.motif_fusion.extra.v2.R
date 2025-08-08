
library(GetoptLong)
library(Seurat)
library(scapekit)
library(ggplot2)


motif_id <- 'MA1615.1'
show_celltypes <- c('RGC', 'CycProg', 'Astrocyte', 'Epend', 'Neuron')

fct <- 'cluster_ct.xls'
fmeta <- 'fus_meta.xls'

dir.create(outdir)

info <- fread(fct, data.table=F)
d_fus <- read.table(fmeta, sep='\t', header=T, row.names=1)

seu <- readRDS(rds)


# seurat_clusters  cell_type2
seu$cell_type2 <- info$group[match(seu$seurat_clusters), info$seurat_clusters)]


d_motif <- as.data.frame(seu@assays$chromvar@data[motif_id,])
colnames(d_motif) <- motif_id
d_motif$cell_type2 <- seu$cell_type2[match(rownames(d_motif), rownames(seu@meta.data))]
d_motif$fus_score1 <- d_fus$fus_score1[match(row.names(d_motif), row.names(d_fus))]


d_motif2 <- filter(d_motif, cell_type2 %in% show_celltypes)
p <- ScatterPerGroupWithCor(d=d_motif2, group='cell_type2', x='fus_score1', y=motif_id,  
        x_lab='ZR-Fusion signal', y_lab=paste0(gene, ' motif score'))


ggsave(paste0(motif_id, '.cor.motif_fus.perCellType.pdf'), w=2.7, h=2.5)








