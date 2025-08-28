
library(Seurat)
library(Signac)
library(seuextra)
library(ggplot2)
library(patchwork)

rds <- 'multiome_integrated_plus.seurat5.rds'
fmeta <- 'out_celltype/cluster_cellType.corrected.tsv'
fragment <- "combined_fragments/combined_fragments_sorted.tsv.gz"

celltype <- c('RGC', 'Neuron')

outdir <- 'out_peakRegion'



dir.create(outdir)

 <- readRDS(rds)
DefaultAssay(seu) <- "peaks"

metadata <- fread(fmeta, data.table=F)
seu$cell_type2 <- metadata$cell_type2[match(seu$seurat_clusters, metadata$seurat_clusters)]

seu <- subset(seu, subset=cell_type2 %in% celltype)

Idents(seu) <- 'cell_type2'
levels(seu) <- celltype

# plot
gene <- "Ephb2"
outpdf <- paste0(outdir, "/fragments.", gene, ".pdf")
p <- FragmentDensityForGene(obj=seu, fragment=fragment, gene=gene)
ggsave(outpdf, width = 3, height = 2)

# plot
gene <- "Notch1"
outpdf <- paste0(outdir, "/fragments.", gene, ".pdf")
p <- FragmentDensityForGene(obj=seu, fragment=fragment, gene=gene)
ggsave(outpdf, width = 3, height = 2)




