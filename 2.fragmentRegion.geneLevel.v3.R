
library(Seurat)
library(Signac)
library(seuextra)
library(ggplot2)
library(patchwork)

rds <- 'snatac_integ/snatac_integ.withumap.rds'
fmeta <- 'out_celltype/metaData.cellType.xls'
fragment <- "combined_fragments/combined_fragments_sorted.tsv.gz"

celltype <- c('RGC', 'Neuron')

outdir <- 'out_peakRegion'



dir.create(outdir)

combined_snatac <- readRDS(rds)
DefaultAssay(combined_snatac) <- "peaks"

metadata <- read.table(fmeta, sep='\t', header=T, row.names=1)
combined_snatac$cell_type2 <- metadata$cell_type2[match(rownames(combined_snatac@meta.data), rownames(metadata))]

combined_snatac <- subset(combined_snatac, subset=cell_type2 %in% celltype)
combined_snatac <- subset(combined_snatac, subset=orig.ident %in% c('ZR0', 'ZR1', 'ZR2'))

Idents(combined_snatac) <- 'cell_type2'
levels(combined_snatac) <- celltype

# plot
gene <- "Ephb2"
outpdf <- paste0(outdir, "/fragments.", gene, ".pdf")
p <- FragmentDensityForGene(obj=combined_snatac, fragment=fragment, gene=gene)
ggsave(outpdf, width = 3, height = 2)

# plot
gene <- "Notch1"
outpdf <- paste0(outdir, "/fragments.", gene, ".pdf")
p <- FragmentDensityForGene(obj=combined_snatac, fragment=fragment, gene=gene)
ggsave(outpdf, width = 3, height = 2)




