
library(Seurat)
library(seuextra)
library(ggplot2)
library(GetoptLong)


rds <- 'seurat5.1_v6.2/multiome_integrated.plus.rds'
umap <- 'umap'
gene <- 'Ptprc'
min_cutoff <- NA
max_cutoff <- NA
outdir <- 'out_geneExp'

GetoptLong(
    "rds=s",    ".rds seurat file",
    "umap=s",   "umap",
    "gene=s",  "gene name",
    "min_cutoff=s",  "min_cutoff",
    "max_cutoff=s",  "max_cutoff",
    "outdir=s",  "outdir"
)


dir.create(outdir)


seu <-readRDS(rds)

# colors=c('#E0E0E0', 'red')
#p <- Seurat_FeaturePlot(obj=seu, features=gene, reduction=umap, min_cutoff = 'q10', max_cutoff = 'q90')
p <- Seurat_FeaturePlot(obj=seu, features=gene, reduction=umap, min_cutoff = min_cutoff, max_cutoff = max_cutoff)
ggsave(paste0(outdir, '/geneExp.', gene, '.umap.pdf'), w=2.4, h=2)


#p <- scCustom_FeaturePlot(obj=seu, features='Plagl1', reduction='wnn.umap')
#ggsave('geneExp.plagl1.umap.sccustom.pdf', w=2.4, h=2)




