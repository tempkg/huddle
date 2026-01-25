
library(Seurat)
library(scatlasR)
library(scapekit)
library(ggplot2)
library(GetoptLong)


rds <- ''
fdb <- ''
name <- 'zr-fus'
assay <- 'SCT'
outdir <- 'out_zrFusSig'

GetoptLong(
    "rds=s",         "matrix file",
    "fdb=s",         "fdb",
    "name=s",        "name",
    "assay=s",       "assay",
    "outdir=s",      "output path"
)


seu <- readRDS(rds)
gene_list <- readLines(fdb)

seu <- SeuratAddScore(obj=seu, name='zr_fus', features=gene_list)

df <- seu@meta.data[,c('orig.ident', 'seurat_clusters', paste0(name, '1'))]
colnames(df) <- c('sample', 'cluster', 'signal')
df$UMAP_1 <- seu@reductions[[umap]]@cell.embeddings[,1]
df$UMAP_2 <- seu@reductions[[umap]]@cell.embeddings[,2]

df$group[df$signal < 0.1 ] <- "Low"
df$group[df$signal > 0.2 ] <- "High"
df$group[df$signal >= 0.1 & df$signal <= 0.2] <- "Middle"

# output
write.table(df, file=paste0(outdir, '/seuscore.', name, '.xls'), sep='\t', quote=F, col.names=NA)


# plot
p <- UMAPSignal(data=df, signal='signal', group='group', breaks=c("High", "Middle", "Low"),
    colors=c("High"="#CB4335", "Middle"="#2E86C1", "Low"="#D7DBDD"), 
    legend_title=label_name)
ggsave(paste0(outdir, '/umap.', name, '.pdf'), width=5.1, height=4)






