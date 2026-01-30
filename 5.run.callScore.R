
library(Seurat)
library(scatlasR)
library(ggplot2)
library(GetoptLong)



rds <- 'multiome.rds'
f <- 'hs.EPN-ZR.txt'
name <- 'zr_fus'
umap <- 'wnn.umap'
cutoff <- 0.1
assay <- 'SCT'
w <- 3.5
h <- 3
outdir <- 'out_score'

GetoptLong(
    "f=s",      "file",
    "rds=s",    "rds file",
    "name=s",   "name score",
    "cutoff=f",   "cutoff",
    "assay=s",    "assay",
    "umap=s",   "umap",
    "w=i",      "width",
    "h=i",      "height",
    "outdir=s",    "out pdf"
)

dir.create(outdir)

features <- readLines(f)

seu <- readRDS(rds)
seu <- SeuratAddScore(obj=seu, assay=assay, name=name, features=features)


# show features in UMAP
p <- FeaturePlot(seu, features=paste0(name, '1'), reduction=umap, order=T, min.cutoff='q1')
# 3x3 image
ggsave(paste0(outdir, '/', name, '.score_umap.pdf'), w=w, h=h)




df <- seu@meta.data[,c('orig.ident', 'seurat_clusters', paste0(name, '1'))]
colnames(df) <- c('sample', 'cluster', 'signal')
df$UMAP_1 <- seu@reductions[[umap]]@cell.embeddings[,1]
df$UMAP_2 <- seu@reductions[[umap]]@cell.embeddings[,2]

high <- cutoff + 0.1
df$group[df$signal < cutoff ] <- "Low"
df$group[df$signal > high ] <- "High"
df$group[df$signal >= cutoff & df$signal <= high] <- "Middle"
# output
write.table(df, file=paste0(outdir, '/', name, '.seuscore.xls'), sep='\t', quote=F, col.names=NA)


# number of cells
n_sig <- nrow(df[df$signal >= cutoff, ])
n_total <- nrow(df)
ratio_sig <- round(n_sig/n_total * 100, 1)
writeLines(c(n_sig, paste0(ratio_sig, '%')), file=paste0(outdir, '/', name, '.cellRatio.log'))

# plot
p <- UMAPSignal(data=df, signal='signal', group='group', breaks=c("High", "Middle", "Low"),
        colors=c("High"="#CB4335", "Middle"="#2E86C1", "Low"="#D7DBDD"), 
        legend_title='')
ggsave(paste0(outdir, '/', name, '.umap.signal.pdf'), width=2.6, height=1.8)








