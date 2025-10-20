
library(GetoptLong)
library(Seurat)
library(seuextra)
library(scapekit)
library(dplyr)
library(ggplot2)
library(data.table)
library(this.path)
library(RColorBrewer)
library(yaml)

path <- dirname(this.path())


rds <- ''
assay <- 'SCT'

umap <- 'wnn.umap'
marker <- 'zr-fus'

ref <- 'mouse'     # human/mouse

sort_info <- ''    # sample group (with header)

outdir <- 'out_zrFusSig'

GetoptLong(
    "rds=s",         "matrix file",
    "ref=s",         "ref",
    "marker=s",      "marker",
    "assay=s",       "assay",
    "umap=s",        "reduction",
    "sort_info=s",   "info for sample",
    "outdir=s",      "output path"
)


dyml <- yaml.load_file(paste0(path, '/config.yml'))
name <- dyml[[marker]]['name']
label_name <- dyml[[marker]]['label']
prefix <- dyml[[marker]]['db_prefix']
cutoff <- dyml[[marker]][[ref]]
fgene <- paste0(path, '/db/', prefix, '.', ref, '.gene')



dir.create(outdir)

geneList <- readLines(fgene)
print(paste('input genes:', length(geneList)))

seu <- readRDS(rds)

seu <- SeuratAddScore(obj=seu, assay='SCT', name=name, features=geneList)
# umap
seu$UMAP_1 <- seu@reductions[[umap]]@cell.embeddings[,1]
seu$UMAP_2 <- seu@reductions[[umap]]@cell.embeddings[,2]

metadata <- seu@meta.data
write.table(metadata, file = paste0(outdir,'/metadata_with_featuresSig.xls'), sep = "\t", quote=FALSE, col.names = NA)



df <- metadata[,c('orig.ident', 'UMAP_1', 'UMAP_2', paste0(name,'1'))]
colnames(df) <- c('sample', 'UMAP_1', 'UMAP_2', 'signal')
#UMAPPlot(df, outdir)
#UMAPPlot(df=df, min=-0.2, max=0.2, outdir=outdir)


# if set max & min
df$group <- ''

low <- cutoff - 0.1

df$group[df$signal < low ] <- "Low"
df$group[df$signal > cutoff ] <- "High"
df$group[df$signal >= low & df$signal <= cutoff] <- "Middle"



p <- UMAPSignal(data=df, signal='signal', group='group', breaks=c("High", "Middle", "Low"),
    colors=c("High"="#CB4335", "Middle"="#2E86C1", "Low"="#D7DBDD"), 
    legend_title=label_name)
ggsave(paste0(outdir, '/umap.featuresSig.pdf'), width=5.1, height=4)


# sort sample
sorted_sample <- ''
if (nchar(sort_info) > 1){
    d_info <- fread(sort_info, select=c('V1', 'V2'), data.table=F, header=F)
    sorted_sample <- d_info$V1
}
#print(sorted_sample)


# number of cells
dtemp <- df[df$signal >= low, ]
ratio_sig <- round(nrow(dtemp)/nrow(df) * 100, 1)
writeLines(c(nrow(dtemp), paste0(ratio_sig, '%')), paste0(outdir, '/cellRatio.positiveSignal.log'))







