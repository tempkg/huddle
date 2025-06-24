
# https://ajitjohnson.com/tsne-for-biologist-tutorial/
# tSNE vs PCA - tSNE is better
# https://distill.pub/2016/misread-tsne/

# info
# sample group

library(Seurat)
library(Rtsne)
library(ggplot2)
require(ggrepel)
library(RColorBrewer)
library(GetoptLong)

# integrated scRNA
rds <- 'snrna_harmony.rds'
fcell <- ''
perplexity <- 3
gene <- ''
info <- ''
outdir <- 'out_tsne.subgroup'

GetoptLong(
    "rds=s",         "matrix file",
    "fcell=s",       "barcode file",
    "perplexity=i",  "tsne perplexity",
    "gene=s",        "gene list file",
    "info=s",        "info",
    "outdir=s",      "output path"
)

dir.create(outdir)

seurat_obj <- readRDS(rds)
metadata <- seurat_obj@meta.data
print(nrow(metadata))

# extract target cells
if (nchar(fcell)>2){
	target_cells <- readLines(fcell)
	metadata <- metadata[target_cells,]
	print(nrow(metadata))
}



obj.list <- SplitObject(seurat_obj, split.by = "Sample")

features <- c()
if (gene == ''){
	features <- SelectIntegrationFeatures(obj.list, nfeatures = 2000)
} else {
	#geneList <- readLines(gene)
	pipe_g <- pipe(paste0('cut -f1 ', gene))
	geneList <- readLines(pipe_g)
	close(pipe_g)

	all_featuers <- row.names(seurat_obj@assays$SCT@data)
	features <- intersect(all_featuers, geneList)
}


# exp matrix
exp_mtx <- as.matrix(seurat_obj@assays$SCT@data[features,])


df <- c()
for (sample in unique(metadata$Sample)){
	print(sample)
	barcode <- row.names(metadata[metadata$Sample == sample,])

	tar_mtx <- exp_mtx[,barcode]
	
	df <- cbind(df, sample = rowMeans(tar_mtx))
}

colnames(df) = unique(metadata$Sample)
converted <- t(df)



## Run the t-SNE algorithm and store the results into an object called tsne_results
#sample_size <- length(unique(metadata$Sample))
#perplexity <- 1
#if (sample_size > 10){ perplexity <- 3 }

set.seed(42)
# dims should be either 1, 2 or 3
tsne_res <- Rtsne(converted, dims = 2, perplexity=perplexity, check_duplicates = FALSE)
df_tsne <- as.data.frame(tsne_res$Y)
row.names(df_tsne) <- row.names(converted)
colnames(df_tsne) <- c('tSNE_1', 'tSNE_2')

if (info != ''){
	d_info <- read.table(info, sep='\t', header=T)
	colnames(d_info) <- c('sample', 'group')
	df_tsne <- merge(df_tsne, d_info, by.x='row.names', by.y='sample')
	rownames(df_tsne) <- df_tsne[,1]
	df_tsne[,1] <- NULL
}

write.table(df_tsne, paste0(outdir, "/tsne.samples.xls"), sep='\t', quote=F, col.name=NA)


p <- ''
if (info != ''){
	getPalette = colorRampPalette(brewer.pal(9, "Set1"))
	colourCount = length(unique(df_tsne$group))

	p <- ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=group)) + 
		geom_point(size=3) + 
		theme_classic(base_line_size=0.2) +
		theme(legend.title=element_blank()) +
		scale_color_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(colourCount)) 
} else {
	p <- ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2)) + 
		geom_point(color='#C0392B') + 
		theme_classic(base_line_size=0.2)
}

p <- p + labs(x='tSNE 1', y='tSNE 2')


# without label
pdf(paste0(outdir, "/tsne.samples.pdf"), width=5, height=3.5)
print(p)
dev.off()

# with label
p <- p + geom_text_repel(aes(label = rownames(df_tsne)), size=3, color='#808B96', min.segment.length = 0, box.padding = 0.2)
pdf(paste0(outdir, "/tsne.samples.withlabel.pdf"), width=5, height=3.5)
print(p)
dev.off()




#plot(tsne_res$Y, col = "black", bg= tsne_res$origD, pch = 21, cex = 1)





# PCA plot 
#data.pca <- princomp(df)
#plot(data.pca$loadings[, 1:2])





