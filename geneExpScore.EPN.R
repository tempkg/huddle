
# Calculate module scores for feature expression programs in single cells
# AddModuleScore (Seurat function)
	
# Set Library
library(GetoptLong)
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(this.path)


script_dir <- dirname(this.path())
fgene <- paste0(script_dir, '/db/EPN.2015cc.hs_gene.info')

rds <- ''
umap <- 'wnn.umap'
cutoff <- 0.2
info <- ''    # sample group (with header)

label_name <- 'EPN signal'  # legend title
outdir <- '.'

GetoptLong(
    "rds=s",         "matrix file",
    "umap=s",        "reduction",
    "cutoff=f",      "cutoff",
    "fgene=s",       "gene info",
    "info=s",        "sample info",
    "label_name=s",  "legend title",
    "split",         "split gene group",
    "outdir=s",      "output path"
)





######################################
##            Set Func.
######################################

UMAPPlot_Score <- function(df=NULL, cutoff=0.2, outdir='.'){
    df <- df[order(df$EPN_score, decreasing=FALSE), ]
    
    # if set max & min
    df$signal <- ''
    cutoff1 <- paste0(">= ", cutoff)
    cutoff2 <- paste0("< ", cutoff)
    
    df$signal[df$EPN_score >= cutoff ] <- cutoff1
    df$signal[df$EPN_score < cutoff ] <- cutoff2
    
    p <- ggplot(df,aes(x=UMAP_1, y=UMAP_2, color=signal)) + 
            geom_point(size=0.3) +
            theme_classic() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
            labs(x = "UMAP 1", y="UMAP 2") +
            scale_color_manual(name = label_name, breaks = c(cutoff1, cutoff2), values=c("#CB4335", "#D7DBDD")) +
            guides(color = guide_legend(override.aes = list(size = 2))) +
            theme(legend.title = element_text(size=9))


    pdf(paste0(outdir, '/umap.epn_signal.pdf'), width=4, height=3, useDingbats=FALSE)
    print(p)
    dev.off()
}



## UMAP per sample
SignalUMAP_perSample <- function(df, cutoff, sorted_sample='', outfile='')
{
    df <- df[order(df$EPN_score, decreasing=FALSE), ]
    
     # if set max & min
    df$signal <- ''
    cutoff1 <- paste0(">= ", cutoff)
    cutoff2 <- paste0("< ", cutoff)
    
    df$signal[df$EPN_score >= cutoff ] <- cutoff1
    df$signal[df$EPN_score < cutoff ] <- cutoff2

    if (length(sorted_sample) > 0){ 
        df$Sample <- factor(df$Sample, levels=sorted_sample) 
    }


    p <- ggplot(df,aes(x=UMAP_1, y=UMAP_2, color=signal)) + 
            geom_point(size=0.2) +
            theme_classic() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
            labs(x = "UMAP 1", y="UMAP 2") +
            scale_color_manual(name = label_name, breaks = c(cutoff1, cutoff2), values=c("#CB4335", "#D7DBDD")) +
            guides(color = guide_legend(override.aes = list(size = 2))) +
            theme(legend.title = element_text(size=9))

    # per sample
    p <- p + facet_wrap(~Sample, ncol = 5) + 
        theme(strip.background = element_blank(), strip.text = element_text(size=10))
    

    height = 2
    sample_size <- length(unique(df$Sample))
    if (sample_size > 4){ height <- height * round(sample_size/5 + 0.4) * 0.7 }
    pdf(paste0(outdir, '/umap.epn_signal.splited.pdf'), width = 6, height = height, useDingbats=FALSE)
    print(p)
    dev.off()
}




## gradient
UMAPPlot_Gradient <- function(df=NULL, sorted_sample='', outdir='.')
{
    df <- df[order(df$EPN_score, decreasing=FALSE), ]

    if (length(sorted_sample) > 0){ 
        df$Sample <- factor(df$Sample, levels=sorted_sample) 
    }
    
    #minVal <- min(df$signal)
    df$EPN_score[df$EPN_score < 0] <- 0

    p <- ggplot(df,aes(x=UMAP_1, y=UMAP_2, color=EPN_score)) + 
            geom_point(size=0.3) +
            theme_classic() +
            scale_colour_viridis_c(option = "plasma") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
            theme(legend.key.width=unit(5,"mm")) +
            labs(x = "UMAP 1", y="UMAP 2") 
            #+
            #scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), name = "Score", limits = c(minVal, maxVal))

    pdf(paste0(outdir, '/umap.epn_signal.gradient.pdf'), width = 4, height = 3, useDingbats=FALSE)
    print(p)
    dev.off()


    # per sample
    p <- p + facet_wrap(~Sample, ncol = 5) + 
            theme(strip.background = element_blank(), strip.text = element_text(size=10))

    height = 2
    sample_size <- length(unique(df$Sample))
    if (sample_size > 4){ height <- height * round(sample_size/5 + 0.4) * 0.7 }
    pdf(paste0(outdir, '/umap.epn_signal.gradient.splited.pdf'), width = 6.5, height = height, useDingbats=FALSE)
    print(p)
    dev.off()
}




######################################
##               Main
######################################


dir.create(outdir)

seurat_obj <- readRDS(rds)

# umap
seurat_obj$UMAP_1 <- seurat_obj@reductions[[umap]]@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions[[umap]]@cell.embeddings[,2]
# out umap
write.table(seurat_obj@meta.data, file = paste0(outdir,'/metadata_with_featuresSig.xls'), sep = "\t", quote=FALSE, col.names = NA)


# feature
gene_info <- read.table(fgene, sep='\t', header=F)
colnames(gene_info) <- c('gene', 'group')

print(paste('input genes:', nrow(gene_info)))

all_featuers <- row.names(seurat_obj[['SCT']]@data)

df <- c()

if (split){
    gene_group <- unique(gene_info$group)
    for (subgroup in gene_group){
        geneset <- unique(gene_info$gene[gene_info$group==subgroup])
        featuers <- intersect(all_featuers, geneset)
        print(length(featuers))
        if (length(featuers) == 0){ 
            gene_group <- gene_group[ ! gene_group == subgroup ]
            next
        }

        seurat_obj <- AddModuleScore(
            object = seurat_obj,
            assay = 'SCT',
            features = list(featuers),
            name = subgroup
        )
    }
    score_col <- paste0(gene_group, 1)
    df <- seurat_obj@meta.data[,c('UMAP_1', 'UMAP_2', 'seurat_clusters', score_col)]
    df$EPN_score <- apply(data.frame(df[,4:ncol(df)]), 1, max, na.rm=TRUE)
} else {
    featuers <- intersect(all_featuers, unique(gene_info$gene))
    seurat_obj <- AddModuleScore(
            object = seurat_obj,
            assay = 'SCT',
            features = list(featuers),
            name = 'EPN_score'
        )
    df <- seurat_obj@meta.data[,c('UMAP_1', 'UMAP_2', 'seurat_clusters', 'EPN_score1')]
    colnames(df) <- c('UMAP_1', 'UMAP_2', 'seurat_clusters', 'EPN_score')
}


# umap score
UMAPPlot_Score(df, cutoff, outdir)



# sort sample
sorted_sample <- 'na'
if (info != ''){
    d_info <- read.table(info, sep='\t', header=T)
    colnames(d_info) <- c('sample', 'group')
    sorted_sample <- d_info$sample
}


# show split per sample 
if ('Sample' %in% colnames(seurat_obj@meta.data)){
    df$Sample <- seurat_obj$Sample
    SignalUMAP_perSample(df, cutoff, sorted_sample, outdir)
} else {
    if (length(unique(seurat_obj@meta.data$orig.ident)) > 1){
        df$Sample <- seurat_obj$orig.ident
        SignalUMAP_perSample(df, cutoff, sorted_sample, outdir)
    }
}

# plot with gradient
UMAPPlot_Gradient(df, sorted_sample, outdir)





