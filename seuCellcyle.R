# Hua Sun
# v2
	
# Set Library
library(GetoptLong)
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(this.path)
library(RColorBrewer)

path <- dirname(this.path())
fpath <- paste0(path, '/src/')
r_source <- list.files(fpath, recursive = T, full.names = T, pattern = ".R")
invisible(lapply(r_source, source))


rds <- ''
fcell <- ''

umap <- 'wnn.umap'
#seq <- 'multiome'   # multiome/scrna
name <- 'cc_score'
label_name <- ''   # legend title
sort_info <- ''    # sample group (with header)

ref <- 'human'     # human/mouse
cutoff <- '0.2'   # '0.2' or '0.3'
assay <- 'SCT'
point_size <- 0.1 

sample_info <- ''  # sample group info
sgroup <- ''     # subgroup 

outdir <- 'out_cellcycle'

GetoptLong(
    "rds=s",         "matrix file",
    "ref=s",         "ref",
    "cutoff=s",      "cutoff",
    "assay=s",       "assay",
    "point_size=f",  "point size",
    "fcell=s",       "target cell list",
    "umap=s",        "reduction",
    "name=s",        "name of the score",
    "label_name=s",  "legend title",
    "sort_info=s",        "info for sample",
    "sample_info=s", "sample info for sample",
    "sgroup=s",      "subgroup",
    "outdir=s",      "output path"
)



# cell cycle
fg2m <- paste0(path, '/public_db/', 'cellcyle/CellCycle.seurat.', ref, '.2019.g2m.gene')
fs <- paste0(path, '/public_db/', 'cellcyle/CellCycle.seurat.', ref, '.2019.s.gene')

if (label_name == ''){
    label_name <- 'Cell cycle signal'
}

g2m.genes <- readLines(fg2m)
s.genes <- readLines(fs)



######################################
##               Main
######################################

dir.create(outdir)

seu <- readRDS(rds)

# cell barcode file
if (nchar(fcell)>2){
    barcodes <- readLines(fcell)
    seu <- subset(seu, cells=barcodes)
}

# sample info file
if (nchar(sample_info)>2){
    sample_info <- read.table(sample_info, sep='\t', header=T)
    colnames(sample_info) <- c('sample', 'group', 'subgroup')
    # add subgroup
    seu@meta.data$subgroup <- sample_info$subgroup[match(seu@meta.data$Sample, sample_info$sample)]
    print(table(seu@meta.data$subgroup))

    # extract subgroup
    if (nchar(sgroup) > 1){ seu <- subset(x = seu, subset = subgroup==sgroup) }
    print(table(seu@meta.data$subgroup))
}


# CellCycleScoring
seu <- CellCycleScoring(
            object = seu,
            s.features = s.genes,
            g2m.features = g2m.genes
        )


# CellCycleScoring2 for all cell cycle genes
cc.genes <- c(s.genes, g2m.genes)
matchedFeatuers <- intersect(row.names(seu[[assay]]@data), cc.genes)

seu <- AddModuleScore(
    object = seu,
    assay = assay,
    slot = "data",
    features = list(matchedFeatuers),
    name = name
)

# umap
seu$UMAP_1 <- seu@reductions[[umap]]@cell.embeddings[,1]
seu$UMAP_2 <- seu@reductions[[umap]]@cell.embeddings[,2]

metadata <- seu@meta.data
write.table(metadata, file = paste0(outdir,'/metadata_with_cellcycleSig.xls'), sep = "\t", quote=FALSE, col.names = NA)



df <- metadata[,c('UMAP_1', 'UMAP_2', paste0(name,'1'))]
colnames(df) <- c('UMAP_1', 'UMAP_2', 'signal')

# if set max & min
df$group <- ''
breaks <- NULL
pos_cutoff <- NULL  # for counting singal cells

# 0.2
if (cutoff=='0.2'){
    df$group[df$signal < 0.1 ] <- "< 0.1"
    df$group[df$signal > 0.2 ] <- "> 0.2"
    df$group[df$signal >= 0.1 & df$signal <= 0.2] <- "0.1-0.2"
    breaks = c("> 0.2", "0.1-0.2", "< 0.1")

    pos_cutoff = .1
}

# 0.3
if (cutoff=='0.3'){
    df$group[df$signal < 0.2 ] <- "< 0.2"
    df$group[df$signal > 0.3 ] <- "> 0.3"
    df$group[df$signal >= 0.2 & df$signal <= 0.3] <- "0.2-0.3"
    breaks = c("> 0.3", "0.2-0.3", "< 0.2")

    pos_cutoff = .2
}

# 0.4
if (cutoff=='0.4'){
    df$group[df$signal < 0.3 ] <- "< 0.3"
    df$group[df$signal > 0.4 ] <- "> 0.4"
    df$group[df$signal >= 0.3 & df$signal <= 0.4] <- "0.3-0.4"
    breaks = c("> 0.4", "0.3-0.4", "< 0.3")

    pos_cutoff = .3
}


p <- UMAPPlot_Score(df=df, breaks=breaks, point_size=point_size)
ggsave(paste0(outdir, '/umap.featuresSig.pdf'), width=5.1, height=4)


# sort sample
sorted_sample <- ''
if (nchar(sort_info) > 1){
    d_info <- fread(sort_info, select=c('sample', 'group'), data.table=F)
    sorted_sample <- d_info$sample
}
#print(sorted_sample)



# show split per sample 
if ('Sample' %in% colnames(metadata)){
    df$Sample <- metadata$Sample
    p <- SignalUMAP_perSample(df=df, breaks=breaks, point_size=point_size, sorted_sample=sorted_sample, cutoff=cutoff)
    ggsave(paste0(outdir, '/umap.featuresSig.splited.pdf'), width = 6, height = 2)
} else {
    if (length(unique(metadata$orig.ident)) > 1){
        df$Sample <- metadata$orig.ident
        p <- SignalUMAP_perSample(df=df, breaks=breaks, point_size=point_size, sorted_sample=sorted_sample, cutoff=cutoff)
        ggsave(paste0(outdir, '/umap.featuresSig.splited.pdf'), width = 6, height = 2)
    }
}



# number of cells
dtemp <- df[df$signal >= pos_cutoff, ]
ratio_sig <- round(nrow(dtemp)/nrow(df) * 100, 1)
writeLines(c(nrow(dtemp), paste0(ratio_sig, '%')), paste0(outdir, '/cellRatio.positiveSignal.log'))




