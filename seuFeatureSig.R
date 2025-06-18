# Hua Sun
# v3
	
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
assay <- 'SCT'

umap <- 'wnn.umap'
marker <- 'zr'

ref <- 'mouse'     # human/mouse
cutoff <- 0.5
point_size <- 0.1

fcell <- ''
sample_info <- ''  # sample group info
sgroup <- ''     # subgroup 
sort_info <- ''    # sample group (with header)

outdir <- '.'

GetoptLong(
    "rds=s",         "matrix file",
    "ref=s",         "ref",
    "cutoff=f",      "cutoff",
    "point_size=f",  "point_size",
    "marker=s",      "marker",
    "assay=s",       "assay",
    "fcell=s",       "target cell list",
    "umap=s",        "reduction",
    "sort_info=s",   "info for sample",
    "sample_info=s", "sample info for sample",
    "sgroup=s",      "subgroup",
    "outdir=s",      "output path"
)


fgene <- NULL
name <- NULL
label_name <- NULL



# ST-EPN-ZFTA-RELA FUS
if (marker == 'zr'){
    cutoff <- 0.2
    name = 'zr_score'

    fgene <- paste0(path, '/public_db/', 'ZRFus.', ref, '.gene')
    label_name <- 'ZR-Fus signal'
}


# ST-EPN-YAP1 FUS
if (marker == 'yap1'){
    cutoff <- 0.4
    name = 'zr_score'

    fgene <- paste0(path, '/public_db/', 'YAP1Fus.', ref, '.gene')
    label_name <- 'YAP1-Fus signal'
}



# GBM
if (marker == 'gbm'){
    cutoff <- 0.5
    name = 'gbm_score'

    fgene <- paste0(path, '/public_db/', 'GBM.manual.', ref, '.gene')
    label_name <- 'GBM signal'
}



#-------- Customized markers

# ST-EPN-ZFTA-RELA FUS
if (marker == 'custom.zr'){
    cutoff <- 0.6
    name = 'zr_score'

    fgene <- paste0(path, '/public_db/customized/', 'ST-EPN-ZR-Fus.20250513.', ref, '.gene')
    label_name <- 'ZR-Fus signal'
}

# ST-EPN-YAP1 FUS
if (marker == 'custom.yap1'){
    cutoff <- 0.5
    name = 'yap1_score'

    fgene <- paste0(path, '/public_db/customized/', 'ST-EPN-YAP1-Fus.20250513.', ref, '.gene')
    label_name <- 'YAP1-Fus signal'
}

# GBM
if (marker == 'custom.gbm'){
    cutoff <- 0.5
    name = 'gbm_score'

    fgene <- paste0(path, '/public_db/customized/', 'GBM.20250513.', ref, '.gene')
    label_name <- 'GBM signal'
}




######################################
##               Main
######################################

dir.create(outdir)

geneList <- readLines(fgene)
print(paste('input genes:', length(geneList)))

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


matchedFeatuers <- intersect(row.names(seu[[assay]]@data), geneList)
print(paste('matched genes:', length(matchedFeatuers)))

#DefaultAssay(seu) <- assay
# set slot as 'data'
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
write.table(metadata, file = paste0(outdir,'/metadata_with_featuresSig.xls'), sep = "\t", quote=FALSE, col.names = NA)



df <- metadata[,c('UMAP_1', 'UMAP_2', paste0(name,'1'))]
colnames(df) <- c('UMAP_1', 'UMAP_2', 'signal')
#UMAPPlot(df, outdir)
#UMAPPlot(df=df, min=-0.2, max=0.2, outdir=outdir)


# if set max & min
df$group <- ''
breaks <- NULL
low <- cutoff - 0.1

df$group[df$signal < low ] <- "Low"
df$group[df$signal > cutoff ] <- "High"
df$group[df$signal >= low & df$signal <= cutoff] <- "Middle"
breaks = c("High", "Middle", "Low")



p <- UMAPPlot_Score(df=df, breaks=breaks, label_name=label_name, point_size=point_size)
ggsave(paste0(outdir, '/umap.featuresSig.pdf'), width=5.1, height=4)


# sort sample
sorted_sample <- ''
if (nchar(sort_info) > 1){
    d_info <- fread(sort_info, select=c('sample', 'group'), data.table=F)
    sorted_sample <- d_info$sample
}
#print(sorted_sample)


# number of cells
dtemp <- df[df$signal >= low, ]
ratio_sig <- round(nrow(dtemp)/nrow(df) * 100, 1)
writeLines(c(nrow(dtemp), paste0(ratio_sig, '%')), paste0(outdir, '/cellRatio.positiveSignal.log'))















