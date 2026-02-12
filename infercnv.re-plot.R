library(infercnv)
library(GetoptLong)


rds <- 'run.final.infercnv_obj'
title <- 'InferCNV'
outdir <- "./custom_ordered"


GetoptLong(
    "rds=s",     "seurat obj",
    "title=s",   "infercnv obj",
    "outdir=s",  "output path"
)


# Read your infercnv object
infercnv_obj <- readRDS(rds)

# Define your custom chromosome order
custom_chr_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                      "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                      "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", 
                      "chr21", "chr22")

# Reorder - THIS IS THE KEY LINE
infercnv_obj@gene_order <- infercnv_obj@gene_order[order(
  factor(infercnv_obj@gene_order$chr, levels = custom_chr_order), 
  infercnv_obj@gene_order$start
), ]

# Reorder expression data to match
infercnv_obj@expr.data <- infercnv_obj@expr.data[rownames(infercnv_obj@gene_order), ]


#infercnv_obj@expr.data[infercnv_obj@expr.data < 0.8] <- 0.8
#infercnv_obj@expr.data[infercnv_obj@expr.data > 1.2] <- 1.2


# https://www.rdocumentation.org/packages/infercnv/versions/1.3.3/topics/plot_cnv
plot_cnv(infercnv_obj,
        hclust_method='ward.D2',
        cluster_by_groups=TRUE,
        color_safe_pal=FALSE,
        title = title,
        out_dir = outdir,
        output_filename="custom.infercnv",
        output_format="pdf"
)





