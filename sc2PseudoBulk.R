# Hua Sun

library(Seurat)
library(data.table)
library(seuproc)
library(scqc)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(EnsDb.Hsapiens.v86)

set.seed(42)
options(future.globals.maxSize=100000000000)



finfo <- 'sample.info'
rds <- 'seurat5.1_v6.2/multiome_integrated.plus.rds'
outdir <- 'out_deseq2'
dir.create(outdir)

info <- fread(finfo)
# sample_id  group  type

info <- info[,c("id", "group")]
info$group <- factor(info$group)
# set DEG condition in 'group' - Relapse vs Primary
contrast_by <- c("group", "Relapse", "Primary")


seu <- readRDS(rds)
# Seurat 5.1.0
seu[['RNA']] <- split(seu[['RNA']], f=seu$orig.ident)


d_count <- NULL
for (x in info$id){
    print(x)

    # counts.sample
    x_count <- as.matrix(seu@assays$RNA@layers$counts.[[x]])
    x_sumcount <- as.data.frame(rowSums(x_count))
    colnames(x_sumcount) <- x

    d_count <- cbind(d_count, x_sumcount)
}


# ENSG & GeneSymbol mixed
features <- rownames(seu@assays$RNA@features)
rownames(d_count) <- features
# index(gene)  sample1 sample2 sample3 ....

# genesymbol <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembl_ids, keytype = "GENEID", columns = c("geneSymbol", "ensembl"))
# geneSymbol  ensembl


# check
all(rownames(info) %in% colnames(d_count))
# TRUE


# DESeq2
# https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
# https://hds-sandbox.github.io/bulk_RNAseq_course/develop/05c_count_normalization.html
dds <- DESeq2::DESeqDataSetFromMatrix(countData=d_count, colData=info, design=~group)

# Pre-filtering
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
dds <- DESeq2::estimateSizeFactors(dds)

# Normalization
normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
write.table(normalized_counts, file=paste0(outdir, "/normalized_counts.txt"), sep="\t", quote=F, col.names=NA)



# Differential expression analysis
dds <- DESeq2::DESeq(dds)
# Relapse vs Primary
res <- DESeq2::results(dds, contrast=contrast_by)
res <- res[order(res$padj),]

write.table(res, file=paste0(outdir, "/diffExp_results.xls"), sep="\t", quote=F, row.names=F)









