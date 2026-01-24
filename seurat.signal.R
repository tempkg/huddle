library(GetoptLong)
library(Seurat)
library(seuextra)
library(dplyr)
library(ggplot2)
library(data.table)

rds <- 'multiome_integrated.plus.rds'
gene <- 'PTPRC'
assay <- 'SCT'

GetoptLong(
    "rds=s",       "matrix file",
    "gene=s",      "gene",
    "assay=s",      "assay"
)


seu <- readRDS(rds)


vec <- GeneSignalPropotion(obj=seu, assay=assay, idents='Sample', gene=gene)
df <- data.frame(vec)
colnames(df) <- gene

print(df)



