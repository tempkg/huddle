
library(dplyr)
library(data.table)
library(stringr)
library(readxl)
library(GetoptLong)


f <- 'cell_type2.findallmarker.geneExp.xls'
db <- 'publicdb.human.20250802.xlsx'
out <- 'confirmMarkers.xls'

GetoptLong(
    "f=s",       "file",
    "db=s",     "database",
    "out=s",  "out"
)

dir.create(outdir)

diff <- fread(f, data.table=F)
d_db <- read_excel(db)


df <- NULL
for (x in unique(diff$cluster)){
    geneset <- diff$gene[diff$cluster == x]

    for (ct in unique(d_db$cellName)){
        markers <- d_db$geneSymbolmore1[d_db$cellName == ct]
        markers <- str_split(markers, ',')[[1]]
        overlap_genes <- intersect(geneset, markers)
        geneSymbol <- paste(overlap_genes, collapse = ",")
        avg_log2fc <- round(mean(diff$avg_log2FC[diff$cluster == x & diff$gene %in% overlap_genes]), 3)
        n <- length(overlap_genes)
        if (n == 0){ next }
        new_row <- c(cluster=x, cellType=ct, geneCount=n, avg_log2fc=avg_log2fc, geneSymbol=geneSymbol)

        df <- rbind(df, new_row)
    }
}

df <- as.data.frame(df)
df_sorted <- df[order(df$cluster, as.numeric(df$geneCount), decreasing = TRUE), ]
rownames(df_sorted) <- NULL

write.table(df, file=out, sep='\t', quote=F, row.names=F)







