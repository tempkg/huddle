
library(dplyr)
library(data.table)
library(stringr)
library(readxl)
library(this.path)
library(GetoptLong)

path <- dirname(this.path())
fpath <- paste0(path, '/db/')
db_path <- paste0(path, '/db')

f <- 'seurat_clusters.findallmarker.geneExp.xls'
db <- ''
out <- 'confirmMarkers.publicdb_20250802.xls'

GetoptLong(
    "f=s",       "file",
    "db=s",     "database",
    "out=s",  "out"
)


diff <- fread(f, data.table=F)

if (db == ''){ db <- paste0(db_path, '/publicdb.human.20250802.xlsx')}
d_db <- readxl::read_excel(db)


df <- NULL
for (x in unique(diff$cluster)){
    geneset <- diff$gene[diff$cluster == x]

    for (ct in unique(d_db$cellName)){
        markers <- d_db$geneSymbol[d_db$cellName == ct]
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

# show top10
d_top10 <- df %>% group_by(cluster) %>% arrange(desc(as.numeric(geneCount))) %>% slice_head(n = 10)


write.table(d_top10, file=out, sep='\t', quote=F, row.names=F)







