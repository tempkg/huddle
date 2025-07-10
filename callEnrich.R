# Hua Sun


library(GetoptLong)
library(data.table)
library(estEnrich)
library(dplyr)
library(this.path)
path <- dirname(this.path())


f <- 'out_expMarkers.seurat_clusters/cell_type2.findallmarker.geneExp.xls'
gmt <- paste0(path, '/db/mm.test.gmt')
method <- 'ora'
prefix <- 'test'
outdir <- 'out_enrich'

GetoptLong(
    "f=s",       "data file",
    "method=s",  "method",
    "gmt=s",     "gmt db",
    "prefix=s",  "prefix",
  "outdir=s",  "out dir"
)



dir.create(outdir)

expMarker <- read.table(f, sep='\t', header=T)


df <- NULL
for (clus in unique(expMarker$cluster)){
    dexp <- expMarker %>% filter(cluster == clus)

    d_gene_symbol <- dexp[,c('gene', 'avg_log2FC')]

    if (method=='ora'){
        d_sig <- RunORA(genes=d_gene_symbol$gene, fgmt=gmt)
    } else {
        d_sig <- RunGSEA(data=d_gene_symbol, fgmt=gmt)
            
    }

    if (nrow(d_sig) == 0){ next }

    d_sig$cluster <- clus
    df <- rbind(df, d_sig)
}


if (method=='ora'){
    write.table(as.data.frame(df), paste0(outdir, '/', prefix, '.enrichment.ora.xls'), sep='\t', quote=F, row.names=F)
} else {
    fwrite(df, file=paste0(outdir, '/', prefix, '.enrichment.gsea.xls'), sep="\t", sep2=c("", " ", ""))
}













