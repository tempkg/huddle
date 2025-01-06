# boxplot for merged raw data metadata 


library(data.table)
library(ggplot2)

fmeta <- 'merged.rawdata.metadata.xls'

meta <- fread(fmeta, data.table=F)

columns <- c('nCount_ATAC', 'nucleosome_signal', 'TSS.enrichment', 'blacklist_ratio', 'pct_reads_in_peaks', 'nCount_RNA', 'nFeature_RNA', 'percent.mt')



# remove one max data
DelMax <- function(df, group, col){
    df_filtered <- df %>% group_by(get(group)) %>% filter(get(col) < max(get(col)))

    return(df_filtered)
}



# plot
pdf('rawqc.pdf', w=2, h=3)
for (x in columns){
    print(x)
    df2 <- DelMax(meta, 'orig.ident', x)
    
    p <- ggplot(df2, aes_string(x='orig.ident', y=x)) + geom_boxplot(outlier.size = 0.1) + labs(x='') + theme(text = element_text(size = 6)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

    print(p)
    
}
dev.off()


