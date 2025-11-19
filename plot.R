
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)



# barplot - summmary MalignantCellPercentagePerSample
MalignantCellPercentagePerSample <- function(df){
    d_cellstatus <- df %>% dplyr::group_by(sample) %>% dplyr::count(cell_status) %>% dplyr::mutate(percent = n/sum(n) * 100)
    # Sample   cell_status   n   percent
    d_cellstatus$percent <- round(d_cellstatus$percent, digits=1)

    #write.table(df, paste0(outdir, '/malignantCells.sample.xls'), sep='\t', quote=F, col.names=NA)


    p <- ggplot(d_cellstatus, aes(fill=cell_status, y=percent, x=sample)) +
          geom_bar(stat="identity")+
          geom_text(aes(label = percent), position = position_stack(vjust = 0.5), color="white", size=2) +
          #scale_fill_brewer(palette="Paired")+
          scale_fill_manual(values=c('#6C3483', '#E1C2ED'))+
          theme_minimal() +
          xlab("") + ylab('Percentage of cells (%)') +
          theme(legend.title=element_blank(), 
                legend.key.size = unit(3, 'mm')) +
          theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=8))

    p
}



