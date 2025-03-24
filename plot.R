
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)



# barplot - summmary MalignantCellPercentagePerSample
MalignantCellPercentagePerSample <- function(df, outdir){
    d_cellstatus <- df %>% dplyr::group_by(Sample) %>% dplyr::count(cell_status) %>% dplyr::mutate(percent = n/sum(n) * 100)
    # Sample   cell_status   n   percent
    d_cellstatus$percent <- round(d_cellstatus$percent, digits=1)

    #write.table(df, paste0(outdir, '/malignantCells.sample.xls'), sep='\t', quote=F, col.names=NA)


    p <- ggplot(d_cellstatus, aes(fill=cell_status, y=percent, x=Sample)) +
          geom_bar(stat="identity")+
          geom_text(aes(label = percent), position = position_stack(vjust = 0.5), color="white", size=2) +
          #scale_fill_brewer(palette="Paired")+
          scale_fill_manual(values=c('#6C3483', '#E1C2ED'))+
          theme_minimal() +
          xlab("") + ylab('Percentage of cells (%)') +
          theme(legend.title=element_blank(), 
                legend.key.size = unit(3, 'mm')) +
          theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size=8))
           
    sample_size <- length(unique(df$Sample))
    width <- 3.5
    if (sample_size > 10){ width <- 6 }
    pdf(paste0(outdir, '/malignantCells.sample.pdf'), width = width, height = 3, useDingbats=FALSE)
    print(p)
    dev.off()
}


# barplot - summmary MalignantCellPercentage per cell type2
MalignantCellPercentageForCellType <- function(df, outdir){
    d_cellstatus <- df %>% dplyr::group_by(cell_type2) %>% dplyr::count(cell_status) %>% dplyr::mutate(percent = n/sum(n) * 100)
    # Sample   cell_status   n   percent
    d_cellstatus$percent <- round(d_cellstatus$percent, digits=1)

    #write.table(df, paste0(outdir, '/malignantCells.celltype.xls'), sep='\t', quote=F, col.names=NA)


    p <- ggplot(d_cellstatus, aes(fill=cell_status, y=n, x=reorder(cell_type2, n), n)) +
          geom_bar(stat="identity")+
          coord_flip() +
          geom_text(aes(label = n), position = position_stack(vjust = 0.6), color="black", size=2.5) +
          scale_fill_manual(values=c('#6C3483', '#E1C2ED'))+
          theme_classic() +
          xlab("") + ylab('The number of cells') +
          theme(legend.title=element_blank(), 
                legend.key.size = unit(3, 'mm'),
                legend.position="top") +
          theme(axis.text.x = element_text(size=8))
           

    pdf(paste0(outdir, '/malignantCells.celltype.pdf'), width = 4, height = 2.5, useDingbats=FALSE)
    print(p)
    dev.off()  
}



