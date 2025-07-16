# Hua Sun

library(data.table)
library(dplyr)

library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)

library(GetoptLong)




fmd <- 'seurat.metaData.xls'
flb <- 'out_lineage/cell_lineageBarcodes.filtered.txt'
outdir <- 'out_clean_lb_cells'


GetoptLong(
    "fmd=s",   "metadata",
    "flb=s",    "barcode file",
    "outdir=s", "outdir"
)

dir.create(outdir)

metadata <- read.table(fmd, sep='\t', header=T, row.names=1)
lineageBarcode <- read.table(flb, sep='\t', header=F)
lineageBarcode$V1 <- paste0(lineageBarcode$V1, '-1')
lineageBarcode <- aggregate(V2 ~ V1, data = lineageBarcode, FUN = function(x) paste(x, collapse = ","))

metadata$lineageBarcode <- lineageBarcode$V2[match(rownames(metadata), lineageBarcode$V1)]


# use with metadata
clean_data <- data.frame('cell'=rownames(metadata), 'lineagebarcode'=metadata$lineageBarcode)
clean_data <- na.omit(clean_data)
write.table(clean_data, file=paste0(outdir, '/1.cleancells_lineagebarcode.tsv'), sep='\t', quote=F, row.names=F)



clean_data$comma_count <- sapply(clean_data$lineagebarcode, function(x){ sum(unlist(strsplit(x, '')) == ',') })
clean_data$group <- 'Lineage'
clean_data$group[clean_data$comma_count > 0] <- 'Mixed'



CalFreqLB <- function(data=NULL){
     
    dfreq <- as.data.frame(table(data$lineagebarcode))
    colnames(dfreq) <- c('lineagebarcode', 'ncell')

    # ordered index
    dfreq <- dfreq[order(dfreq$ncell, decreasing = T),]
    row.names(dfreq) <- NULL
    dfreq$index_lb <- paste0('LB', rownames(dfreq))

    # ratio
    total_n <- sum(dfreq$ncell)
    dfreq$ratio <- round(dfreq$ncell/total_n * 100, 1)

    return(dfreq)
}

# format multiple lineage barcode
d_lb <- data.frame(clean_data %>% separate_rows(lineagebarcode, sep=','))

d_rank <- CalFreqLB(d_lb)


write.table(d_rank, paste0(outdir, '/2.cleanLB.freq.tsv'), sep='\t', quote=F, row.names=F)




split_lb <- strsplit(clean_data$lineagebarcode, ',')
clean_data$rankLB <- sapply(split_lb, function(x){ 
        matched <- d_rank$index_lb[match(x, d_rank$lineagebarcode)]
        paste(matched, collapse = ',')
    })



clean_data$rankLineage <- clean_data$rankLB
clean_data$rankLineage[clean_data$comma_count > 0] <- 'Mixed'
clean_data$rankLineage <- gsub("LB", "Lineage", clean_data$rankLineage)

# lineage ratio
d_freq <- as.data.frame(table(clean_data$rankLineage))
colnames(d_freq) <- c('lineage', 'ncell')
total_n <- nrow(clean_data)
d_freq$ratio <- round(d_freq$ncell/total_n * 100, 2)
d_freq <- d_freq[order(d_freq$ncell, decreasing = T), ]
clean_data$ratioLineage <- d_freq$ratio[match(clean_data$rankLineage, d_freq$lineage)]


write.table(clean_data, paste0(outdir, '/3.summary.lineage_cells.withMixed.tsv'), sep='\t', quote=F, row.names=F)
write.table(d_freq, paste0(outdir, '/3.summary_frequence.withMixed.log'), sep='\t', quote=F, row.names=F)



clean_plus <- clean_data
clean_plus$tempLineage <- clean_plus$rankLineage
clean_plus$tempLineage[clean_plus$comma_count > 0] <- NA



d_freq2 <- as.data.frame(table(clean_plus$tempLineage))
colnames(d_freq2) <- c('lineage', 'ncell')
total2 <- sum(d_freq2$ncell)
d_freq2$ratio <- round(d_freq2$ncell/total2 * 100, 2)
d_freq2 <- d_freq2[order(d_freq2$ncell, decreasing = T), ]
row.names(d_freq2) <- NULL
d_freq2$pureLineage <- paste0('Lineage-', rownames(d_freq2))
clean_plus$pureLineage <- d_freq2$pureLineage[match(clean_plus$rankLineage, d_freq2$lineage)]
clean_plus$pureLineage_ratio <- d_freq2$ratio[match(clean_plus$rankLineage, d_freq2$lineage)]
clean_plus$tempLineage <- NULL



clean_pureOnly <- clean_plus[clean_plus$comma_count == 0,]
write.table(clean_pureOnly, paste0(outdir, '/4.summary.lineage_cells.pure.tsv'), sep='\t', quote=F, row.names=F)
write.table(d_freq2, paste0(outdir, '/4.summary_frequence.pure.log'), sep='\t', quote=F, row.names=F)




PieChart_pureLineageAndCellCounts <- function(data='')
{
    unique_lineages <- unique(data$pureLineage)
    n <- length(unique_lineages)

    mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(n)
    names(mycolors) <- unique_lineages

    mycolors["Lineage1"] <- "#e21f26"

    p <- ggplot(data = data, aes(x = '', y = ncell, fill = pureLineage)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      theme_void() +
      theme(legend.position = "none") +
      scale_fill_manual(values = mycolors)

    return(p)
}
p <- PieChart_pureLineageAndCellCounts(d_freq2)
ggsave(paste0(outdir, '/piechart.pureLineage_cellCount.pdf'), w=0.9, h=0.9)




Barplot_pureLineageAndCellCounts <- function(data='')
{
    top5_dsub <- head(data, n=5)

    # barplot
    p <- ggplot(data=top5_dsub, aes(y=reorder(pureLineage, ncell), x=ncell)) +
        geom_bar(stat="identity", fill='#b25694', width=0.65) +
        geom_text(
        aes(label = ncell), size = 2.5, fontface = "bold", hjust = -0.2) + 
        theme_classic(base_line=0.47) +
        labs(y='', x='The number of cells in the lineage') +
        theme(plot.title=element_text(hjust=0.5, size=8, face = "bold")) +
        theme(text = element_text(face="bold", colour="black", size = 7),
            axis.text.y = element_text(face="bold", colour="black", size = 7),
            axis.text.x = element_text(face="bold", colour="black"),
            axis.line = element_line(colour='black'),
            axis.ticks = element_line(colour='black')
        )

    return(p)
}

p <- Barplot_pureLineageAndCellCounts(d_freq2)
ggsave(paste0(outdir, '/barplot.pureLineage_cellCount.pdf'), w=2, h=1.6)






