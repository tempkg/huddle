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
    # cell with lineage barcodes (one or multiple lineage barcodes) 
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

d_lb <- data.frame(clean_data %>% separate_rows(lineagebarcode, sep=','))

d_rank <- CalFreqLB(d_lb)
write.table(d_rank, paste0(outdir, '/2.cleanLB.freq.tsv'), sep='\t', quote=F, row.names=F)



split_lb <- strsplit(clean_data$lineagebarcode, ',')
clean_data$rankLB <- sapply(split_lb, function(x){ 
        matched <- d_rank$index_lb[match(x, d_rank$lineagebarcode)]
        paste(matched, collapse = ',')
    })


clean_data$rankLineage <- clean_data$rankLB
clean_data$tempLID[clean_data$comma_count > 0] <- 'Mixed'
clean_data$tempLID <- gsub("LB", "Lineage", clean_data$tempLID)



clean_plus <- clean_data
clean_plus$tempLineage <- clean_plus$tempLID
clean_plus$tempLineage[clean_plus$comma_count > 0] <- NA

d_freq_pure <- as.data.frame(table(clean_plus$tempLineage))
colnames(d_freq_pure) <- c('lineage', 'ncell')
total2 <- sum(d_freq_pure$ncell)
d_freq_pure <- d_freq_pure[order(d_freq_pure$ncell, decreasing = T), ]
row.names(d_freq_pure) <- NULL
d_freq_pure$pureLineage <- paste0('Lineage-', rownames(d_freq_pure))
clean_plus$pureLineage <- d_freq_pure$pureLineage[match(clean_plus$tempLID, d_freq_pure$lineage)]
clean_plus$tempLineage <- NULL




final_meta <- clean_data[,c('cell', 'lineagebarcode', 'comma_count', 'group', 'rankLB')]
final_meta$rankLineage <- clean_plus$pureLineage[match(final_meta$cell, clean_plus$cell)]
final_meta$rankLineage[is.na(final_meta$rankLineage)] <- 'Non-Lineage'

write.table(final_meta, paste0(outdir, '/3.metadata.clean_lineage_cells.tsv'), sep='\t', quote=F, row.names=F)

df <- as.data.frame(table(final_meta$rankLineage))
colnames(df) <- c('lineage', 'ncell')
df <- df[order(df$ncell, decreasing = T), ]


total_n <- sum(df$ncell)
df$ratio <- round(df$ncell/total_n * 100, 2)
write.table(df_pure, paste0(outdir, '/4.ratio.allLineage.out'), sep='\t', quote=F, row.names=F)

df_pure <- df[df$lineage != 'Non-Lineage', ]
total_n <- sum(df_pure$ncell)
df_pure$ratio <- round(df_pure$ncell/total_n * 100, 2)
write.table(df_pure, paste0(outdir, '/5.ratio.pureLineage.out'), sep='\t', quote=F, row.names=F)





#------------ PLOT

PieChart_LineageBarcodeAndCellCounts <- function(df='')
{

    unique_lineages <- unique(df$lineage)
    n <- length(unique_lineages)

    mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(n)
    names(mycolors) <- unique_lineages

    # Manually set "Mixed" to gray
    mycolors["Non-Lineage"] <- "#d8d8d8"
    mycolors["Lineage-1"] <- "#e21f26"

    p <- ggplot(data = df, aes(x = '', y = ncell, fill = lineage)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      theme_void() +
      theme(legend.position = "none") +
      scale_fill_manual(values = mycolors)

    return(p)
}

p <- PieChart_LineageBarcodeAndCellCounts(df)
ggsave(paste0(outdir, '/piechart.lineage_cellCount.all.pdf'), w=0.9, h=0.9)

p <- PieChart_LineageBarcodeAndCellCounts(df_pure)
ggsave(paste0(outdir, '/piechart.Lineage_cellCount.pure.pdf'), w=0.9, h=0.9)






# barplot
Barplot_LineageBarcodeAndCellCounts <- function(df='')
{

    top5_dsub <- head(df, n=5)

    # barplot
    p <- ggplot(data=top5_dsub, aes(y=reorder(lineage, ncell), x=ncell)) +
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

p <- Barplot_LineageBarcodeAndCellCounts(df)
ggsave(paste0(outdir, '/barplot.lineage_cellCount.all.pdf'), w=2, h=1.6)

p <- Barplot_LineageBarcodeAndCellCounts(df_pure)
ggsave(paste0(outdir, '/barplot.lineage_cellCount.pure.pdf'), w=2, h=1.6)







