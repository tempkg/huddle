# input table
# motif observed    background  percent.observed    percent.background  fold.enrichment pvalue  motif.name  p.adjust


# Set Library
library(GetoptLong)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)


f <- ''
title <- 'Human ST2 Group'
color <- 'red'
prefix <- 'ST2.rankedMotifs'
outdir <- '.'

GetoptLong(
    "f=s",         "rds file",
    "title=s",     "title for plot",
    "color=s",     "color",
    "show_label",  "show label",
    "prefix=s",    "prefix",
    "outdir=s",    "outdir"
)



###############################################
##                    Main
###############################################

dir.create(outdir)


enrichedMotif <- read.table(f, sep='\t', header=T)
enrichedMotif$motif_name_id <- paste0(enrichedMotif$motif.name, ' (', enrichedMotif$motif, ')')
enrichedMotif$probability <- -log10(enrichedMotif$pvalue)
enrichedMotif <- enrichedMotif[order(enrichedMotif$probability, decreasing=TRUE),]
enrichedMotif$rank <- 1:nrow(enrichedMotif)
# add adjusted rank value
enrichedMotif$rank.adjust <- nrow(enrichedMotif):1
enrichedMotif$rank.adjust <- enrichedMotif$rank.adjust/nrow(enrichedMotif)

write.table(enrichedMotif, paste0(outdir, '/', prefix, '.xls'), sep='\t', quote=F, row.names=F)

# show top 10
top100 <- enrichedMotif$motif_name_id[1:100]
top10 <- enrichedMotif$motif_name_id[1:10]

p <- ggplot(enrichedMotif, aes(x = rank.adjust, y = probability)) +
     #geom_point(size=0.2) + 
     geom_line(linewidth=0.2) +
     theme_bw(base_line_size=0.2) +
     #theme_classic(base_line_size=0.3) +
     theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth = 0.2)) +
     labs(title=title) +
     theme(plot.title = element_text(hjust = 0.5, size=7)) +
     theme(text = element_text(size = 6)) +
     labs(x='Ranked motifs', y='Motif probability') +
     geom_point(data = enrichedMotif[enrichedMotif$motif_name_id %in% top100,], size=0.5, color=color)


if (show_label){
    p <- p + geom_text_repel(
            data = subset(enrichedMotif, motif_name_id %in% top10),
            aes(label = motif_name_id),
            size = 3,
            box.padding = unit(0.3, "lines"),
            point.padding = unit(0.3, "lines"),
            max.overlaps = 100,
            min.segment.length = 0,
            segment.angle = 10, segment.size = 0.2
          )
}


pdf(paste0(outdir, '/', prefix, '.pdf'), w=1.8, h=1.5)
print(p)
dev.off()




