
library(stringr)
library(Seurat)
library(Signac)
library(seuextra)
library(scapekit)
library(dplyr)

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(data.table)
library(GenomicRanges)


rds <- 'multiome_integrated.plus.rds'
meta <- 'out_celltype/metaData.cellType.xls'
y <- 'ZR'
x <- 'EWP'
fbed <- '/db/ZRFUS.2021.cd/enhancer_peaks.mm10/fusionBindingSite.M21_HA_HOMER_peaks_merged_DIFF_IgG.simple.bed'
cutoff <- '/db/ZRFUS.2021.cd/2021.cd.93genes/93genes/ZRFUS1_signature_genes.93.mouse.gene'
outdir <- 'out_zrfus_binding.zr_ewp'

dir.create(outdir)



FusionBindingPeaks <- function(
      seurat_obj=NULL, 
      metadata=NULL, 
      fdb_bed=NULL, 
      cutoff=0.3
){
        # cell type comparison
      barcode.y_x <- row.names(metadata[metadata$cell_type2 %in% c(y, x),])
      barcode.y <- row.names(metadata[metadata$cell_type2==y,])
      barcode.x <- row.names(metadata[metadata$cell_type2==x,])

      # overlap peaks
      peaks <- row.names(seurat_obj@assays$peaks@data)
      overlap_peaks <- CalIntersectPeaks(peaks=peaks, fbed=fdb_bed, cutoff=cutoff)
      writeLines(overlap_peaks, paste0(outdir, '/overlap.zrfus_binding_peaks.out'))

      tar_counts <- as.matrix(seurat_obj@assays$peaks@counts[overlap_peaks, barcode.y_x])

      # based on ZR-Fus genes
      filtered_peaks <- row.names(tar_counts[rowSums(tar_counts>1)>9, ])
      print(length(filtered_peaks))

      writeLines(filtered_peaks, paste0(outdir, '/overlap_filtered.zrfus_binding_peaks.out'))


      # calculate peak signal by average
      tar_peakdata_fus_bind_data <- as.matrix(seurat_obj@assays$peaks@data[filtered_peaks, barcode.y_x])

      mean_peak_sig.y = apply(tar_peakdata_fus_bind_data[,barcode.y], 1, mean, na.rm=TRUE)
      mean_peak_sig.x = apply(tar_peakdata_fus_bind_data[,barcode.x], 1, mean, na.rm=TRUE)

      df <- data.frame('peaksig_y'=mean_peak_sig.y, 'peaksig_x'=mean_peak_sig.x)
      df <- df[rowSums(df[])>0,]

      return(df)
}



PlotTwoGroups <- function(
      df=NULL,
      x=NULL, 
      y=NULL, 
      max_sig=1.3,
      point_size=0.5,
      outdir='.'
){
      # plot 2d - scatter plot
      p <- ggplot(df, aes(x=peaksig_x, y=peaksig_y) ) +
            #geom_point(shape = 21, size=point_size, stroke=.5, color = "white", fill = "#EC5228", alpha=.8) + 
            geom_point(shape = 16, size=point_size, color = "#EC5228", alpha=.6) + 
            xlim(0, max_sig) + ylim(0, max_sig) +
            geom_abline(intercept = 0, slope = 1, linewidth = 0.3, linetype = "dashed") +
            theme_linedraw() + 
            theme(axis.ticks = element_line(size = 0.5)) +
            theme(panel.border = element_rect(size = 1)) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            labs(title='Accessibility of ZFTA-RELA\n binding sites', x=paste0(x,' (Mean peak signal)'), y=paste0(y,' (Mean peak signal)')) +
            theme(plot.title = element_text(size = 7, face='bold', hjust = 0.5)) +
            theme(text = element_text(size = 6, face='bold')) 
            
      # standard size
      ggsave(paste0(outdir, '/comp.zrfus_binding_peaksig.pdf'), width = 1.8, height = 1.8)



      # box plot
      df$group_x <- x
      df$group_y <- y
      df_group <- df[,c('peaksig_x', 'group_x')]
      colnames(df_group) <- c('peak_signal', 'group')
      temp <- df[,c('peaksig_y', 'group_y')]
      colnames(temp) <- c('peak_signal', 'group')
      df_group <- rbind(df_group, temp)

      p <- ggboxplot(df_group, x = "group", y = "peak_signal",
                color = "group", palette = "jco") +
            stat_compare_means(size = 3) +
            labs(x='', y='Peak signal') +
            theme(plot.title = element_text(size = 9, hjust = 0.5)) +
            theme(text = element_text(size = 9)) +
            theme(legend.position = "none")

      ggsave(paste0(outdir, '/comp.zrfus_binding_peaksig.boxplot.pdf'), width = 2, height = 2.5)
}


#------------------------- Main

# read obj
seu <- readRDS(rds)

# subset 
info <- read.table(finfo, sep='\t', header=F)
colnames(info) <- c('sample', 'group')
seu$group <- info$group[match(seu$Sample, info$sample)]

# subgroup
subseu <- subset(x=seu, subset = group %in% c('ZR', 'EWP'))

peaks <- row.names(subseu@assays$peaks@data)
overlap_peaks <- CalIntersectPeaks(peaks=peaks, fbed=fdb_bed, cutoff=cutoff)
writeLines(overlap_peaks, paste0(outdir, '/overlap.zrfus_binding_peaks.out'))



# closet genes
DefaultAssay(subseu) <- 'peaks'
closet_genes_peaks <- Signac::ClosestFeature(subseu, rownames(df))
write.table(closet_genes_peaks, paste0(outdir, '/comp.zrfus_binding_peaksig.closetgene.tsv'), sep='\t', quote=F, row.names=F, col.names=T)




# check overlap with 93 ZR-Fus genes
#zrfus_genes <- readLines(dyml$zrgenes)
#overlap_closet_genes_peaks <- closet_genes_peaks %>% filter(gene_name %in% zrfus_genes)
#write.table(overlap_closet_genes_peaks, paste0(outdir, '/closetgene.overlap_with_93zrfus_genes.tsv'), sep='\t', quote=F, col.names=T)




ScatterPlotWithLabel <- function(
      df=NULL,
      x=NULL, 
      y=NULL, 
      query_genes=NA,
      show_genes=NA,
      max_sig=1.3,
      point_size=0.5
){
      # plot 2d - scatter plot
      p <- ggplot(df, aes(x=peaksig_x, y=peaksig_y) ) +
            #geom_point(shape = 21, size=point_size, stroke=.5, color = "white", fill = "#EC5228", alpha=.8) + 
            geom_point(shape = 16, size=point_size, color = "#cdcdcd", alpha=.6) + 
            xlim(0, max_sig) + ylim(0, max_sig) +
            geom_abline(intercept = 0, slope = 1, linewidth = 0.3, linetype = "dashed") +
            theme_linedraw() + 
            theme(axis.ticks = element_line(size = 0.5)) +
            theme(panel.border = element_rect(size = 1)) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            labs(title='Accessibility of ZFTA-RELA\n binding sites', x=paste0(x,' (Mean peak signal)'), y=paste0(y,' (Mean peak signal)')) +
            theme(plot.title = element_text(size = 7, face='bold', hjust = 0.5)) +
            theme(text = element_text(size = 6, face='bold')) 

      p <- p + geom_point(data=subset(df, gene_name %in% query_genes), 
                  shape = 16, size=point_size, color = "#EC5228")

      p <- p + geom_text_repel(data = subset(df, gene_name %in% show_genes),
                     aes(label = gene_name),
                     size = 2,
                     colour = "black", fontface='bold',
                     min.segment.length = 0,
                     segment.size  = 0.3,
                     box.padding = unit(0.2, "lines")
                  )
      
      p 
}

df$gene_name <- closet_genes_peaks$gene_name[match(rownames(df), closet_genes_peaks$query_region)]

#features <- c('Aif1l','Axl','Jag1','Gli2','Lfng','Rela')
#features <- c('Aif1l','Ephb2','Jag1','Gli2','Rela')

#p <- ScatterPlotWithLabel(df=df, x=x, y=y, query_genes=zrfus_genes, show_genes=features, max_sig=1.3, point_size=0.5)
# standard size
#ggsave(paste0(outdir, '/comp.zrfus_binding_peaksig.withFusGenes.pdf'), width = 1.8, height = 1.8)







