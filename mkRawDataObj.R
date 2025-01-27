# Hua Sun


library(Seurat)
library(Signac)
library(dplyr)
library(stringr)
library(stringi)
library(GetoptLong)
library(SeuMultiome)
library(data.table)
library(scapekit)
library(ggplot2)


options(future.globals.maxSize=1000000000000)


rdir <- 'out_cellranger_arc'
ref <- 'mm10'  
min_cells <- 3
width <- 2
height <- 3
outdir <- 'out_rawobject'

GetoptLong(
    "rdir=s",      "dir",
    "ref=s",       "ref",
    "min_cells=i", "min_cells",
    "width=f",     "width",
    "height=f",    "height",
    "outdir=s",    "outdir"
)


dir.create(outdir)


# set
annotation <- Annotations(ref)


# make object list
print('[INFO] make object ...')
multiome_list <- DirCreateMultiomeObject(rdir, ref, annotation, min_cells, outdir)
print(multiome_list)
saveRDS(multiome_list, file=paste0(outdir, '/rawdata.list.rds'))


cutoff_v2 <- CalMedianMaxCutoff(multiome_list)
print(cutoff_v2)
cutoff_x1 <- CalMedianMaxCutoff(multiome_list, .985, .85)
d_cutoff <- data.frame('v2'=cutoff_v2, 'x1'=cutoff_x1)
write.csv(d_cutoff, file = paste0(outdir, '/predicted_cutoff.csv'), row.names = FALSE) 


merged_meta <- NULL
i = 1
for (obj in multiome_list){
    meta <- obj@meta.data
    name <- unique(meta$orig.ident)
    rownames(meta) <- paste0(name, '_', rownames(meta))
    if (i==1){
        merged_meta <- meta
        i = 2
    } else {
        merged_meta <- rbind(merged_meta, meta)
    }
}
write.table(merged_meta, paste0(outdir,'/merged_metadata.rawdata.xls'), sep='\t', quote=F, col.names=NA)


columns <- c('nCount_ATAC', 'nucleosome_signal', 'TSS.enrichment', 'blacklist_ratio', 'pct_reads_in_peaks', 'nCount_RNA', 'nFeature_RNA', 'percent.mt')

# remove one max data
DelMax <- function(df, group, col){
    df_filtered <- df %>% dplyr::group_by(get(group)) %>% dplyr::filter(get(col) < max(get(col)))

    return(df_filtered)
}

# plot
pdf(paste0(outdir,'/rawdata.qc.pdf'), w=width, h=height)
for (x in columns){
    print(x)
    df2 <- DelMax(merged_meta, 'orig.ident', x)
    
    p <- ggplot(df2, aes_string(x='orig.ident', y=x)) + geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.1) + labs(x='') + theme(text = element_text(size = 6)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

    print(p)
    
}
dev.off()


# extra
# Feature QC
merged_meta$gene_group <- ''
merged_meta$gene_group[merged_meta$nFeature_RNA <= 500] <- '<=500'
merged_meta$gene_group[merged_meta$nFeature_RNA > 500 & merged_meta$nFeature_RNA <= 600] <- '500-600'
merged_meta$gene_group[merged_meta$nFeature_RNA > 600 & merged_meta$nFeature_RNA <= 700] <- '600-700'
merged_meta$gene_group[merged_meta$nFeature_RNA > 700 & merged_meta$nFeature_RNA <= 800] <- '700-800'
merged_meta$gene_group[merged_meta$nFeature_RNA > 800 & merged_meta$nFeature_RNA <= 900] <- '800-900'
merged_meta$gene_group[merged_meta$nFeature_RNA > 900 & merged_meta$nFeature_RNA <= 1000] <- '900-1k'
merged_meta$gene_group[merged_meta$nFeature_RNA > 1000 & merged_meta$nFeature_RNA <= 2000] <- '1k-2k'
merged_meta$gene_group[merged_meta$nFeature_RNA > 2000 & merged_meta$nFeature_RNA <= 3000] <- '2k-3k'
merged_meta$gene_group[merged_meta$nFeature_RNA > 3000 & merged_meta$nFeature_RNA <= 4000] <- '3k-4k'
merged_meta$gene_group[merged_meta$nFeature_RNA > 4000 & merged_meta$nFeature_RNA <= 5000] <- '4k-5k'
merged_meta$gene_group[merged_meta$nFeature_RNA > 5000 & merged_meta$nFeature_RNA <= 6000] <- '5k-6k'
merged_meta$gene_group[merged_meta$nFeature_RNA > 6000 & merged_meta$nFeature_RNA <= 7000] <- '6k-7k'
merged_meta$gene_group[merged_meta$nFeature_RNA > 7000 & merged_meta$nFeature_RNA <= 8000] <- '7k-8k'
merged_meta$gene_group[merged_meta$nFeature_RNA > 8000] <- '>8k'

d_meta <- data.frame(table(merged_meta[,c('orig.ident', 'gene_group')])) 
p <- DotLinePlotGroup(data=d_meta, x='gene_group', y='Freq', group='orig.ident', title='nFeature_RNA', y_lab='# Cells', hline='100', hcol='blue')
ggsave(paste0(outdir,'/rawdata.qc.nfeature_rna.pdf'), w=4, h=3)







