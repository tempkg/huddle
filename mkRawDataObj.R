# Hua Sun


library(Seurat)
library(Signac)
library(dplyr)
library(stringr)
library(stringi)
library(GetoptLong)
library(SeuMultiome)
library(data.table)
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
saveRDS(multiome_list, file=paste0(outdir,'/rawdata.list.rds'))

merged_meta <- c()
i <- 1
for (obj in multiome_list){
    meta <- obj@meta.data
    name <- unique(meta$orig.ident)
    meta <- paste0(name, '_', rownames(meta))
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
    df_filtered <- df %>% group_by(get(group)) %>% filter(get(col) < max(get(col)))

    return(df_filtered)
}

# plot
pdf(paste0(outdir,'/rawdata.qc.pdf'), w=width, h=height)
for (x in columns){
    print(x)
    df2 <- DelMax(merged_meta, 'orig.ident', x)
    
    p <- ggplot(df2, aes_string(x='orig.ident', y=x)) + geom_boxplot(outlier.size = 0.1) + labs(x='') + theme(text = element_text(size = 6)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

    print(p)
    
}
dev.off()








