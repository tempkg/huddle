# input: signac closest gene result
#tx_id   gene_name   gene_id gene_biotype    type    closest_region  query_region    distance
#ENST00000635509 RP5-857K21.4    ENSG00000230021 lincRNA gap chr1-711923-827670  chr1-826709-827679  0

library(Seurat)
library(Signac)

library(GenomicRanges)
library(S4Vectors)

library(dplyr)
library(stringr)
library(GetoptLong)

library(ggplot2)
library(ggrepel)

rds <- 'celltype_anno.mm_l4/multiome_celltype.rds'
ident <- 'subgroup'
fanno <- ''

label <- ''

outdir <- 'out_openregion'
GetoptLong(
    "rds=s",     ".rds file",
    "ident=s",   "ident",
    "fanno=s",     "closest gene info",
    "label=s",     "label gene",
    "outdir=s",    "outdir"
)






###########################
##          Main
###########################

dir.create(outdir)


seurat_obj <- readRDS(rds)
metadata <- seurat_obj@meta.data


anno_info <- read.table(fanno, sep='\t', header=T)
target_regions <- unique(anno_info$query_region)
tar_region_peak_counts <- as.matrix(seurat_obj@assays$peaks@counts[target_regions,])


i = 1
df <- c()
for (group in unique(metadata$subgroup)){
    sub_cell <- rownames(metadata[metadata$subgroup==group,])
    n_cell <- length(sub_cell)
    sub_mtx <- tar_region_peak_counts[,sub_cell]

    open_cell_ratio <- round(rowSums(sub_mtx>0)/n_cell * 100, 2)
    df_sub <- data.frame('region'=target_regions, 'open_ratio'=open_cell_ratio)

    if (i == 1){
        df <- df_sub
        colnames(df) <- c('region', group)
        i = 2
        } else {
            df[[group]] <- df_sub$open_ratio[match(df$region, df_sub$region)]
        }
}


#df$region <- NULL
df <- merge(df, anno_info, by.x='region', by.y='query_region')
write.table(df, paste0(outdir, '/out_table.xls'), sep='\t', quote=F, row.names=F)


#---- plot

df$ST2_high <- 'No'
df$ST2_high[df$ST2/df$ST1 > 2] <- 'Yes'

write.table(df, paste0(outdir, '/out_table.plus.xls'), sep='\t', quote=F, row.names=F)


p <- ggplot(df, aes(x=ST1, y=ST2)) + geom_point(shape=16, size=1, alpha=.5, color='gray') + 
        geom_point(data=filter(df, ST2_high=='Yes'), shape=16, size=1, alpha=.5, color='#D25827') +  
        theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(text=element_text(size=8), plot.title=element_text(size=9, hjust=0.5)) +
        labs(title='ZR-Fus binding regions', x="ST1 group\n(Ratio of cells in the peaks)", y="ST2 group\n(Ratio of cells in the peaks)") 


# output - fitted-line
p <- p + geom_smooth(method='lm', formula= y~x)
# output - 45C line
# color = "#628cff"
p <- p + geom_abline(slope=1, intercept=0) + lims(x = c(0, 100), y = c(0, 100))


# label
label <- str_split(label, ',')[[1]]
if (length(label) > 1){
    df <- df[order(-df$ST2),]
    top_gene <- head(df$gene_name, n=1)
    label <- c(label, top_gene)

    p <- p + geom_text_repel(
               data = subset(df, gene_name %in% label),
               aes(label = gene_name),
               size = 2, colour = "black", fontface='bold',
               point.padding = 0, # additional padding around each point
               min.segment.length = 0, # draw all line segments
               max.time = 1, max.iter = 1e5,
               max.overlaps = 20, 
               seed=47
               #segment.color = 'grey80'
               ) + geom_point(data = subset(df, gene_name %in% label), shape=1, size=1, lwd=0.05)
}


# width=2.5, height=2.5
pdf(paste0(outdir, '/fus_binding_site.open_ratio.subgroup.st2high.label.pdf'), width=2.3, height=2.3)
print(p)
dev.off()






################################
##  Gene expression
################################

cloest_genes <- anno_info$gene_name
tar_genes <- intersect(cloest_genes, rownames(seurat_obj@assays$SCT@data))

i = 1
df <- c()
for (x in unique(metadata$subgroup)){
    sub_cell <- rownames(metadata[metadata$subgroup==x,])
    n_cell <- length(sub_cell)
    sub_mtx <- seurat_obj@assays$SCT@data[tar_genes, sub_cell]
    df_sub <- data.frame(group=rowMeans(sub_mtx))
    colnames(df_sub) <- x
    print(head(df_sub))
    if (i==1){
        df <- df_sub
        i = 2
    } else {
        df[[x]] <- df_sub[[x]][match(rownames(df), rownames(df_sub))]
    }
}


print(head(df))

max_exp <- max(df)
df$gene_name <- rownames(df)
write.table(df, paste0(outdir, '/out_table.closest_gene.xls'), sep='\t', quote=F, row.names=F)


p <- ggplot(df, aes(x=ST1, y=ST2)) + geom_point(shape=16, size=1, alpha=.5, color='gray') +  
        theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(text=element_text(size=8), plot.title=element_text(size=9, hjust=0.5)) +
        labs(title='Closest genes to ZR-Fus binding regions', x="ST1 group\n(Average gene expression)", y="ST2 group\n(Average gene expression)") +
        lims(x = c(0, max_exp), y = c(0, max_exp))


# output - fitted-line
p <- p + geom_smooth(method='lm', formula= y~x)
# output - 45C line
# color = "#628cff"
p <- p + geom_abline(slope=1, intercept=0)


# label
p <- p + geom_text_repel(
               data = subset(df, gene_name %in% label),
               aes(label = gene_name),
               size = 2, fontface='bold',
               point.padding = 0, # additional padding around each point
               min.segment.length = 0, # draw all line segments
               max.time = 1, max.iter = 1e5,
               max.overlaps = 20, 
               seed=47
               #segment.color = 'grey80'
               ) + geom_point(data = subset(df, gene_name %in% label), shape=1, size=1, lwd=0.05)



# width=2.5, height=2.5
pdf(paste0(outdir, '/fus_binding_peak.cloest_genes.pdf'), width=2.3, height=2.3)
print(p)
dev.off()




