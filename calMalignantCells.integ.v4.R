
# Hua Sun
# 2025-11-19 v4.1


library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(this.path)
library(yaml)

#path <- dirname(this.path())
#source(paste0(path, '/plot.R'))

args <- commandArgs(trailingOnly = TRUE)
print(args)
f <- args[1]
dyml <- yaml.load_file(f)

fmd <- dyml$input$meta
control <- dyml$input$ctl

cutoff_region_size <- dyml$parameter$min_region
min_ratio <- dyml$parameter$min_ratio
min_cnv <- dyml$parameter$min_cnv
min_gene <- dyml$parameter$min_gene

cutoff_gain <- dyml$parameter$cutoff_gain
cutoff_loss <- dyml$parameter$cutoff_loss

rm_clu <- dyml$parameter$rm_clu
rm_chr <- dyml$parameter$rm_chr

icnv_exp <- dyml$data$icnv_exp
hmm_gene <- dyml$data$hmm_gene
hmm_region <- dyml$data$hmm_region

outdir <- dyml$out$outdir

dir.create(outdir)



#-------- 1. select cnv regions and genes
# read hmm-cnv-region file
df_cnv_region <- read.table(hmm_region, sep='\t', header=T)
# cell_group_name   cnv_name    state   chr start   end

# filter region
# calculate region unit as 'Mb'
df_cnv_region$region_size <- (df_cnv_region$end - df_cnv_region$start)/1000000

# cutoff region size >= x Mb

df_cnv_region <- df_cnv_region %>% filter(region_size >= cutoff_region_size)
filtered_cnv_region <- unique(df_cnv_region[['cnv_name']])
#print(filtered_cnv_region)

# read hmm-cnv-gene file
df_cnv_gene <- read.table(hmm_gene, sep='\t', header=T)
# cell_group_name   gene_region_name    state   gene    chr start   end

df_cnv_gene <- df_cnv_gene %>% filter(gene_region_name %in% filtered_cnv_region)
print(dim(df_cnv_gene))

# filter chr, if set 'rm_chr'
if (nchar(rm_chr) > 0){
    remove_chr <- str_split(rm_chr, ',')[[1]]
    df_cnv_gene <- df_cnv_gene %>% filter(!(chr %in% remove_chr))
}
print(unique(df_cnv_gene$chr))



#-------- 2. read infercnv exp-obj & cnv genes
filtered_cnv_gene_set <- unique(df_cnv_gene$gene)

d_icnv_exp <- readRDS(icnv_exp)
# cell(rowname) gene1 gene2 gene3 ...
used_cells <- rownames(d_icnv_exp)

d_icnv_exp_filtered <- d_icnv_exp[, filtered_cnv_gene_set]
print(dim(d_icnv_exp_filtered))




#-------- 3. calculate malignant/non-malignant per sample cell based on cnv region
# sample name
df_cnv_gene$sample <- gsub("\\.\\w+", "", df_cnv_gene$cell_group_name, perl=TRUE)




#-------- 4. metadata from seurat
# read metadata
metadata <- read.table(fmd, sep='\t', header=T)
rownames(metadata) <- metadata$cell
metadata$cell <- NULL
print(dim(metadata))
# used cells
metadata_filtered <- metadata[used_cells,]

# filter clusters
if (nchar(rm_clu) > 0){
    del_clusters <- str_split(rm_clu, ',')[[1]]
    metadata_filtered <- metadata_filtered %>% filter(!(seurat_clusters %in% del_clusters))
}

print(unique(metadata_filtered$seurat_clusters))
print(nrow(metadata_filtered))


df_merged <- ''
i <- 0
for (s in unique(df_cnv_gene$sample)){
    ds <- df_cnv_gene[df_cnv_gene$sample==s,]

    # only use metadata barcodes
    ds_report <- metadata_filtered[metadata_filtered$sample==s, c('sample', 'seurat_clusters')]
    cells <- row.names(ds_report)

    s_regions <- c()
    s_cnv_gene_count <- c()
    s_cnv_ratio <- c()

    # region cnv
    for (r in unique(ds$gene_region_name)){

        genes <- ds$gene[ds$gene_region_name==r]
        
        n_gene_region <- length(genes)
        # skip regions if less than 20 genes
        if (n_gene_region < min_gene){
            print(paste('skip ', r, '-', n_gene_region))
            next
        }
        #print(paste(r, '-', n_gene_region))

        ds_exp <- d_icnv_exp_filtered[cells, genes]
        total_genes <- length(genes)
        gain <- rowSums(ds_exp >= cutoff_gain)
        loss <- rowSums(ds_exp <= cutoff_loss)

        s_regions <- cbind(s_regions, r)

        # region cnv gene counts
        cnv_genes <- gain + loss
        s_cnv_gene_count <- cbind(s_cnv_gene_count, cnv_genes)
        colnames(s_cnv_gene_count) <- as.vector(s_regions)

        # region cnv ratio
        cnv_ratio <- (gain+loss)/total_genes
        s_cnv_ratio <- cbind(s_cnv_ratio, cnv_ratio)
        colnames(s_cnv_ratio) <- as.vector(s_regions)
    }
    
    #ds_report$regions <- toString(as.vector(s_regions))
    ds_report$total_cnv.regions <- length(as.vector(s_regions))
    ds_report$cnv_gene_count <- as.vector(unite(col="merged_cnv_genes", as.data.frame(s_cnv_gene_count), sep=","))
    ds_report$cnv_ratio <- as.vector(unite(col="merged_cnv_ratio", as.data.frame(s_cnv_ratio), sep=","))
    
    ds_report$cnv_gene_count.mean <- round(apply(s_cnv_gene_count, 1, mean, na.rm=TRUE), 0)
    ds_report$cnv_gene_count.max <- apply(s_cnv_gene_count, 1, max, na.rm=TRUE)
    ds_report$cnv_ratio.mean <- apply(s_cnv_ratio, 1, mean, na.rm=TRUE)
    ds_report$cnv_ratio.max <- apply(s_cnv_ratio, 1, max, na.rm=TRUE)

    ds_report$cnv_ratio.max <- round(ds_report$cnv_ratio.max, 4)
    ds_report$cell_status <- 'Non-malignant'
    ds_report$cell_status[ds_report$cnv_ratio.max >= min_ratio | ds_report$cnv_gene_count.max >= min_cnv] <- 'Malignant'

    i <- i + 1
    if (i==1) {
        df_merged <- ds_report
    } else {
        df_merged <- rbind(df_merged, ds_report)
    }
    
}


print(dim(df_merged))
df_merged <- df_merged %>% select(sample, seurat_clusters, total_cnv.regions, cnv_gene_count.mean, cnv_gene_count.max, cnv_ratio.mean, cnv_ratio.max, cell_status)



#-------- 5. plot
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


# malignant cell percentage per sample
p <- MalignantCellPercentagePerSample(df_merged)
sample_size <- length(unique(df_merged$sample))
width <- 3.5
if (sample_size > 10){ width <- 6 }
ggsave(paste0(outdir, '/malignantCells.sample.pdf'), width = width, height = 3)




#-------- [Extra] 6. calculate malignant/non-malignant per cell based on all gene values (simple approach)

# cell g1 g2 g3 ...
total_genes <- ncol(d_icnv_exp)
print(total_genes)

gain <- rowSums(d_icnv_exp >= cutoff_gain)
loss <- rowSums(d_icnv_exp <= cutoff_loss)
total_cnv_genes <- gain + loss

gain_ratio.by_allgenes <- round(gain/total_genes * 100, 2)
loss_ratio.by_allgenes <- round(loss/total_genes * 100, 2)
cnv_ratio.by_allgenes <- round(total_cnv_genes/total_genes * 100, 2)

d_cnv_ratio2 <- data.frame(gain_ratio.by_allgenes, loss_ratio.by_allgenes, cnv_ratio.by_allgenes)
rownames(d_cnv_ratio2) <- row.names(d_icnv_exp)

df_merged <- merge(df_merged, d_cnv_ratio2, by='row.names', all.x=TRUE)

# output
write.table(df_merged, paste0(outdir, '/infercnv.estimate_malignant_cells.', cutoff_region_size, 'Mb.xls'), sep='\t', quote=FALSE, row.names=FALSE)





#-------- [Extra] 7. clean reference cells
ref_clu <- stringr::str_split(control, ',')[[1]]
d_poor_ref <- df_merged %>% 
            filter(seurat_clusters %in% ref_clu) %>% 
            filter(cell_status=='Malignant')

poor_ref_cells <- d_poor_ref[,1]
writeLines(poor_ref_cells, paste0(outdir, '/poor_ref_cells.out'))




