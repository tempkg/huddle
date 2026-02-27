
# Hua Sun
# 2026-2-26 v5


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
min_cnv_gene <- dyml$parameter$min_cnv_gene

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


d_merged <- ''
i <- 0
for (s in unique(df_cnv_gene$sample)){

    # per sample
    ds <- df_cnv_gene[df_cnv_gene$sample==s,]

    # only use metadata barcodes
    ds_report <- metadata_filtered[metadata_filtered$sample==s, c('sample', 'seurat_clusters')]
    cells <- row.names(ds_report)

    ds_exp <- d_icnv_exp_filtered[cells,]


    # calculate
    total_genes <- ncol(ds_exp)
    print(total_genes)

    gain <- rowSums(ds_exp >= cutoff_gain)
    loss <- rowSums(ds_exp <= cutoff_loss)
    cnv_gene_count <- gain + loss

    gain_ratio <- round(gain/total_genes * 100, 2)
    loss_ratio <- round(loss/total_genes * 100, 2)

    # global cnv-ratio
    total_cnv_genes <- gain + loss
    cnv_ratio <- round(total_cnv_genes/total_genes * 100, 2)

    d_cnv_ratio <- data.frame(cnv_gene_count, cnv_ratio, gain_ratio, loss_ratio)
    rownames(d_cnv_ratio) <- row.names(ds_exp)

    # define malignant
    d_cnv_ratio$cell_status <- 'Non-malignant'
    d_cnv_ratio$cell_status[d_cnv_ratio$cnv_gene_count >= min_cnv_gene & d_cnv_ratio$cnv_ratio >= min_ratio] <- 'Malignant'


    i <- i + 1
    if (i==1) {
        d_merged <- d_cnv_ratio
    } else {
        d_merged <- rbind(d_merged, d_cnv_ratio)
    }

}


print(dim(d_merged))
# regions, n_gene_count, cnv_regions, no need them
d_merged <- d_merged %>% select(sample, seurat_clusters, cnv_gene_count, cnv_ratio, gain_ratio, loss_ratio, cell_status)
write.table(d_merged, paste0(outdir, '/infercnv.estimate_malignant_cells.global.xls'), sep='\t', quote=FALSE, col.names=NA)



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
p <- MalignantCellPercentagePerSample(d_merged)
sample_size <- length(unique(d_merged$sample))
width <- 3.5
if (sample_size > 10){ width <- 6 }
ggsave(paste0(outdir, '/malignantCells.sample.pdf'), width = width, height = 3)




#-------- [Extra] 7. clean reference cells
ref_clu <- stringr::str_split(control, ',')[[1]]
d_poor_ref <- d_merged %>% 
            filter(seurat_clusters %in% ref_clu) %>% 
            filter(cell_status=='Malignant')

poor_ref_cells <- d_poor_ref[,1]
writeLines(poor_ref_cells, paste0(outdir, '/poor_ref_cells.out'))




