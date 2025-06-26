# 2023-08-21 v0.6


# Set Library
library(GetoptLong)
library(Seurat)
library(Signac)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(data.table)
library(this.path)



rds <- ''

idents <- 'cell_type2'
umap <- 'wnn.umap'
split_by <- ''

title <- 'PLAG family (PLAGL1/2,PLAG1)'
motif_id <- 'MA1615.1,MA1548.1,MA0163.1'

outdir <- ''

GetoptLong(
    "func=s",         "analysis item name",
    "rds=s",          ".rds file",
    "title=s",        "title",
    "motif_id=s",     "motif_id",
    "idents=s",       "idents",
    "umap=s",         "reduction",
    "split_by=s",     "split_by",
    "outdir=s",       "output path"
)


#script_dir <- dirname(this.path())
#source(paste0(script_dir,'/src/func.seurat_motif.R'))





##################################
##        Set Func.
##################################

col_group <- c(
        "ST1" = "#EFC41C",                  
        "ST2" = "#D25727",
        "PF" = "#425FAC"
    )
SortedGroupList <- c("ST1", "ST2", "PF")


theme_text_size <- theme(axis.ticks = element_line(size = 0.2)) +
     theme(plot.title = element_text(size = 7, hjust = 0.5),
                text=element_text(size=6), 
                axis.text=element_text(size=6)) +
     theme(axis.ticks = element_line(linewidth = rel(0.5)),
            axis.ticks.length=unit(1, "mm"),
            axis.line = element_line(colour = 'black', size = 0.3)) +
     theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

#theme(legend.position = "none")

theme_legend <- theme(
                legend.title = element_blank(),
                legend.text = element_text(),
                legend.key.size=unit(2,"mm")
            ) 


##################################
##        Main   
##################################

dir.create(outdir)

seurat_obj <- readRDS(rds)


DefaultAssay(seurat_obj) <- 'chromvar'
metadata <- seurat_obj@meta.data


# sum_motif_set & mean_motif_set patterns are same. so, use sum_motif_set

motif_set <- str_split(motif_id,',')[[1]]
tar_motif_data <- as.matrix(seurat_obj@assays$chromvar@data[motif_set,])
d_tar_motif_data <- data.frame(t(tar_motif_data))

# <0 to zero
d_tar_motif_data[d_tar_motif_data <= 0] <- 0

# check zero signal % per group
zero_rows <- d_tar_motif_data[rowSums(d_tar_motif_data[])==0,]
zero_rows$subgroup <- metadata$subgroup[match(rownames(zero_rows), rownames(metadata))]
d_zero_ratio <- data.frame(table(zero_rows$subgroup))
d_total_cells <- data.frame(table(metadata$subgroup))
d_zero_ratio$total_cells <- d_total_cells$Freq[match(d_zero_ratio$Var1,d_total_cells$Var1)]
d_zero_ratio$zero_ratio <- d_zero_ratio$Freq/d_zero_ratio$total_cells
d_zero_ratio

# only use >0 signal cells
d_tar_motif_data <- d_tar_motif_data[rowSums(d_tar_motif_data[])>0,]
#d_tar_motif_data$subgroup <- metadata$subgroup[match(rownames(d_tar_motif_data), rownames(metadata))]




# I) sum PLAG family per cell
report <- data.frame(rowSums(d_tar_motif_data))
#report <- data.frame(rowMeans(d_tar_motif_data))
colnames(report) <- 'sum_motif_sig'
report$subgroup <- metadata$subgroup[match(rownames(report), rownames(metadata))]
report$subgroup <- factor(report$subgroup, levels=SortedGroupList)

# mean is better than median in this data
d_sum_family <- report %>%
        group_by(subgroup) %>%
        summarise_at(vars(sum_motif_sig), list(name = mean))
d_sum_family <- data.frame(d_sum_family)
# bar plot
d_sum_family$subgroup <- factor(d_sum_family$subgroup, levels=c('PF', 'ST1', 'ST2'))
p <- ggplot(d_sum_family, aes(x=name, y=subgroup, fill=subgroup)) + 
    geom_bar(stat='identity') +
    theme_classic(base_line_size=0.2) +
    labs(title='PLAG family', x='Motif activity\n(Mean)', y='') +
    scale_fill_manual(values=col_group) +
    theme(legend.position = "none") +
    theme_text_size

outfile <- paste0(outdir, '/plag_family.motifSig.subgroup.barplot.pdf')
pdf(outfile, width = 1, height = 1.2, useDingbats=FALSE)
print(p)
dev.off()




# II) sum each motif per subgroup
d_tar_motif_data$subgroup <- metadata$subgroup[match(rownames(d_tar_motif_data), rownames(metadata))]

avg_motif_set <- d_tar_motif_data %>%
                group_by(subgroup) %>%
                summarise_at(motif_set, list(motif = mean))
avg_motif_set <- data.frame(avg_motif_set)
colnames(avg_motif_set) <- c('subgroup', motif_set)
rownames(avg_motif_set) <- avg_motif_set[,1]
avg_motif_set[,1] <- NULL
colnames(avg_motif_set) <- c('PLAGL1', 'PLAGL2', 'PLAG1')


library(ComplexHeatmap)

NR=dim(avg_motif_set)[1]
NC=dim(avg_motif_set)[2]

ht <- Heatmap(as.matrix(avg_motif_set),     
                col = c('#3498DB','#F9E79F','#E74C3C'),
                show_row_dend = F,
                show_column_dend = F,
                row_names_side = "left",
                column_names_side = "top",

                # heatmap cell size
                width = unit(5, "mm")*NC,
                height = unit(4, "mm")*NR,

                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 5),

                # legend
                heatmap_legend_param = list(
                        title = 'Motif activity (Mean)',
                        #title_position = "leftcenter-rot",
                        direction = "horizontal",
                        title_gp = gpar(fontsize = 6), 
                        labels_gp = gpar(fontsize = 5),
                        grid_width = unit(1.5, "mm"),
                        grid_height = unit(1.5, "mm")
                        )
            )
                    

outfile <- paste0(outdir, '/plag_family.heatmap.motifSig.subgroup.pdf')
pdf(outfile, width = 1.5, height = 1.5, useDingbats=FALSE)
draw(ht, heatmap_legend_side = "bottom")
dev.off()






VlnPlot_Test <- function(){
    # sum_motif_set
    p <- ggviolin(report, x='subgroup', y='sum_motif_sig', fill = "subgroup",
                palette = col_group, size=0.1) +
        geom_boxplot(width=0.2, size=0.2, outlier.size = 0.2) +
            theme_classic(base_line_size=0.2) +
            labs(title=title, x="", y='Motif activity') +
            theme_text_size + theme(legend.position = "none") +
        geom_hline(yintercept=0, linetype="dashed", color = "blue", linewidth=0.3)
        
    out_name <- str_replace_all(motif_id, ',', '_')
    pdf(paste0(outdir, '/', out_name, '.vln.pdf'), w=2.5, h=2)
    print(p)
    dev.off() 
}





