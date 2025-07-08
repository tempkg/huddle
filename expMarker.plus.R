# Hua Sun



library(Seurat)
library(data.table)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(RColorBrewer)
library(randomcoloR)
library(circlize)
library(GetoptLong)



rds <- 'sc_celltype_anno.rds'
fmarker <- ''
tf <- '~/Documents/database/cistarget/tf_lists/allTFs_hg38.txt'
group_by <- 'subgroup'
width <- 3.5
height <- 1.5
style <- '1'
outdir <- 'out_expMarker'

GetoptLong(
    "rds=s",        ".rds seurat file",
    "fmarker=s",    "file from FindMarker output or peaks",
    "tf=s",         "TF data",
    "style=s",      "style",
    "group_by=s",   "group by",
    "width=f",      "width",
    "height=f",     "height",
    "outdir=s",     "output path"
)




#################################
##           Set Func.
#################################

color_group <- function()
{
    col_group <- c(
            "ST1" = "#EFC41C",                  
            "ST2" = "#D25727",
            "PF" = "#425FAC"
        )

    col_group
}


color_sample <- function()
{
    col_sample <- c(
        "H_02_2138" = "#DBA97F",
        "H_03_1276" = "#7FC4DA",
        "H_04_1478" = "#C6E251",
        "H_05_0665" = "#DA979E",
        "H_06_1225" = "#759B80",
        "H_09_2459" = "#8869D6",
        "H_09_2461" = "#E35FD9",
        "H_15_4174" = "#78E18E",
        "H_16_05327" = "#6CE652",
        "H_17_07683" = "#D8E1DA",
        "H_93_0050" = "#DE4C92",
        "H_93_1095" = "#CFE4A4",
        "H_94_0046" = "#C9B5DC",
        "H_96_0032" = "#6996D7",
        "H_96_0119" = "#A53CE4",
        "H_96_933" = "#E15A55",
        "H_97_0095" = "#7EE7D0",
        "H_98_0130" = "#CC86C4",
        "H_99_0398" = "#DFBD55"
    )

    col_sample
}


HeatMap_ExpGroups_Vertical <- function(
    expMat=NULL, metadata=NULL, metagroup='subgroup', 
    gene_info=NULL, show_gene=NULL, width = 3.5, height = 1.5,
    outfile='heatmap.pdf')
{   
    # z-score: row z-score
    scaled_mat = t(scale(t(expMat)))
    dim(scaled_mat)

    SortedGroupList <- c("ST1", "ST2", "PF")
    
    if (length(unique(gene_info$cluster))==2){
        col_group <- c(                 
            "ST" = "#D25727",
            "PF" = "#425FAC"
        )
        SortedGroupList <- c("ST", "PF")
    }


    n = length(unique(metadata$Sample))
    set.seed(2000 + n)
    palette <- distinctColorPalette(n)
    names(palette) <- sort(unique(metadata$Sample))
    
    top_anno <- HeatmapAnnotation(
                #Group = anno_simple(x = metadata$subgroup, simple_anno_size = unit(2, "mm"), col=col_group),
                Sample = anno_simple(x = metadata$Sample, simple_anno_size = unit(2, "mm"), col=color_sample),
                annotation_name_side = "right",
                annotation_name_gp = gpar(fontsize = 6)
        )

    # gene info
    gene_at <- gene_info$index[gene_info$gene %in% show_gene]
    gene_anno <- rowAnnotation( link = anno_mark(at=gene_at, labels=show_gene, which="bottom", link_width=unit(2,"mm"), labels_gp=gpar(fontsize = 5)) )

    #col_zscore = colorRamp2(c(-2, -1, 0, 1, 2), c("#09103b", "#5f79cf", "white", "#eb6565", "#540506"))
    col_zscore = colorRamp2(c(-2, 0, 2), c("#2E86C1", "white", "#CB4335"))
    title_ht = "Row Z-Score"

    ht <- Heatmap(
                as.matrix(scaled_mat),
                    
                col = col_zscore,
                    
                show_row_names = F,
                show_row_dend = F,
                cluster_rows = T,

                show_column_names = F,
                show_column_dend = F,
                cluster_columns = T,
                
                # heatmap cell size
                #width = unit(4, "mm")*NC,
                #height = unit(0.05, "mm")*NR,


                # split
                #column_title = NULL,
                column_title_gp = gpar(fontsize = 7),
                column_split = factor(metadata[[metagroup]], levels=SortedGroupList), 
                cluster_column_slices = F,
                column_gap = unit(0.5, "mm"),
                
                
                # split rows to three groups
                row_title = 'Differentially expressed genes', 
                row_title_gp = gpar(fontsize = 7),
                row_split = factor(gene_info$cluster, levels=SortedGroupList),
                cluster_row_slices = F,
                row_gap = unit(0.5, "mm"),

                # top-anno
                top_annotation = top_anno,
                right_annotation = gene_anno,

                # legend
                heatmap_legend_param = list(
                        title = 'Expression\nZ-score',
                        title_gp = gpar(fontsize = 5), 
                        labels_gp = gpar(fontsize = 5),
                        grid_width = unit(2, "mm"),
                        grid_height = unit(2, "mm")
                        )
            )

    pdf(outfile, width=width, height=height, useDingbats=FALSE)
    draw(ht)
    dev.off()
}



## Heatmap - horizontal
HeatMap_ExpGroups_Horizontal_2 <- function(
    expMat=NULL, metadata=NULL, metagroup='subgroup', 
    gene_info=NULL, show_gene=NULL, width = 3.5, height = 1.5,
    outfile='heatmap.pdf')
{   
    # z-score: row z-score
    scaled_mat = scale(t(expMat))
    dim(scaled_mat)

    SortedGroupList <- c("ST1", "ST2", "PF")

    if (length(unique(gene_info$cluster))==2){
        col_group <- c(                 
            "ST" = "#D25727",
            "PF" = "#425FAC"
        )
        SortedGroupList <- c("ST", "PF")
    }
    
    n = length(unique(metadata$Sample))
    set.seed(2000 + n)
    palette <- distinctColorPalette(n)
    names(palette) <- sort(unique(metadata$Sample))
    #print(palette)
    # rowAnnotation can't set color
    row_ha = rowAnnotation(df=metadata$Sample, simple_anno_size=unit(2,"mm"), show_legend=FALSE, show_annotation_name=FALSE)

    col_zscore = colorRamp2(c(-2, 0, 2), c("#3A7EB8", "white", "#FF0000"))


    # gene info
    gene_at <- gene_info$index[gene_info$gene %in% show_gene]
    gene_anno <- HeatmapAnnotation( link = anno_mark(at=gene_at, labels=show_gene, which="bottom", link_width=unit(2,"mm"), labels_gp=gpar(fontsize = 5)) )

    ht <- Heatmap(
                as.matrix(scaled_mat),
                    
                col = col_zscore,
                    
                show_row_names = F,
                show_row_dend = F,
                cluster_rows = T,
                row_dend_reorder = T,

                show_column_names = F,
                show_column_dend = F,
                cluster_columns = T,
                
                # split
                column_title = NULL,
                column_names_gp = gpar(fontsize = 6),
                column_split = factor(gene_info$cluster, levels=SortedGroupList), 
                cluster_column_slices = F,
                column_gap = unit(0.5, "mm"),
                

                # split rows to three groups
                row_title_gp = gpar(fontsize = 7),
                row_names_gp = gpar(fontsize = 6),
                row_split = factor(metadata[[metagroup]], levels=SortedGroupList),
                cluster_row_slices = F,
                row_gap = unit(0.5, "mm"),

                # show text
                top_annotation = gene_anno,
                left_annotation = row_ha,

                # legend
                heatmap_legend_param = list(
                        title = 'Expression\nZ-score',
                        title_gp = gpar(fontsize = 5), 
                        labels_gp = gpar(fontsize = 5),
                        grid_width = unit(2, "mm"),
                        grid_height = unit(2, "mm")
                        )
            )

    pdf(outfile, width=width, height=height, useDingbats=FALSE)
    draw(ht)
    dev.off()
}





## Heatmap - horizontal
HeatMap_ExpGroups_Horizontal_3 <- function(
    expMat=NULL, metadata=NULL, metagroup='subgroup', 
    gene_info=NULL, show_gene=NULL, width = 3.5, height = 1.5,
    outfile='heatmap.pdf')
{   

    # z-score: row z-score
    scaled_mat = scale(t(expMat))
    dim(scaled_mat)

    col_group <- c(
        "ST1" = "#EFC41C",                  
        "ST2" = "#D25727",
        "PF" = "#425FAC"
    )
    SortedGroupList <- c("ST1", "ST2", "PF")

    if (length(unique(gene_info$cluster))==2){
        col_group <- c(                 
            "ST" = "#D25727",
            "PF" = "#425FAC"
        )
        SortedGroupList <- c("ST", "PF")
    }
    
    # gene info
    gene_at <- gene_info$index[gene_info$gene %in% show_gene]
    gene_anno <- HeatmapAnnotation( link = anno_mark(at=gene_at, labels=show_gene, which="bottom", link_width=unit(2,"mm"), labels_gp=gpar(fontsize = 5)) )


    n = length(unique(metadata$Sample))
    set.seed(2000 + n)
    palette <- distinctColorPalette(n)
    names(palette) <- unique(metadata$Sample)
    # rowAnnotation can't set color
    row_ha = rowAnnotation(df=metadata$Sample, simple_anno_size=unit(2,"mm"), show_legend=FALSE, show_annotation_name=FALSE, col=color_sample)

    col_zscore = colorRamp2(c(-2, 0, 2), c("#3A7EB8", "white", "#FF0000"))

    ht <- Heatmap(
                as.matrix(scaled_mat),
                    
                col = col_zscore,
                    
                show_row_names = F,
                show_row_dend = F,
                cluster_rows = T,
                #row_dend_reorder = FALSE,

                show_column_names = F,
                show_column_dend = F,
                cluster_columns = T,
                
                # split
                column_title = NULL,
                column_names_gp = gpar(fontsize = 6),
                column_split = factor(gene_info$cluster, levels=SortedGroupList), 
                cluster_column_slices = F,
                column_gap = unit(0.5, "mm"),
                

                # split rows to three groups
                row_title_gp = gpar(fontsize = 7),
                row_names_gp = gpar(fontsize = 6),
                row_split = factor(metadata[[metagroup]], levels=SortedGroupList),
                cluster_row_slices = F,
                row_gap = unit(0.5, "mm"),

                # show text
                top_annotation = gene_anno,
                left_annotation = row_ha,

                # legend
                heatmap_legend_param = list(
                        title = 'Expression Z-score',
                        title_position = "lefttop",
                        direction = "horizontal",
                        title_gp = gpar(fontsize = 5), 
                        labels_gp = gpar(fontsize = 5),
                        grid_width = unit(2, "mm"),
                        grid_height = unit(2, "mm")
                        )
            )

    pdf(outfile, width=width, height=height, useDingbats=FALSE)
    draw(ht, heatmap_legend_side = "bottom")
    dev.off()
}




#################################
##            Main
#################################

dir.create(outdir)


# 1.read object
print('[INFO] reading .rds ...')
seurat_obj <- readRDS(rds)


d_exp_markers <- read.table(fmarker, sep='\t', header=T)
# "p_val"  "avg_log2FC" "pct.1"  "pct.2"  "p_val_adj"  "cluster"  "gene"  "diff_pct"  

# read tf list
tf_list <- readLines(tf)
# add TF to marker
d_exp_markers$TF <- tf_list[match(d_exp_markers$gene, tf_list)]

d_exp_marker_filter <- d_exp_markers
# filter markers
#d_exp_marker_filter <- d_exp_markers %>% 
#            filter(p_val_adj < 0.05) %>%
#            filter(avg_log2FC > 0.585) %>% 
#            filter(pct.1 > 0.3) %>% 
#            filter(pct.2 < 0.15)
#write.table(d_exp_marker_filter, paste0(outdir, '/expMarker.filter.xls'), sep='\t', quote=F, row.names=F)




# target genes for exp matrix

# label genes
# 'RELA','ZFTA','C11orf95'
gene_info <- d_exp_marker_filter %>% select(gene, cluster)
gene_info$index <- 1:nrow(gene_info)


seurat_obj$group <- seurat_obj$subgroup
seurat_obj$group[seurat_obj$group!='PF'] <- 'ST'
Idents(seurat_obj) <- group_by
seurat_obj_sub <- subset(x=seurat_obj, downsample=500)

metadata <- seurat_obj_sub@meta.data

d_matrix <- as.matrix(seurat_obj_sub$SCT@data[gene_info$gene,])

#show_gene <- intersect(gene_info$gene, c(tf_list, 'C11orf95'))
show_gene <- intersect(gene_info$gene, tf_list)

# best for z-score
if (style == '1'){
    outfile <- paste0(outdir,'/heatmap.expMarker.filtered.1.pdf')
    HeatMap_ExpGroups_Vertical(expMat=d_matrix, metadata=metadata, metagroup=group_by, 
        gene_info=gene_info, show_gene=show_gene, width = width, height = height, outfile=outfile)
}
if (style == '2'){
    outfile <- paste0(outdir,'/heatmap.expMarker.filtered.2.pdf')
    HeatMap_ExpGroups_Horizontal_2(expMat=d_matrix, metadata=metadata, metagroup=group_by, 
        gene_info=gene_info, show_gene=show_gene, width = width, height = height, outfile=outfile)
}
if (style == '3'){
    outfile <- paste0(outdir,'/heatmap.expMarker.filtered.3.pdf')
    HeatMap_ExpGroups_Horizontal_3(expMat=d_matrix, metadata=metadata, metagroup=group_by, 
            gene_info=gene_info, show_gene=show_gene, width = width, height = height, outfile=outfile)
}







