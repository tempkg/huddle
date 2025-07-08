# Hua Sun

# fpeak (no header)
# chr20-24744751-24746693   PF


library(Seurat)
library(Signac)
library(data.table)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(RColorBrewer)
library(randomcoloR)
library(circlize)
library(this.path)
library(GetoptLong)
library(ggplot2)



func <- ''
rds <- 'celltype_hs.l5.renamed/sc_celltype_anno.rds'

sample_info <- ''
cell_info <- ''
fcell <- ''

fpeak <- ''

outdir <- 'out_peaks'

GetoptLong(
    "func=s",          "function",
    "rds=s",           ".rds seurat file",
    "sample_info=s",   "sample info",
    "cell_info=s",     "cell info",
    "fcell=s",         "cell list",
    "fpeak=s",         "file from FindMarker output or peaks",
    "outdir=s",        "output path"
)




#################################
##           Set Func.
#################################

# Find closest genes for peaks
FindClosestFeatureForPeaks <- function(obj=NULL, peak_set=NULL, outfile=NULL){
    print('[INFO] find closest genes ...')
    closest_genes <- ClosestFeature(obj, regions = peak_set)
    # save peak matched gene name
    write.table(closest_genes, file = outfile, sep = "\t", row.names = F)

    return(closest_genes)
}



# only use for check since output pdf not work well
ShowSeuratHeatmap_Group <- function(obj=NULL, assay='peaks', features=NULL, outfile='out.pdf'){
    col_group <- c(
            "ZR" = "#E21F26",                  
            "EWP" = "#F57F20",
            "PF" = "#7C50A0"
        )

    p <- DoHeatmap(object=obj, 
        features=features, 
        group.by='subgroup', assay=assay, slot='data', angle=0,
        size=4, disp.min=0, disp.max=2, group.bar.height = 0.05, 
        group.colors=col_group ) + 
    guides(color="none") +
    guides(fill = guide_colourbar(barwidth = 0.8, barheight = 5)) +
    scale_fill_gradientn(colours = c("#010087", "#FCF3CF", "red"))

    #pdf(outfile, width = 6, height = 3)
    ggplot2::ggsave( plot=p, filename=outfile, width = 10, height = 8)
    
}



# complexHeatmap
HeatMap_Peaks_Groups <- function(obj=NULL, peaks=NULL, outfile='heatmap.peak.pdf')
{
    d_peak <- obj$peaks@data[peaks,]
    metadata <- obj@meta.data

    col_group <- c(
        "ZR" = "#E21F26",                  
        "EWP" = "#F57F20",
        "PF" = "#7C50A0"
    )
    SortedGroupList <- c("ZR", "EWP", "PF")

    top_anno <- HeatmapAnnotation(
                Group = anno_simple(x = metadata$subgroup, simple_anno_size = unit(2, "mm"), col=col_group),
                annotation_name_side = "left",
                annotation_name_gp = gpar(fontsize = 6)
        )


    color_codes <- c('#2E317E','#FCF3CF','#EA2626')
    color_heatmap <- colorRamp2(seq(0, max(d_peak), length = length(color_codes)), color_codes)
    
    NR=dim(d_peak)[1]
    NC=dim(d_peak)[2]

    ht <- Heatmap(
                as.matrix(d_peak),
                    
                #col = color_heatmap,
                    
                show_row_names = F,
                show_row_dend = F,
                cluster_rows = T,

                show_column_names = F,
                cluster_columns = T,
                
                # heatmap cell size
                #width = unit(4, "mm")*NC,
                #height = unit(0.05, "mm")*NR,

                
                row_title = paste0('Peaks'), 
                row_title_gp = gpar(fontsize = 8),

                # split
                #column_title = NULL,
                column_title_gp = gpar(fontsize = 8),
                column_split = factor(metadata$subgroup, levels=SortedGroupList), 
                cluster_column_slices = F,
                column_gap = unit(0.5, "mm"), 
                
                # top-anno
                top_annotation = top_anno,
                
                # legend
                heatmap_legend_param = list(
                        title = "Mean motif activity score",
                        #title_position = "leftcenter-rot",
                        direction = "horizontal",
                        title_gp = gpar(fontsize = 7), 
                        labels_gp = gpar(fontsize = 6),
                        grid_width = unit(2.5, "mm"),
                        grid_height = unit(2.5, "mm")
                        )
            )

    pdf(outfile, width = 3, height = 3, useDingbats=FALSE)
    draw(ht, heatmap_legend_side = "bottom")
    dev.off()
}




## Heatmap - z-score
HeatMap_Peaks_Groups_rowZscore <- function(data=NULL, metadata=NULL, 
    row_title='Peaks', outfile='heatmap.peak.pdf')
{   
    col_group <- c(
        "ZR" = "#E21F26",                  
        "EWP" = "#F57F20",
        "PF" = "#7C50A0"
    )
    SortedGroupList <- c("ZR", "EWP", "PF")

    top_anno <- HeatmapAnnotation(
                Group = anno_simple(x = metadata$subgroup, simple_anno_size = unit(2, "mm"), col=col_group),
                annotation_name_side = "left",
                annotation_name_gp = gpar(fontsize = 6)
        )


    # z-score: row z-score
    scaled_mat = t(scale(t(data)))
    #col_zscore = colorRamp2(c(-2, -1, 0, 1, 2), c("#09103b", "#5f79cf", "white", "#eb6565", "#540506"))
    col_zscore = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
    title_ht = "Row Z-Score"

    ht <- Heatmap(
                as.matrix(scaled_mat),
                    
                col = col_zscore,
                    
                show_row_names = F,
                #row_names_gp = gpar(fontsize = 6),
                show_row_dend = F,
                cluster_rows = T,
                row_dend_reorder = FALSE,

                show_column_names = F,
                show_column_dend = F,
                cluster_columns = T,
                
                # heatmap cell size
                #width = unit(4, "mm")*NC,
                #height = unit(0.05, "mm")*NR,

                row_title = row_title, 
                row_title_gp = gpar(fontsize = 7),

                # split
                #column_title = NULL,
                column_title_gp = gpar(fontsize = 8),
                column_split = factor(metadata$subgroup, levels=SortedGroupList), 
                cluster_column_slices = F,
                column_gap = unit(0.5, "mm"),
                
                # split rows to three groups
                #row_split = row_split,

                # top-anno
                top_annotation = top_anno,

                # legend
                heatmap_legend_param = list(
                        title = title_ht,
                        #title_position = "leftcenter-rot",
                        direction = "horizontal",
                        title_gp = gpar(fontsize = 6), 
                        labels_gp = gpar(fontsize = 5),
                        grid_width = unit(2.5, "mm"),
                        grid_height = unit(2, "mm")
                        )
            )

    pdf(outfile, width = 2, height = 3, useDingbats=FALSE)
    draw(ht, heatmap_legend_side = "bottom")
    dev.off()
}






## Heatmap - z-score
HeatMap_Peaks_Groups_rowZscore_Plus <- function(data=NULL, metadata=NULL, row_info=NULL,
    row_title='Peaks', outfile='heatmap.peak.pdf')
{   
    col_group <- c(
        "ZR" = "#E21F26",                  
        "EWP" = "#F57F20",
        "PF" = "#7C50A0"
    )
    SortedGroupList <- c("ZR", "EWP", "PF")

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


    n = length(unique(metadata$Sample))
    set.seed(2000 + n)
    palette <- distinctColorPalette(n)
    names(palette) <- sort(unique(metadata$Sample))
    
    top_anno <- HeatmapAnnotation(
                Group = anno_simple(x = metadata$subgroup, simple_anno_size = unit(2, "mm"), col=col_group),
                Sample = anno_simple(x = metadata$Sample, simple_anno_size = unit(2, "mm"), col=col_sample),
                annotation_name_side = "left",
                annotation_name_gp = gpar(fontsize = 6)
        )


    # z-score: row z-score
    scaled_mat = t(scale(t(data)))
    #col_zscore = colorRamp2(c(-2, -1, 0, 1, 2), c("#09103b", "#5f79cf", "white", "#eb6565", "#540506"))
    col_zscore = colorRamp2(c(-2, 0, 2), c("#2E86C1", "white", "#CB4335"))
    title_ht = "Row Z-Score"

    ht <- Heatmap(
                as.matrix(scaled_mat),
                    
                col = col_zscore,
                    
                show_row_names = F,
                #row_names_gp = gpar(fontsize = 6),
                show_row_dend = F,
                cluster_rows = T,
                row_dend_reorder = FALSE,

                show_column_names = F,
                show_column_dend = F,
                cluster_columns = T,
                
                # heatmap cell size
                #width = unit(4, "mm")*NC,
                #height = unit(0.05, "mm")*NR,


                # split
                #column_title = NULL,
                column_title_gp = gpar(fontsize = 8),
                column_split = factor(metadata$subgroup, levels=SortedGroupList), 
                cluster_column_slices = F,
                column_gap = unit(0.5, "mm"),
                

                row_title = row_title, 
                row_title_gp = gpar(fontsize = 7),

                # split rows to three groups
                row_split = factor(row_info$group, levels=SortedGroupList),
                cluster_row_slices = F,
                row_gap = unit(0.5, "mm"),

                # top-anno
                top_annotation = top_anno,

                # legend
                heatmap_legend_param = list(
                        title = title_ht,
                        #title_position = "leftcenter-rot",
                        title_position = "lefttop",
                        direction = "horizontal",
                        title_gp = gpar(fontsize = 6), 
                        labels_gp = gpar(fontsize = 5),
                        grid_width = unit(2.5, "mm"),
                        grid_height = unit(2, "mm")
                        )
            )

    pdf(outfile, width = 2, height = 3, useDingbats=FALSE)
    draw(ht, heatmap_legend_side = "bottom")
    dev.off()
}





#################################
##            Main
#################################


##----------- prepare -----------##

dir.create(outdir)


# 1.read object
print('[INFO] reading .rds ...')
seurat_obj <- readRDS(rds)


# 2.sample info
if (nchar(sample_info)>1){
    sample_info <- read.table(sample_info, sep='\t', header=T)
    colnames(sample_info) <- c('sample', 'group', 'subgroup')

    # add subgroup
    seurat_obj@meta.data$subgroup <- sample_info$subgroup[match(seurat_obj@meta.data$Sample, sample_info$sample)]
}


# 3.cell info for malignant cells
if (nchar(cell_info)>1){
    cell_info <- read.table(cell_info, sep='\t', header=T, row.names=1)
    seurat_obj@meta.data$cell_status <- cell_info$cell_status[match(rownames(seurat_obj@meta.data), rownames(cell_info))]

    # only use malignant cells
    seurat_obj <- subset(x = seurat_obj, subset = cell_status == 'Malignant')
}


# 4.extract target cells
if (nchar(fcell)>1){
    barcodes <- readLines(fcell)
    print(length(barcodes))
    # only use malignant cells
    seurat_obj <- subset(x = seurat_obj, cells = barcodes)
    print(nrow(seurat_obj@meta.data))
}

# only 500 cells
Idents(seurat_obj) <- 'subgroup'
seurat_obj <- subset(x=seurat_obj, downsample=500)

DefaultAssay(seurat_obj) <- 'peaks'
metadata <- seurat_obj@meta.data


# read & filter dup. peaks
peak_info <- read.table(fpeak, sep='\t', header=F)
colnames(peak_info) <- c('peak', 'group')
# filter duplicated peak
print(nrow(peak_info))

p_freq <- data.frame(table(peak_info$peak))
filtered_peaks <- as.vector(p_freq[p_freq$Freq==1,]$Var1)
peak_info <- peak_info %>% filter(peak %in% filtered_peaks)
peak_info <- peak_info[order(peak_info$group),]
write.table(peak_info, paste0(outdir, '/filteredDup.peaks.out'), sep='\t', quote=F, row.names=F)
print(nrow(peak_info))




##----------- Peaks heatmap -----------## 
# Heatmap for the peaks per group
if (func=='heatmap-peak'){ 
    d_matrix <- as.matrix(seurat_obj$peaks@data[peak_info$peak,])
    #HeatMap_Peaks_Groups_rowZscore(data=d_matrix, metadata=metadata, row_title='Peaks', outfile=paste0(outdir,'/heatmap.peaks.group.zscore.plus.pdf'))
    HeatMap_Peaks_Groups_rowZscore_Plus(data=d_matrix, metadata=metadata, row_info=peak_info, row_title='Top 200 peaks per group', outfile=paste0(outdir,'/heatmap.peaks.group.zscore.plus.pdf'))
}




##----------- Peaks + Gene Exp -----------## 
# Heatmap for the peaks per group
if (func=='peak-exp'){ 
    
    d_closest_genes <- FindClosestFeatureForPeaks(seurat_obj, as.vector(peak_info$peak), paste0(outdir, '/peaks.closestGenes.xls'))
    d_closest_genes$group <- peak_info$group[match(d_closest_genes$query_region, peak_info$peak)]

    # gene list
    print('Unique genes in the closest genes')
    closest_genes <- as.vector(unique(d_closest_genes$gene_name))
    print(length(closest_genes))

    all_genes <- rownames(seurat_obj$SCT@data)
    # exp >= 10 cells with count > 3
    #cutoff_cell <- as.integer(ncol(seurat_obj$SCT@data)*0.1)
    exp_genes <- rownames(seurat_obj$SCT@data[rowSums(as.matrix(seurat_obj$SCT@counts)>3)>9,])
    closest_genes_exp <- intersect(closest_genes, exp_genes)
    print(length(closest_genes_exp))
    write.table(d_closest_genes[d_closest_genes$gene_name %in% closest_genes_exp,], file = paste0(outdir, '/peaks.closestGenes.withExpGenes.xls'), sep = "\t", row.names = F)

    d_matrix <- as.matrix(seurat_obj$SCT@data[closest_genes_exp,])
    gene_info <- d_closest_genes %>% filter(gene_name %in% closest_genes_exp) %>% select(gene_name, group)
    gene_info <- gene_info[!duplicated(gene_info$gene_name),]

    # best for z-score
    #HeatMap_Peaks_Groups_rowZscore(data=d_matrix, metadata=metadata, row_title='Closest gene expression', outfile=paste0(outdir,'/heatmap.peaks_closestGenes.zscore.pdf'))
    HeatMap_Peaks_Groups_rowZscore_Plus(data=d_matrix, metadata=metadata, row_info=gene_info, row_title='Closest gene expression', outfile=paste0(outdir,'/heatmap.peaks_closestGenes.zscore.plus.pdf'))
}





