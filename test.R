
library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(dplyr)
library(data.table)

set.seed(42)

outdir <- 'out_monocle3'
dir.create(outdir)

rds <- 'seurat5.1_v6.2/multiome_integrated.plus.rds'
seu <- readRDS(rds)

fmeta <- 'out_celltype/metaData.cellType.xls'
meta <- fread(fmeta, data.table=F)
seu$cell_type2 <- meta$cell_type2[match(rownames(seu@meta.data), meta$V1)]


# https://stuartlab.org/signac/articles/monocle


Monocle3MakeCDS <- function(
    obj=NULL, 
    subgroup=NULL,
    assay='SCT', 
    ident='cell_type2', 
    umap='wnn.umap'
){
    DefaultAssay(obj) <- assay
    Idents(obj) <- ident

    if (length(subgroup) > 0){
      obj <- subset(obj, idents = subgroup)
    }

    obj[["UMAP"]] <- obj[[umap]]

    # make cds
    cds <- as.cell_data_set(obj)
    cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
    cds <- learn_graph(cds, use_partition = TRUE)
    
    return(cds)
}



# make cds
#cds <- Monocle3MakeCDS( obj=seu, assay='SCT', ident='seurat_clusters', subgroup=c(8,6,3,19,1,13,7,9,18,10), umap='wnn.umap')
#cds <- Monocle3MakeCDS( obj=seu, assay='SCT', ident='seurat_clusters', subgroup=c(8,11,4,5,12,15,14,2,0), umap='wnn.umap')
#cds <- Monocle3MakeCDS( obj=seu, assay='SCT', ident='cell_type2', subgroup=c('RGC', 'CycProg', 'Neuron'), umap='wnn.umap')
cds <- Monocle3MakeCDS( obj=seu, assay='SCT', ident='cell_type2', umap='wnn.umap')
saveRDS(cds, paste0(outdir, '/monocle3_cds.rds'))


# choose cells for a subset ?
# cds_subset <- choose_cells(cds)



# method 1: set root cells
#meta_filtered <- meta %>% filter(orig.ident=='E12' & cell_type2=='RGC')
#root_cells <- rownames(meta_filtered)
#cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)


# method 2: set root cell from map to call pseudotime
cds <- order_cells(cds)

#saveRDS(cds, 'monocle3_cds.ordered.rds')


# 1. export pseudotime
pseudotime <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
d_pseudotime <- data.frame(pseudotime)
write.table(d_pseudotime, paste0(outdir, '/pseudotime.out'), sep='\t', quote=F, col.names=NA)


# 2. plot trajectories colored by pseudotime
plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

p <- plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_leaves=FALSE,
  label_branch_points=FALSE
)
pdf(paste0(outdir, '/pseudotime.umap.pdf'), w=6.5, h=5)
p
dev.off()



# plot cell type
plot_cells(cds,
           color_cells_by = "cell_type2",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)




# Perform differential expression analysis
# With regression:
#gene_fits <- fit_models(cds, model_formula_str = "~orig.ident")
#fit_coefs <- coefficient_table(gene_fits)
#emb_time_terms <- fit_coefs %>% filter(term == "orig.ident")
#emb_time_terms <- emb_time_terms %>% mutate(q_value = p.adjust(p_value))
#sig_genes <- emb_time_terms %>% filter (q_value < 0.05) %>% pull(gene_short_name)

# With graph autocorrelation:
#pr_test_res <- graph_test(cds,  neighbor_graph="principal_graph", cores=4)
#pr_deg_ids <- row.names(subset(pr_test_res, q_value < 0.05))




# 3. heatmap plot for pseudotime and gene/cluster
# https://github.com/cole-trapnell-lab/monocle-release/issues/295
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

# Identifying genes that vary in expression over a trajectory
#https://cole-trapnell-lab.github.io/monocle-release/monocle3/
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
# In Monocle3, "Morans_I" refers to a statistical measure called Moran's I, which is used to assess the spatial autocorrelation of gene expression along a cell trajectory, essentially indicating how similar the expression levels of a gene are between neighboring cells on the trajectory; a high Moran's I value suggests that nearby cells tend to have similar expression levels for that gene
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.5))

# used pseudotime
pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes


# predict the best k-mean cutoff
# https://rpkgs.datanovia.com/factoextra/reference/fviz_nbclust.html
#gap_stat <- cluster::clusGap(pt.matrix, FUN = kmeans, nstart = 25, K.max = 10, B = 50)
#factoextra::fviz_gap_stat(gap_stat)




# plot heatmap - k-mean based (*****)
# set km=6 as 6 groups
km <- 4
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = circlize::colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = km,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)


#print(htkm)

#---- Extract htkm cluster matched genes
HM <- draw(htkm)  #Show the heatmap

r.dend <- row_dend(HM)  #If needed, extract row dendrogram
rcl.list <- row_order(HM)  #Extract clusters (output is a list)
  
lapply(rcl.list, function(x) length(x))  #check/confirm size gene clusters

library(magrittr) # needed to load the pipe function '%%'
 
clu_df <- lapply(names(rcl.list), function(i){
            out <- data.frame(GeneID = rownames(pt.matrix[rcl.list[[i]],]),
                     Cluster = paste0("cluster", i),
                     stringsAsFactors = FALSE)
            return(out)
        }) %>% do.call(rbind, .)

#export
write.table(clu_df, file= paste0(outdir, '/gene_clusters.kmean_heatmap.txt'), sep="\t", quote=F, row.names=FALSE)

# save plot
pdf(paste0(outdir, '/gene_clusters.kmean_heatmap.pdf'), w=5, h=6)
HM
dev.off()







# plot Ward.D2 Hierarchical Clustering (****)
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = circlize::colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

print(hthc)









