
library(Seurat)
library(seumultiome)
library(dplyr)
library(ggplot2)


rds <- 'malignant_cells.sct_lsi.plus2.rds'
finfo <- '21sample.group.txt'
outdir <- 'out_subgroup.method3'

dim_rna <- 30
dim_atac <- 20
seed <- 47



dir.create(outdir)


seu <- readRDS(rds)

info <- read.table(finfo, sep='\t', header=F)
colnames(info) <- c('sample', 'group', 'newid')
info$group[info$group == 'PLAG/L'] <- 'PLAGL'
seu$group <- info$group[match(seu$orig.ident, info$sample)]
Idents(seu) <- 'group'


# extract subgroup
#for (g in info$group){
for (g in c('ZR', 'PF', 'PLAGL')){
    subseu <- NULL
    subseu <- subset(seu, subset=group == g)

    DefaultAssay(subseu) <- 'SCT'
    subseu <- RunPCA(subseu, seed.use=seed, verbose = F)

    original_tfidf <- GetAssayData(seu, assay = "peaks", slot = "data")
    filtered_tfidf <- original_tfidf[, colnames(subseu)]

    # Put it back into your filtered object
    subseu[["peaks"]] <- SetAssayData(
      seurat_filtered[["peaks"]], 
      slot = "data", 
      new.data = filtered_tfidf
    )

    # Then run only SVD
    subseu <- RunSVD(subseu)
    subseu <- harmony::RunHarmony(object=subseu, group.by.vars = 'orig.ident', reduction = 'lsi', assay.use = 'peaks', reduction.save='harmony_atac', project.dim = FALSE)
    
    subseu <- RunWNNPlus(obj=subseu, reduc_rna='pca', reduc_atac='harmony_atac', dim_max_rna=dim_rna, dim_max_atac=dim_atac, res_clus=0.4, seed=seed)

    saveRDS(subseu, file=paste0(outdir,'/', g, '.multiome_integrate.', dim_rna, '_', dim_atac, '.rds'))

    p <- DimPlot(subseu, reduction = 'wnn.umap', group.by = 'Sample', pt.size = 0.1)
    ggsave(paste0(outdir, '/', g, '.multiome_integrate.wnnUMAP.', dim_rna, '_', dim_atac, '.pdf'), width = 6.5, height = 5)

}





