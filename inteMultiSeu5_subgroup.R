# Hua Sun
# Seurat v5.1
# v6.5  7/3/25


library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(seumultiome)
library(GetoptLong)


set.seed(47)
options(future.globals.maxSize=2000000000000)


rds <- 'multiome_integrated_plus.seurat5.rds'
finfo <- 'sample.group.txt'  # sample group (no header)
outdir <- 'out_subgroup.harmony'
group <- 'ZR'

seed <- 47
dim_rna <- 30
dim_atac <- 20


GetoptLong(
    "rds=s",       "rds",
    "finfo=s",     "info",
    "group=s",     "group",
    "dim_rna=s",     "dim rna",
    "dim_atac=s",     "dim atac",
    "outdir=s",    "outdir"
)

dir.create(outdir)


seu <- readRDS(rds)

info <- read.table(finfo, sep='\t', header=F)
colnames(info) <- c('sample', 'group')
seu$group <- info$group[match(seu$orig.ident, info$sample)]




RunHarmonyForSubgroup <- function(obj, group, dim_rna, dim_atac, seed, outdir){
    # snRNA
    obj <- harmony::RunHarmony(object=obj, group.by.vars = 'orig.ident', reduction = 'pca', assay.use = 'SCT', reduction.save='harmony_rna')

    # snATAC
    #obj <- ATACIntegrationHarmonyViaObject(obj=obj, assay_use='peaks', new_reduc='harmony_atac', max_dim=dim_atac, seed=seed)
    obj <- harmony::RunHarmony(object=obj, group.by.vars = 'orig.ident', reduction = 'lsi', assay.use = 'peaks', reduction.save='harmony_atac', project.dim = FALSE)

    # wnn
    obj <- RunWNNPlus(obj=obj, reduc_rna='harmony_rna', reduc_atac='harmony_atac', dim_max_rna=dim_rna, dim_max_atac=dim_atac, weight_name='SCT.weight', res_clus=0.4, seed=seed)

    saveRDS(obj, file=paste0(outdir,'/', group, '.multiome_integrate.', dim_rna, '_', dim_atac, '.rds'))

    p <- DimPlot(obj, reduction = 'wnn.umap', group.by = 'seurat_clusters', pt.size = 0.1)
    ggsave(paste0(outdir, '/', group, '.multiome_integrate.wnnUMAP.', dim_rna, '_', dim_atac, '.pdf'), width = 6.5, height = 5)
}




#for (g in unique(info$group)){
#    tarcells <- rownames(seu@meta.data)[seu$group == g]
#    subseu <- subset(seu, cells=tarcells)   
#    RunHarmonyForSubgroup(obj=subseu, group=g, dim_rna=dim_rna, dim_atac=dim_atac, seed=seed, outdir=outdir)
#}

tarcells <- rownames(seu@meta.data)[seu$group == group]
subseu <- subset(seu, cells=tarcells) 
RunHarmonyForSubgroup(obj=subseu, group=group, dim_rna=dim_rna, dim_atac=dim_atac, seed=seed, outdir=outdir)





