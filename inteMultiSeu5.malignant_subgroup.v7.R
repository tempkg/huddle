
library(Seurat)
library(seumultiome)
library(dplyr)
library(ggplot2)


rds <- 'malignant_cells.sct_lsi.plus2.rds'
finfo <- '21sample.group.txt'
outdir <- 'out_subgroup.method7'

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

    DefaultAssay(subseu) <- 'RNA'
    subseu <- SCTransform(seu, vst.flavor = "v2", vars.to.regress = 'percent.mt', variable.features.n = 3000)
    subseu <- RunPCA(subseu, seed.use=seed, verbose = F)
    subseu <- harmony::RunHarmony(object=subseu, group.by.vars = 'orig.ident', reduction = 'pca', assay.use = 'SCT', reduction.save='harmony_rna')

    DefaultAssay(subseu) <- 'peaks'
    subseu <- ATACIntegrationHarmonyViaObject(obj=NULL, assay_use='peaks', group='orig.ident', new_reduc='harmony_atac', min_cutoff=10, max_dim=30, res_clus=0.4, seed=seed)

    subseu <- RunWNNPlus(obj=subseu, reduc_rna='harmony_rna', reduc_atac='harmony_atac', dim_max_rna=dim_rna, dim_max_atac=dim_atac, res_clus=0.4, seed=seed)

    saveRDS(subseu, file=paste0(outdir,'/', g, '.multiome_integrate.', dim_rna, '_', dim_atac, '.rds'))

    p <- DimPlot(subseu, reduction = 'wnn.umap', group.by = 'Sample', pt.size = 0.1)
    ggsave(paste0(outdir, '/', g, '.multiome_integrate.wnnUMAP.', dim_rna, '_', dim_atac, '.pdf'), width = 6.5, height = 5)

}





