# Hua Sun
# 2025-08-04

library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(seumultiome)
library(GetoptLong)


set.seed(47)
options(future.globals.maxSize=2000000000000)


rds <- 'multiome_integrated.plus.v2.malignantTumorOnly.rds'
finfo <- 'sample.group.txt'
fpfm <- 'pfm.jaspar2022_vertebrates_core.rds'
ref <- 'hg38'
outdir <- 'out_harmony'

regout <- NULL
#regout <- 'percent.mt'

dim_rna <- 30
dim_atac <- 30


GetoptLong(
    "rds=s",       "rds",
    "finfo=s",      "info group",
    "regout=s",    "regress-out",
    "fpfm=s",      "pfm",
    "ref=s",       "ref",
    "outdir=s",    "outdir"
)


dir.create(outdir)
seed <- 47

seu <- readRDS(rds)

info <- read.table(finfo, sep='\t', header=F)
colnames(info) <- c('sample', 'group')
seu$group <- info$group[match(seu$Sample, info$sample)]


# snATAC
seu <- ATACIntegrationHarmonyViaObject(obj=seu, assay_use='peaks', new_reduc='harmony_atac', max_dim=dim_atac, seed=seed)
# saveRDS(seu, 'snatac_renorm.harmony.rds')


# snRNA
# assay = "RNA" default
seu <- SCTransform(seu, vst.flavor='v2', vars.to.regress = regout, variable.features.n=3000)
seu <- RunPCA(seu, seed.use=seed, verbose = F)
# saveRDS(seu, 'sct2Norm_regressout_mt.rds')
# saveRDS(seu, 'sct2Norm.rds')

seu <- harmony::RunHarmony(object=seu, group.by.vars = 'orig.ident', reduction = 'pca', assay.use = 'SCT', reduction.save='harmony_rna')

# wnn
seu <- RunWNNPlus(obj=seu, reduc_rna='harmony_rna', reduc_atac='harmony_atac', dim_max_rna=dim_rna, dim_max_atac=dim_atac, weight_name='SCT.weight', res_clus=0.4, seed=seed)


# save integrated data
saveRDS(seu, file=paste0(outdir,'/multiome_integrate.harmony.', dim_rna, '_', dim_atac, '.rds'))

p <- DimPlot(seu, reduction = 'wnn.umap', group.by = 'group', pt.size = 0.1)
ggsave(paste0(outdir, '/multiome_integrate.harmony.wnnUMAP.', dim_rna, '_', dim_atac, '.pdf'), width = 6.5, height = 5)


# re-call motif 
pfm <- readRDS(fpfm)
genome <- UCSCBSGenome(ref)
seu <- AddGeneActivityAndMotif(obj=seu, genome=genome, pfm=pfm, assay_use='peaks')
saveRDS(seu, file=paste0(outdir,'/multiome_integrate.harmony.plus.', dim_rna, '_', dim_atac, '.rds'))




# extract subgroup
for (g in unique(info$group)){
    #tarcells <- rownames(seu@meta.data)[seu$group == g]
    # subseu <- subset(seu, cells=tarcells)
    samples <- info$sample[info$group == g]
    subseu <- subset(seu, subset = orig.ident %in% samples)
    subseu <- RunWNNPlus(obj=subseu, reduc_rna='harmony_rna', reduc_atac='harmony_atac', dim_max_rna=dim_rna, dim_max_atac=dim_atac, weight_name='SCT.weight', res_clus=0.4, seed=seed)

    saveRDS(subseu, file=paste0(outdir,'/', g, '.multiome_integrate.harmony.', dim_rna, '_', dim_atac, '.rds'))

    p <- DimPlot(subseu, reduction = 'wnn.umap', group.by = 'Sample', pt.size = 0.1)
    ggsave(paste0(outdir, '/', g, '.multiome_integrate.harmony.wnnUMAP.', dim_rna, '_', dim_atac, '.pdf'), width = 6.5, height = 5)
}














