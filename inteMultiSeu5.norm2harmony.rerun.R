# Hua Sun


library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(seumultiome)
library(GetoptLong)



options(future.globals.maxSize=2000000000000)


rna <- 'filtered.remDoublet.rds'
atac <- 'snatac_peaks.list.rds'
fpfm <- 'pfm.jaspar2022_vertebrates_core.rds'
ref <- 'hg38'
outdir <- 'out_re-integration'


dim_rna <- 30
dim_atac <- 30


GetoptLong(
    "rna=s",       "rna",
    "atac=s",       "atac",
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


snrna_list <- readRDS(rna)
snatac_list <- readRDS(atac)


# Seurat 5
integ_rna <- Seurat5SCTIntegration(obj_list=snrna_list, max_dim=dim_rna, regress='percent.mt', integ_method='harmony', new_reduc='harmony_rna', seed=seed)
saveRDS(integ_rna, file=paste0(outdir,'/snRNA.harmony.regout_mt.', dim_rna, '.rds'))


integ_atac <- ATACIntegrationHarmony(snatac_list=snatac_list, max_dim=dim_atac, new_reduc='harmony_atac', seed=seed)
saveRDS(integ_rna, file=paste0(outdir,'/snATAC.harmony.', dim_atac, '.rds'))


multiome_wnn <- RunWNN_Customized(integ_rna=integ_rna, integ_atac=integ_atac, reduc_rna='harmony_rna', reduc_atac='harmony_atac', res_clus=0.6, seed=seed)


# save integrated data
saveRDS(multiome_wnn, file=paste0(outdir,'/multiome_harmony.', dim_rna, '_', dim_atac, '.rds'))

p <- DimPlot(multiome_wnn, reduction = 'wnn.umap', group.by = 'group', pt.size = 0.1)
ggsave(paste0(outdir, '/multiome_harmony.wnnUMAP.', dim_rna, '_', dim_atac, '.pdf'), width = 6.5, height = 5)


# re-call motif 
pfm <- readRDS(fpfm)
genome <- UCSCBSGenome(ref)
integrated_multiome <- AddGeneActivityAndMotif(obj=integrated_multiome, genome=genome, pfm=pfm, assay_use='peaks')

# Seurat v5
integrated_multiome[['RNA']] <- JoinLayers(integrated_multiome[['RNA']])
saveRDS(integrated_multiome, file=paste0(outdir,'/multiome_harmony.', dim_rna, '_', dim_atac, '.plus2.rds'))










