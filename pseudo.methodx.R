
library(Seurat)
library(pseudotime)
library(slingshot)
library(scapekit)
library(ggplot2)
library(dplyr)
library(data.table)


rds <- 'scrna_sct2_harmony.Lineage-1.rds'
fmeta <- 'out_celltype/cluster_cellType.xls'

outdir <- 'out_phate'
dir.create(outdir)

seu <- readRDS(rds)
meta <- fread(fmeta, data.table=F)
seu$cell_type2 <- meta$cell_type2[match(seu$seurat_clusters, meta$seurat_clusters)]

sub_obj <- subset(seu, subset= cell_type2 != 'Unclassified')




sub_obj <- RunPhate(seu_obj=sub_obj, use_py='/Users/hua/opt/miniconda3/envs/phate/bin/python')
saveRDS(sub_obj, paste0(outdir, '/seurat_with_phate.rds'))

meta <- sub_obj@meta.data



start_root <- 'CycProg'
RunSlingshotSeurat(seu=sub_obj, assay='SCT', ident='cell_type2', reduction='PHATE', start_root=start_root,  save=TRUE, outdir=outdir)

# plot pseudotime
fpt <- paste0(outdir, '/2.slingshot.pseudotime.xls')
pt <- fread(fpt, data.table=F)

pt$cell_type2 <- meta$cell_type2[match(pt$V1, rownames(meta))]
write.table(pt, paste0(outdir, '/pseudotime.celltype2.xls'), sep='\t', quote=F, row.names=F)

p <- DensityPlot_Pseudotime(data=pt, x='Lineage1', t1=2, t2=10, xlim=c(-1,15))
ggsave(paste0(outdir, '/pseudotime_density.lineage1.pdf'), w=1.5, h=1.2)





