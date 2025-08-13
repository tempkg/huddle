
library(ggpubr)
library(stringr)
library(scapekit)
library(ggplot2)
library(dplyr)

outdir <- '.'
dir.create(outdir)


f_cc <- '../out_cellcycle.plus/cellcycle_Score.xls'
f_fus <- '../out_zrFusSig93/metadata_with_featuresSig.xls'


d_cc <- read.table(f_cc, sep='\t', header=T, row.names=1)
d_fus <- read.table(f_fus, sep='\t', header=T, row.names=1)

df <- merge(d_fus[c('seurat_clusters', 'cell_type2', 'fus_score1')], d_cc[c('cc_score1', 'UMAP_1', 'UMAP_2')], by='row.names')

#df <- read.table('cc_fus_signal.table', sep='\t', header=T, row.names=1)
#   seurat_clusters cell_type2  cc_score    fus_score   UMAP_1  UMAP_2
# H_09_2459_AAACAGCCATATTGAC-1    4   RGC-Like    -0.0533687560344054 -0.0598604730682405 1.69493194932509    -4.00301053899123

df$group <- 'YES'
df$group[df$fus_score1 < 0 | df$fus_score1 > 0.2 | (df$cc_score1 < 0.2)] <- 'NO'

write.table(df, paste0(outdir, '/cc_fus_signal.table'), sep='\t', quote=F, row.names=F)


# plot
p <- Scatter2DSigSplit(df=df)
ggsave(paste0(outdir, '/celltypes.cc_fus.pdf'), width=7.5, height=2, useDingbats=FALSE)





