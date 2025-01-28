
library(Seurat)
library(seuextra)
library(scapekit)
library(ggplot2)


seu <-readRDS('multiome_integrated.rds')

d_motif <- ExtractMotifSigTranspose(seu, 'MA1615.1,MA1548.1', 'cell_type2')
motifs <- stringr::str_split('MA1615.1,MA1548.1', ',')[[1]]
print(motifs)
d <- CalSignalCountRatioByGroup(df=d_motif, columns=motifs, group='cell_type2', subgroup='orig.ident')
write.table(d, 'motif.plgl.cellcount.tsv', sep='\t', quote=F, row.names=F)

d <- CalCellRatioForMotifSig(obj=seu, motifs='MA1615.1,MA1548.1', group='cell_type2', sample_group='orig.ident')
write.table(d, 'motif.plgl.cellratio.tsv', sep='\t', quote=F, row.names=F)



#d$celltype <- factor(d$celltype, level=c('RGC','CycProg','OPC','Neuron','Microglia','EC'))

p <- BarPlotSplitGroup(data=d, split_group='cell_type2', x='orig.ident', y='MA1615.1', y_lab='Plagl1 motif ratio (%)')
ggsave('motif.plgl.cellratio.pdf', w=4, h=4)








