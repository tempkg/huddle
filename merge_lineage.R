library(dplyr)

dir <- 'out_clean_lb_cells'
samples <- c('C08a', 'C08b', 'C08c', 'K10', 'K96a', 'K96b', 'P4', 'T1', 'Z12a', 'Z12b', 'Z12c')

d_info <- NULL
for (name in samples){
	f <- paste0(dir, '/', name, '/3.metadata.clean_lineage_cells.tsv')
	df <- fread(f, data.table=F)
	df$sample <- name
	df$cell <- paste0(name, '_', df$cell)

	d_info <- rbind(d_info, df)
}

d_info <- d_info[, c('sample', 'cell', 'rankLineage')]
colnames(d_info) <- c('sample', 'cell', 'lineage_group')

write.table(d_info, file='merged.clean_lineage_cells.tsv', sep='\t', quote=F, row.names=F)

