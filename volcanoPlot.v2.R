# Hua Sun

library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(scapekit)


f <- 'diffExp_results.pval.xls'
out <- 'diffExp_results.pval.pdf'

diffExp <- fread(f, data.table=F)


diffExp$signal <- 'NotDiff' 
diffExp$signal[diffExp$log2FoldChange > 1 & diffExp$pvalue < 0.05] <- 'Up' 
diffExp$signal[diffExp$log2FoldChange < -1 & diffExp$pvalue < 0.05] <- 'Down' 

top_gene <- ShowTopNGenes(df=df, n=10, logfc='log2FoldChange', group='signal', gene='V1')



p <- VolcanoPlot(df=diffExp, 
	x='log2FoldChange', 
	y='pvalue', 
	group='signal', 
	gene_col='V1', 
	y_lab='-log10(p-val)',
	x_cutoff=1,
	pt_size=1,
	title='Relapse vs Primary (EPN)',
	label=top_gene
)


ggsave(out, width = 4, height = 3)



