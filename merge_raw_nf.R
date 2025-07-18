

ft <- 'data.table' # sample  path
out <- 'merged_nf.out'

data <- read.table(ft, sep='\t', header=F)
colnames(data) <- c('sample', 'path')


d_inf <- NULL
for (x in data$sample){
	d <- data[$sample == x, ]

	dnf <- read.table(paste0(d$path, '/nuclear_fraction.tsv'), sep='\t', header=F)
	colnames(dnf) <- c('cell', 'nf')
	dnf$cell <- paste0(d$sample, '_', dnf$cell)
	dnf$sample <- x

	d_inf <- rbind(d_inf, dnf)
}
# cell  nf  sample

write.table(d_inf, out, sep='\t', quote=F)

