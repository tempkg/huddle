

ft <- 'groups/hs.multiome.22.out' # sample  path
out <- 'hs22.merged_nf.47.out'

data <- read.table(ft, sep='\t', header=F)
colnames(data) <- c('sample', 'path')


d_inf <- NULL
for (x in data$sample){
	print(x)
	d <- data[data$sample == x, ]

	file <- paste0(d$path, '/nuclear_fraction.tsv')
	if (!file.exists(file)){
		print(paste0(x, ' : No matched file!'))
		next
	}

	dnf <- read.table(file, sep='\t', header=F)
	colnames(dnf) <- c('cell', 'nf')
	dnf$cell <- paste0(d$sample, '_', dnf$cell)
	dnf$sample <- x

	d_inf <- rbind(d_inf, dnf)
}
# cell  nf  sample

write.table(d_inf, out, sep='\t', quote=F)

