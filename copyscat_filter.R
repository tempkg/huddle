
# non-neoplastic barcodes
fbc <- 'immune_cells.txt'
fnc <- 'predicted_normal_cells.out'
baseline1 <- readLines(fbc)
baseline2 <- readLines(fnc)
baseline <- unique(c(baseline1, baseline2))


# read score file
file_path <- Sys.glob(paste0(outdir, "/*_clean_b2_cnv_scores.csv"))
dscore <- read.csv(file_path[1], sep=',', header=T, row.names=1)
dscore_filtered <- dscore[!dscore$rowname %in% baseline,]


all_chromosomes <- c('chr1p', 'chr1q', 'chr2p', 'chr2q', 'chr3p', 'chr3q', 'chr4p', 'chr4q', 'chr5p', 'chr5q', 
    'chr6p', 'chr6q', 'chr7p', 'chr7q', 'chr8p', 'chr8q', 'chr9p', 'chr9q', 'chr10p', 'chr10q', 
    'chr11p', 'chr11q', 'chr12p', 'chr12q', 'chr13p', 'chr13q', 'chr14p', 'chr14q', 'chr15p', 'chr15q', 
    'chr16p', 'chr16q', 'chr17p', 'chr17q', 'chr18p', 'chr18q', 'chr19p', 'chr19q', 'chr20p', 'chr20q', 
    'chr21p', 'chr21q', 'chr22p', 'chr22q')


data2 <- data.frame(rowname = dscore_filtered$rowname)
data2[, all_chromosomes] <- 0

# Match and fill values
matching_cols <- intersect(names(dscore_filtered)[-1], all_chromosomes)
data2[, matching_cols] <- dscore_filtered[, matching_cols]

# View the result
#print(data2)

write.table(data2, file=paste0(outdir, '/clean_b2_cnv_scores.filtered_nonNeoplastic.tsv'), sep='\t', quote=FALSE, row.names=FALSE)






