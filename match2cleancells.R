# Hua Sun

library(data.table)
library(dplyr)
library(GetoptLong)


fmd <- 'filtered_feature_bc_matrix.h5'
flb <- 'cell_lineageBarcodes.filtered.txt'
outdir <- 'out_cleancell_lb'


GetoptLong(
    "meta=s",   "metadata",
    "flb=s",    "barcode file",
    "outdir=s", "outdir"
)


metadata <- read.table(fmd, sep='\t', header=T, row.names=1)
lineageBarcode <- read.table(flb, sep='\t', header=F)
lineageBarcode$V1 <- paste0(lineageBarcode$V1, '-1')
lineageBarcode <- aggregate(V2 ~ V1, data = lineageBarcode, FUN = function(x) paste(x, collapse = ","))

metadata$lineageBarcode <- lineageBarcode$V2[match(rownames(metadata), lineageBarcode$V1)]

write.table(metadata, file=paste0(outdir, '/metadata.added_lineageBarcode.xls'), sep='\t', quote=F, col.names=NA)


clean_data <- data.frame('cell'=rownames(metadata), 'lineageBarcode'=metadata$lineageBarcode)
clean_data <- na.omit(clean_data)
write.table(clean_data, file=paste0(outdir, '/cleancells_lineagebarcode.tsv'), sep='\t', quote=F, row.names=F)



