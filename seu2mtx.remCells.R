# 2025-03-27 v0.1

library(dplyr)
library(data.table)


fcell <- 'poor_ref_cells.out'
fgeneloci <- 'geneLoci.txt'
fcellinfo <- 'cellInfo.txt'
fmtx <- 'expCount.txt'

outdir <- 'out_infercnv_sec'
dir.create(outdir)

poor_cells <- readLines(fcell)
cellinfo <- fread(fcellinfo, data.table=F, header=F)
cellinfo_filtered <- cellinfo %>% filter(!V1 %in% poor_cells)

mtx <- fread(fmtx, data.table=F)
rownames(mtx) <- mtx[,1]
mtx[,1] <- NULL
mtx_filtered <- mtx_filtered[, !colnames(mtx_filtered) %in% poor_cells]


file.copy(from = fgeneloci, to=paste0(outdir, '/geneLoci.txt'))
fwrite(cellinfo_filtered, file=paste0(outdir, '/cellInfo.txt'), sep='\t', row.names=F, col.names=F)
fwrite(mtx_filtered, file=paste0(outdir, '/expCount.txt'), sep = "\t", quote=F, row.names=T, col.names=T)





