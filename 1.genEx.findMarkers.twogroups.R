
library(Seurat)
library(seuextra)
library(stringr)
library(dplyr)
library(GetoptLong)
library(RColorBrewer)
library(ggplot2)
library(data.table)


rds <- '2.multiome_integrated_plus.regout_mt.malignantCellOnly.rds'
idents <- 'group'
finfo <- 'sample.group.text'

outdir <- 'out_expMarkers.zr_ewp'
dir.create(outdir)


info <- fread(finfo, data.table=F)
seu <- readRDS(rds)
seu$group <- info$group[match(seu$orig.ident, info$sample)]


DefaultAssay(seu) <- 'SCT'
Idents(seu) <- 'group'

seu <- PrepSCTFindMarkers(seu) 
diff_exp <- RunFindMarkers(obj=seu, ct_target='ZR', ct_bg='EWP', recor_umi=FALSE, outdir=outdir)

write.table(diff_exp, file = paste0(outdir, '/findallmarker.geneExp.twogroups.xls'), sep = "\t", quote=FALSE, row.names = F)










