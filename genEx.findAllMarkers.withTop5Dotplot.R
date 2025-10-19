
library(Seurat)
library(seuextra)
library(stringr)
library(dplyr)
library(GetoptLong)
library(RColorBrewer)
library(ggplot2)
library(data.table)


rds <- 'ZR.multiome_integrate.30_30.rds'
idents <- 'seurat_clusters'
outdir <- 'out_expMarkers.zr'
dir.create(outdir)


seu <- readRDS(rds)

DefaultAssay(seu) <- 'SCT'

seu <- PrepSCTFindMarkers(seu) 
diff_exp <- RunFindAllMarkers(seu, assay_use='SCT', idents=idents, recorrect_umi=FALSE)

write.table(diff_exp, file = paste0(outdir, '/findallmarker.geneExp.xls'), sep = "\t", quote=FALSE, row.names = F)



diff_exp_filtered <- diff_exp %>% 
                    filter(!grepl('^MT-|^MTRNR', gene)) %>%  
                    filter(!grepl('^RPL|^RPS|^MRP', gene)) %>% 
                    filter(!grepl('^AC[0-9]|^AL[0-9]', gene)) %>% 
                    filter(!grepl('^LINC', gene)) %>%
                    filter(p_val_adj < 0.0001 & avg_log2FC > 1 & pct.1 > 0.3)

write.table(diff_exp_filtered, file = paste0(outdir, '/', idents, '.findallmarker.geneExp.filtered.xls'), sep = "\t", quote=FALSE, row.names = F)


# plot
p <- Seurat_DotPlotMarkers(obj=seu, group=idents, d_markers=diff_exp_filtered)
ggsave(paste0(outdir, '/top5_markers.pval.pdf'), w=5, h=3)


dx_log2fc <- diff_exp_filtered %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = T)
p <- Seurat_DotPlotMarkers(obj=seu, group=idents, d_markers=dx_log2fc)
ggsave(paste0(outdir, '/top5_markers.logfc.pdf'), w=5, h=3)








