
library(Seurat)
library(seuextra)
library(stringr)
library(dplyr)
library(GetoptLong)
library(RColorBrewer)
library(ggplot2)
library(data.table)


rds <- 'multiome_integrated.plus.rds'
idents <- 'seurat_clusters'
fmeta <- ''
ref <- 'human'
prefix <- ''

GetoptLong(
    "rds=s",      ".rds seurat file",
    "ref=s",    "ref",
    "recorrect_umi",    "recorrect_umi",
    "idents=s",   "group",
    "prefix=s",    "prefix",
    "fmeta=s",    "fmeta"
)



outdir <- paste0('out_expMarkers.', idents)
if (prefix != ''){
    outdir <- paste0(prefix, '.out_expMarkers.', idents)
}
dir.create(outdir)

print(outdir)

seu <- readRDS(rds)

# fmeta
if (fmeta != ''){
    meta <- read.table(fmeta, sep='\t', header=T, row.names=1)
    seu[[idents]] <- meta[[idents]][match(rownames(seu@meta.data), rownames(meta))]
}


DefaultAssay(seu) <- 'SCT'


# https://satijalab.org/seurat/articles/integration_introduction
# If running on a subset of the original object after running PrepSCTFindMarkers(), FindMarkers() should be invoked with recorrect_umi = FALSE to use the existing corrected counts:
seu <- PrepSCTFindMarkers(seu) 

if (recorrect_umi){
    diff_exp <- GenxFindAllMarkers(seu, assay_use='SCT', idents=idents)
} else {
    diff_exp <- GenxFindAllMarkers(seu, assay_use='SCT', idents=idents, recorrect_umi=FALSE)
}

write.table(diff_exp, file = paste0(outdir, '/', idents, '.findallmarker.geneExp.xls'), sep = "\t", quote=FALSE, row.names = F)


#//---------- Only use MGI for mouse

ExtractStandardGeneSymbols <- function(genes=NULL, ref='mouse', omit_na=FALSE){
    library(HGNChelper)

    valid_gene <- checkGeneSymbols(genes, species = ref, expand.ambiguous = FALSE)
    valid_gene <- valid_gene$Suggested.Symbol

    # Remove entries without MGI IDs
    if (omit_na) {
        valid_gene <- as.vector(na.omit(valid_gene))
    }

    return(valid_gene)
}


if (ref == 'mouse'){
    valid_gene <- ExtractStandardGeneSymbols(genes=diff_exp$gene, ref='mouse', omit_na=TRUE)
    diff_exp <- diff_exp[diff_exp$gene %in% valid_gene,]
}

#----------\\


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



