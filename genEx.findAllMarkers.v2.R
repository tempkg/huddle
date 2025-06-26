
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

GetoptLong(
    "rds=s",      ".rds seurat file",
    "ref=s",    "ref",
    "recorrect_umi",    "recorrect_umi",
    "idents=s",   "group",
    "fmeta=s",    "fmeta"
)

outdir <- paste0('out_expMarkers.', idents)
dir.create(outdir)

seu <- readRDS(rds)

# fmeta
if (fmeta != ''){
    metadata <- fread(fmeta, data.table=F)
    rownames(metadata) <- metadata[,1]
    metadata[,1] <- NULL
    seu@meta.data <- metadata
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
                    filter(avg_log2FC > 1) %>% 
                    filter(pct.1 > 0.5 & pct.2 < 0.1)
write.table(diff_exp_filtered, file = paste0(outdir, '/', idents, '.findallmarker.geneExp.filtered.xls'), sep = "\t", quote=FALSE, row.names = F)









