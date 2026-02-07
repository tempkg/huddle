# Hua Sun


library(Seurat)
library(Signac)
library(ggplot2)
library(data.table)
library(GetoptLong)


#fct <- 'celltype.info'

#GetoptLong(
#    "fct=s",     "dir outs",
#    "outdir=s",  "outdir"
#)

#dir.create(outdir)


# if use seuproc package
#gene.activity <- GeneActivity(snatac)
#snatac <- AddGeneActivity(obj=snatac, gene_activity=gene.activity, assay='peaks', new_assay='RNA')


fct <- 'cluster_cellType.xls'
snatac <- readRDS('snatac.harmony.with_geneact.rds')
snrna <- readRDS('snrna.harmony.sct_normalized.rds')

info <- fread(fct, data.table=F)
snrna$celltype <- info$cell_type2[match(snrna$seurat_clusters, info$seurat_clusters)]


# FindVariableFeatures in SCT or LogNorm are same
snrna <- FindVariableFeatures(
  object = snrna,
  nfeatures = 5000
)



snrna[["RNA"]] <- split(snrna[["RNA"]], f = snrna$Sample)

DefaultAssay(snrna) <- 'RNA'

# run standard anlaysis workflow
snrna <- NormalizeData(snrna)
snrna <- FindVariableFeatures(snrna)
snrna <- ScaleData(snrna)



DefaultAssay(snatac) <- 'RNA'
transfer.anchors <- FindTransferAnchors(
  reference = snrna,
  query = snatac,
  reduction = 'cca',
  dims = 1:30
)


predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = snrna$celltype,
  weight.reduction = snatac[['lsi']],
  dims = 2:30
)


snatac <- AddMetaData(object = snatac, metadata = predicted.labels)


plot1 <- DimPlot(snrna, group.by = 'celltype', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('snRNA-seq')
plot2 <- DimPlot(snatac, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('snATAC-seq')


pdf('snatac.anno.pdf', w=8, h=4)
plot1 + plot2
dev.off()


# in the snatac, 'predicted.id' will have multiple cell type name in the same 'seurat_clusters'





