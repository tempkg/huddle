# Hua Sun
# 2025-01-13 v3
# seurat 5.10

library(Seurat)
library(Signac)
library(data.table)
library(SeuMultiome)
library(JASPAR2022)
library(GetoptLong)

set.seed(42)



name <- 'sample'
h5 <- 'filtered_feature_bc_matrix.h5'
ref <- 'mm10'
group <- ''
min_cells <- 3
atac_fragment <- 'atac_fragments.tsv.gz'
outdir <- 'out_seurat_signac.single.v3'

macs2 <- '/research/groups/mackgrp/home/common/Software/miniconda3/envs/macs2/bin/macs2'


GetoptLong(
    "name=s",      "sample name",
    "h5=s",        "h5 file",
    "ref=s",       "ref",
    "group=s",     "group",
    "min_cells=i", "min_cells",
    "atac_fragment=s", "atac_fragment",
    "outdir=s",    "outdir"
)


dir.create(outdir)

annotation <- Annotations(ref)
genome <- UCSCBSGenome(ref)


print(name)
#print(h5)

seu <- CreateMultiomeObjectFromH5(h5=h5, min_cells=min_cells, atac_fragment=atac_fragment, ref=ref, name=name, annotation=annotation)
print('[INFO] Created multiome object')
print(dim(seu))

seu$Sample <- name
seu$group <- group
seu$ref <- ref


# QC plot
QCPlotMultiome(seu, name, outdir)


# filter
seu <- FilterCells_multiome(seu)
print('[INFO] Filtered cells ...')
print(dim(seu))

cutoff <- CalMedianMaxCutoffN(multiome_list, 0.986, 0.85)
print(cutoff)
seu <- FilterCellsMultiomeX1(obj=seu, ncount_atac_max=cutoff[1], ncount_rna_max=cutoff[2], nfeature_rna_max=cutoff[3])

seu <- FilterDoubletByDoubletFinder(seu, outdir)
saveRDS(seu, file=paste0(outdir, '/multiome.filtered.rds'))


# normalize
DefaultAssay(seu) <- "RNA"
seu <- SCTransform(seu, vst.flavor = "v2", verbose = FALSE)
seu <- RunPCA(seu, assay = "SCT", npcs=50) %>% 
seu <- RunUMAP(seu, dims = 1:40, reduction.name='rna.umap', reduction.key='rnaUMAP_')
saveRDS(seu, file=paste0(outdir, '/multiome_exp.rds'))


# call peaks using macs2
seu <- CallPeaksUsingMACS2(seu, macs2, annotation)

# call pfm
pfm <- CallPFMFromJASPAR(JASPAR2022)
saveRDS(pfm, file = paste0(outdir, "/pfm.", ref, ".rds"))

# gene activity & motif
print('snatac_geneact_motif')
seu <- AddGeneActivityAndMotif(obj=seu, genome=genome, pfm=pfm)
saveRDS(seu, file=paste0(outdir, '/multiome_exp_act_motif.rds'))

# multiome object
DefaultAssay(seu) <- "peaks"
seu <- RunUMAP(seu, dims = 2:30, reduction='lsi', reduction.name='atac.umap', reduction.key='atacUMAP_')
seu <- CombinationATAC_RNA(seu)
saveRDS(seu, file=paste0(outdir, '/multiome_exp_act_motif.wnn.rds'))
write.table(seu@meta.data, file = paste0(outdir, "/multiome_exp_act_motif.wnn.metaData.xls"), sep="\t", quote=F, col.names = NA)


# out umap
pdf(paste0(outdir, '/umap.pdf'), w=6, h=6, useDingbats=T)
p1 <- DimPlot(seu, reduction='wnn.umap', pt.size=0.1)
p2 <- DimPlot(seu, reduction='rna.umap', pt.size=0.1)
p3 <- DimPlot(seu, reduction='atac.umap', pt.size=0.1)
print(p1)
print(p2)
print(p3)
dev.off()



