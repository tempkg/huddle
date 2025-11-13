# Hua Sun
# Seurat 5.1
# 2025-11-08 v2.3


library(Seurat)
library(harmony)
library(data.table)
library(seuproc)
library(scqc)
library(dplyr)
library(GetoptLong)
library(ggplot2)


options(future.globals.maxSize=100000000000)



table <- 'data.table'
ref <- 'mouse'
min_cells <- 3
max_cells <- 1000000
seed <- 42
method <- 'harmony'  # 'rpca'
percent_mt <- 10
max_dim <- 30

outdir <- 'out_scrna_integrated'

GetoptLong(
    "table=s",      "table",
    "ref=s",       "ref",
    "seed=i",      "seed",
    "min_cells=i", "min_cells",
    "max_cells=i", "set max cell counts",
    "method=s",    "integration method",
    "percent_mt=i", "percent_mt",
    "max_dim=i",   "max_dim",
    "outdir=s",    "outdir"
)

set.seed(seed)
dir.create(outdir)

outdir_process <- paste0(outdir, '/processed')
dir.create(outdir_process)

seurat_list <- NULL

dt <- fread(table, header=FALSE, data.table=F)

for (name in dt$V1){
        print(name)

        rdir <- dt$V2[dt$V1==name]
        h5 <- paste0(rdir, '/filtered_feature_bc_matrix.h5')
        seu <- CreateSeuratObjectFromH5(h5=h5, min_cells=min_cells, name=name, ref=ref)
        print(dim(seu))

        # QC plot
        p <- QCPlotRNA(seu)
        ggsave(paste0(outdir_process, '/', name, '.rawdata.qc.pdf'), w=5, h=2.5, useDingbats=T)

        # filter
        cutoff <- CalMaxCutoffRNA(obj=seu, maxratioCount=0.99, maxratioFeature=0.9, revise_val=TRUE)
        print(cutoff)
        seu <- subset(x = seu,
            subset = nCount_RNA > 2000 &
                nFeature_RNA > 500 &
                nCount_RNA < cutoff[1] &
                nFeature_RNA < cutoff[2] &
                percent.mt < percent_mt
        )
        
        seu <- RunDoubletFinder(obj=seu, singlet=TRUE, outdir=outdir_process)
        #saveRDS(seu, file=paste0(outdir_process, '/', name, '.filtered.remDoublet.rds'))
        
        seu <- RunFilterEmptyDrops(dir=rdir, obj=seu, prefix=name, outdir=outdir_process)
        #saveRDS(seu, file=paste0(outdir_process, '/filtered.remDroplet.rds'))

        # make data list
        seurat_list <- append(seurat_list, seu)
        
}




# extract top x feature cells
if (max_cells < 1000000){
    print('[INFO] extract cells ...')
    seurat_list <- base::lapply(X = seurat_list, FUN = function(x) { x <- ExtractNCells(x, ncells=max_cells) })
}

saveRDS(seurat_list, file=paste0(outdir,'/filtered.remDoublet_emptyDrops.list.rds'))



# Seurat5
seu <- MergeObjects(seurat_list)
DefaultAssay(seu) <- 'RNA'

seurat_list <- NULL

seu <- Seurat5SCTNormalize(obj=seu, regress="percent.mt", variable=3000)
saveRDS(seu, file=paste0(outdir,'/scrna.SCT2Normalize.rds'))

# Integration
integ_rna <- RunHarmonyIntegration(obj=seu, norm_method='SCT', reduc='pca', new_reduc='integrated.rna', umap_reduc='umap', max_dim=30, res_clus=0.4, seed=seed)
saveRDS(integ_rna, file=paste0(outdir,'/scrna_sct2_integrated.harmony.rds'))
write.table(integ_rna@meta.data, file = paste0(outdir, '/scrna_sct2_integ.metadata.xls'), sep = "\t", quote=F, col.names = NA)

seu <- NULL


# plot umap
p <- DimPlot(integ_rna, group.by = 'Sample', reduction = "umap", pt.size = 0.1)
ggsave(paste0(outdir, '/scrna_sct2_integ.umap.pdf'), width = 5.5, height = 4)


pdf(paste0(outdir, '/scrna_sct2_integ.elbowplot.pdf'), w=4, h=2)
ElbowPlot(integ_rna, ndims = max_dim, reduction = "integrated.rna")
dev.off()





