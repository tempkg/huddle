# Hua Sun
# 2026-1-21 v2.4
# seurat v5.4.0


library(Seurat)
library(harmony)
library(data.table)
library(seuproc)
library(scqc)
library(dplyr)
library(GetoptLong)
library(ggplot2)


options(future.globals.maxSize=100000000000)



ids <- 'sample_list.txt'
sdir <- 'out_seurat5'
frds <- 'scrna.remDoublet_emptyDrops.sct2.rds'

seed <- 42
max_cells <- 10000
method <- 'harmony'  # 'rpca'
max_dim <- 30

outdir <- 'out_scrna_integrated'

GetoptLong(
    "ids=s",      "id list",
    "sdir=s",     "seurat result dir.",
    "frds=s",     "read file",
    "seed=i",      "seed",
    "max_cells=i", "set max cell counts",
    "method=s",    "integration method",
    "max_dim=i",   "max_dim",
    "outdir=s",    "outdir"
)

set.seed(seed)

dir.create(outdir)


# read
sample_list <- readLines(ids)

seurat_list <- NULL
for (name in sample_list){
        print(name)

        fseu <- paste0(sdir, '/', name, '/', frds)
        temp_seu <- readRDS(fseu)

        # make data list
        seurat_list <- append(seurat_list, temp_seu)
        
        temp_seu <- NULL
}



# extract top x feature cells
print('[INFO] extract cells ...')
seurat_list <- base::lapply(X = seurat_list, FUN = function(x) { x <- ExtractNCells(x, ncells=max_cells) })

saveRDS(seurat_list, file=paste0(outdir,'/filtered.remDoublet_emptyDrops.list.rds'))



# Seurat5
combined_rna <- MergeObjects(seurat_list)
DefaultAssay(combined_rna) <- 'RNA'

seurat_list <- NULL

combined_rna <- Seurat5SCTNormalize(obj=combined_rna, regress="percent.mt", variable=3000)
saveRDS(combined_rna, file=paste0(outdir,'/scrna.SCT2Normalize.rds'))

# Integration
integ_rna <- RunHarmonyIntegration(obj=combined_rna, norm_method='SCT', reduc='pca', new_reduc='integrated.rna', umap_reduc='rna.umap', max_dim=30, res_clus=0.4, seed=seed)

# Seurat v5
integ_rna[['RNA']] <- JoinLayers(integ_rna[['RNA']])
saveRDS(integ_rna, file=paste0(outdir,'/scrna_sct2_integrated.harmony.rds'))
write.table(integ_rna@meta.data, file = paste0(outdir, '/scrna_sct2_integ.metadata.xls'), sep = "\t", quote=F, col.names = NA)

combined_rna <- NULL


# plot umap
p <- DimPlot(integ_rna, group.by = 'Sample', reduction = "rna.umap", pt.size = 0.1)
ggsave(paste0(outdir, '/scrna_sct2_integ.umap.pdf'), width = 5.5, height = 4)


pdf(paste0(outdir, '/scrna_sct2_integ.elbowplot.pdf'), w=4, h=2)
ElbowPlot(integ_rna, ndims = max_dim, reduction = "integrated.rna")
dev.off()





