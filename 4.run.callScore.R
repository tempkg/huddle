
library(Seurat)
library(scatlasR)
library(ggplot2)
library(GetoptLong)



rds <- 'multiome.rds'

f <- 'hs.EPN-ZR.xlsx'
form <- 'sctype'

out <- 'out_expMarkers.confirm.pdf'

GetoptLong(
    "f=s",      "file",
    "rds=s",    "rds file",
    "form=s",   "form style",
    "outdir=s", "outdir"
)


dir.create(outdir)

seu <- readRDS(rds)
seu <- SeuratAddScorePlus(obj=seu, fdata=f, form=form)


df <- readxl::read_excel(f)
features <- paste0(unique(df$cellName), '1')

# show features in UMAP
p <- FeaturePlot(seu, features=features, reduction='wnn.umap', order=T, min.cutoff='q1')
ggsave(out, w=10, h=8)



