
library(Seurat)
library(scatlasR)
library(ggplot2)
library(GetoptLong)



rds <- 'multiome.rds'

f <- 'hs.EPN-ZR.xlsx'
form <- 'sctype'

# 3x3 w=10,h=8
w <- 10
h <- 8

out <- 'confirm_markers.epn-zr.pdf'

GetoptLong(
    "f=s",      "file",
    "rds=s",    "rds file",
    "form=s",   "form style",
    "w=i",      "width",
    "h=i",      "height",
    "out=s",    "out pdf"
)


seu <- readRDS(rds)
seu <- SeuratAddScorePlus(obj=seu, fdata=f, form=form)


df <- NULL
features <- NULL

if (form == 'sctype'){
    df <- readxl::read_excel(f)
    features <- paste0(unique(df$cellName), '1')
}
if (form == 'seurat'){
    df <- read.table(f, sep='\t', header=T)
    features <- paste0(unique(df$cluster), '1')
}

# show features in UMAP
p <- FeaturePlot(seu, features=features, reduction='wnn.umap', order=T, min.cutoff='q1')
# 3x3 image
ggsave(out, w=w, h=h)



