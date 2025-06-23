
library(GetoptLong)
library(infercnv)



f_obj <- ''
outdir <- "out_infercnv"


GetoptLong(
    "f_obj=s",      "infercnv obj",
    "outdir=s",     "output path"
)


dir.create(outdir)

infercnv_obj = readRDS(f_obj)

# https://www.rdocumentation.org/packages/infercnv/versions/1.3.3/topics/plot_cnv
plot_cnv(infercnv_obj,
         out_dir=outdir,
         cluster_by_groups=TRUE,
         hclust_method='ward.D2',
         color_safe_pal=FALSE,
         x.center=1,
         output_range = c(0.8, 1.2),  # Focus on 0.8 to 1.2
         threshold = c(0.8, 1.2),     # Truncate values outside this range
         output_filename="heatmap_custom_range",
         output_format="pdf"
)



