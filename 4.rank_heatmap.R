

library(ComplexHeatmap)
library(GetoptLong)


input <- ''
out <- 'rank.heatmap.pdf'

GetoptLong(
    "input=s",    "input file",
    "out=s",   "output path"
)




# input <- '/Users/hsun41/Documents/MackLab/projects/1.proj_plagl1/analysis/human.multiome/hs.integ.epn_19s.analysis.v2/infercnv_sample.snatac.hmm/multiome.malignant.All/analysis.All/markers_peaks/out_peakMarkers.st2_vs_rest/ST2.rankedMotifs.xls'
d <- read.table(input, sep='\t', header=T)
d$motif_rank <- -log10(d$pvalue)  # Motif probability



# make label
label_anno <- d[,'motif.name']
names(label_anno) <- rownames(d)
label_anno <- label_anno[1:10]

label_gene_at <- as.numeric(names(label_anno))
label_gene <- as.vector(label_anno)

#
# row_ha = rowAnnotation(link = anno_mark(at=label_gene_at, labels=label_gene, labels_gp=gpar(fontsize = 6), extend=unit(0.5,'mm'), link_width=unit(0.1,'mm')))
row_ha = rowAnnotation( link = anno_mark(at=label_gene_at, labels=label_gene, labels_gp=gpar(fontsize = 6)))

ht <- Heatmap(d$motif_rank, 
            name = "-log10(pvalue)", 
            col = colorRampPalette(c('#2E86C1', '#F9E79F', '#E74C3C'))(100),

            show_row_names = F,
            show_row_dend = F,
            cluster_rows = T,

            show_column_names = F,

            row_title = 'Motif Rank', 
            row_title_gp = gpar(fontsize = 6),
            right_annotation = row_ha,

            heatmap_legend_param = list(
                #title_position = "lefttop-rot",
                #direction = "vertical",
                title_position = "lefttop",
                direction = "horizontal", 
                title_gp = gpar(fontsize = 7), 
                labels_gp = gpar(fontsize = 6),
                grid_width = unit(1.5, "mm"),
                grid_height = unit(1.5, "mm")
            )
        )


#lgd = Legend(col_fun = col_fun, title = "foo", at = c(0, 0.25, 0.5, 0.75, 1))

pdf(out, width = 1.1, height = 1.5, useDingbats=FALSE)
draw(ht, heatmap_legend_side = "bottom", show_heatmap_legend=TRUE)
dev.off()




