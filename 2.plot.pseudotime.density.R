# 9/28/23 v0.2

# R 4.2.1
# monocle v2.26.0

library(ggpubr)
library(ggplot2)
library(ggridges)
library(stringr)
library(dplyr)

library(GetoptLong)



title <- ''
sorted_label <- ''  # 'RGC-Like,OPC-Like,Fibroblast-Like,Neuronal-Like'
fdata <- 'out_monocle2/monocle2_ordercell.root.pseudotime.xls'

height <- 1.5
outdir <- 'out_monocle2plus'

GetoptLong(
    "title=s",         "title",
    "sorted_label=s",  "sorted_label",
    "fdata=s",         "pseudotime data file",
    "height=f",        "height",
    "outdir=s",        "output path"
)




#########################
##       Main
#########################

dir.create(outdir)

d_cds <- read.table(fdata, sep='\t', header=T)
colnames(d_cds) <- c('cells', 'pseudotime', 'cell_type2')


# sorted_label <- 'RGC-Like,OPC-Like,Fibroblast-Like,Neuronal-Like'
if (nchar(sorted_label)>1){ 
    sorted_label <- str_split(sorted_label,',')[[1]]
    d_cds <- filter(d_cds, cell_type2 %in% sorted_label)
    d_cds$cell_type2 <- factor(d_cds$cell_type2, level=rev(sorted_label))
}

p <- ggplot(d_cds, aes(x = pseudotime, y = cell_type2, fill = stat(x))) +
      geom_density_ridges_gradient(lwd = 0.1, scale = 1, gradient_lwd = 1.) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_discrete(expand = expand_scale(mult = c(0.01, 0.25))) +
      scale_fill_viridis_c(name = "Pseudotime", option = "C") +
      labs(title = title) +
      theme_ridges(font_size = 6, grid = FALSE) + 
      theme(axis.title = element_blank()) +
      theme(axis.text=element_text(size=6)) +
      theme(legend.position='bottom',
              legend.text=element_text(size=5),
              legend.key.height=unit(0.4,"line"),
              legend.key.size = unit(0.8, 'lines')) +
      theme(axis.ticks = element_line(linewidth = 0.3),
            axis.ticks.length=unit(1, "mm")) +
      geom_vline(xintercept = 5, linetype="dashed", color = "#696969", size=0.3) +
      geom_vline(xintercept = 15, linetype="dashed", color = "#696969", size=0.3)


pdf(paste0(outdir, '/plot_cell_dencity.pseudotime.pdf'), width = 2, height = height, useDingbats=FALSE)
print(p)
dev.off()








