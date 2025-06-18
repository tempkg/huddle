# Hua Sun

# Calculate module scores for feature expression programs in single cells
# AddModuleScore (Seurat function)

# v2 change cutoff
# human c("> 0.3", "0.2-0.3", "< 0.2")
# mouse c("> 0.2", "0.1-0.2", "< 0.1")

# Set Library
library(dplyr)
library(ggplot2)
library(RColorBrewer)


######################################
##            Set Func.
######################################


UMAPPlot_Score <- function(df=NULL, 
    breaks=NULL, 
    colorset=c("#CB4335", "#2E86C1", "#D7DBDD"), 
    label_name=NULL,
    point_size=0.1
){
    df <- df[order(df$signal, decreasing=FALSE), ]

    p <- ggplot(df,aes(x=UMAP_1, y=UMAP_2, color=group)) + 
            geom_point(size=point_size) +
            theme_void()
            #theme_classic(base_line_size=0.2) + 
            #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            #    axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
            # labs(x = "UMAP 1", y="UMAP 2")

    p <- p + scale_color_manual(name=label_name, breaks=breaks, values=colorset)

    p <- p + guides(color = guide_legend(override.aes = list(size = 2))) +
            theme(legend.title = element_text(size=8, face="bold")) +
            theme(text=element_text(size=8, face="bold"))

    p
}



## UMAP per sample
SignalUMAP_perSample <- function(df, 
    breaks=NULL, 
    colorset=c("#CB4335", "#2E86C1", "#D7DBDD"),
    label_name=NULL,
    point_size=0.1, 
    sorted_sample='', 
    cutoff=''
){
    df <- df[order(df$signal, decreasing=FALSE), ]
    

    if (length(sorted_sample) > 2){ 
        df$Sample <- factor(df$Sample, levels=sorted_sample) 
    }

    p <- ggplot(df,aes(x=UMAP_1, y=UMAP_2, color=group)) + 
            geom_point(size=point_size) +
            theme_classic(base_line_size=0.2) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
            labs(x = "UMAP 1", y="UMAP 2")

    p <- p + scale_color_manual(name=label_name, breaks=breaks, values=colorset)

    p <- p + guides(color = guide_legend(override.aes = list(size = 2))) +
            theme(legend.title = element_text(size=9))

    # per sample
    p <- p + facet_wrap(~Sample, ncol = 5) + 
        theme(strip.background = element_blank(), strip.text = element_text(size=10))
    
    p
}

