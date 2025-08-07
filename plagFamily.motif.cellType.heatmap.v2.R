
# Set Library
library(GetoptLong)
library(Seurat)
library(scapekit)
library(dplyr)
library(stringr)




rds <- 'seurat.rds'
motif_id <- 'MA1615.1,MA1548.1,MA0163.1'
out='out_heatmap.avgMotifSig.plag_family.v2.pdf'


seurat_obj <- readRDS(rds)

#-- motifs
tar_motifs <- str_split(motif_id, ',')[[1]]

df <- data.frame(t(seurat_obj@assays$chromvar@data[tar_motifs,]))
#                               MA1615.1    MA1548.1    MA0163.1
# H_02_2138_AAACAGCCAAGTGAAC.1  0.9646257  0.03137187  0.23462118

# make zero for <0
df[df < 0] <- 0
# only use >0 signal cells
df <- df[rowSums(df[])>0,]



# info
# cell cell_type2  group
# add extra info
df$group <- info$group[match(rownames(df), info$cell)]
df$cell_type2 <- info$cell_type2[match(rownames(df), info$cell)]




# make mean motifs per cell type
df_avg <- df %>% 
            group_by(group, cell_type2) %>% 
            summarize(PLAGL1 = mean(MA1615.1), PLAGL2 = mean(MA1548.1), PLAG1 = mean(MA0163.1), .groups = 'drop')
df_avg <- data.frame(df_avg)

d_mtx <- NULL
for (x in unique(df_avg$group)){
    dtemp <- df_avg %>% filter(group == x)
    rownames(dtemp) <- dtemp$cell_type2
    dtemp[,1:2] <- NULL
    colnames(dtemp) <- paste(colnames(dtemp), x, sep=".")
    ifelse(nchar(d_mtx)==0, d_mtx <- dtemp, d_mtx <- cbind(d_mtx, dtemp))
}

# scale by row 0-1 Not Standardize
#scaled_df <- t(apply(d_mtx, 1, function(x) (x - min(x)) / (max(x) - min(x))))

# scale by row: divides by standard deviation
scaled_df <- as.data.frame(t(scale(t(d_mtx))))

ComplexHeatmap_Motif(data=scaled_df, out=out)






