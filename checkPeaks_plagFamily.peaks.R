
# perl -pe 's/$/\tST2\tMA1615.1/' MA1615.1.targetPeaks.cutoff10percentage.peaks > MA1615.1.targetPeaks.cutoff10percentage.peaks.formted.out
# cat */*.out > merge_with_zr-fus_bindingPeaks/plagFamily.targetPeaks.merged.out

# https://cran.r-project.org/web/packages/bedr/vignettes/Using-bedr.html

library(GenomicRanges)
library(rtracklayer)
library(GetoptLong)


f_peak <- 'suppTable2.peaks.txt'

# PLAG1
f_plag1 <- '/Users/hsun41/Documents/MackLab/projects/1.proj_plagl1/analysis/human.multiome/hs.integ.epn_19s.analysis.v2/infercnv_sample.snatac.hmm/multiome.malignant.All/analysis.All/motif_targetPeaks.PLAG_family/motif_targets_peaks.ST2/MA0163.1.targetPeaks.cutoff10percentage.peaks'
# PLAGL1
f_plagl1 <- '/Users/hsun41/Documents/MackLab/projects/1.proj_plagl1/analysis/human.multiome/hs.integ.epn_19s.analysis.v2/infercnv_sample.snatac.hmm/multiome.malignant.All/analysis.All/motif_targetPeaks.PLAG_family/motif_targets_peaks.ST2/MA1615.1.targetPeaks.cutoff10percentage.peaks'
# PLAGL2
f_plagl2 <- '/Users/hsun41/Documents/MackLab/projects/1.proj_plagl1/analysis/human.multiome/hs.integ.epn_19s.analysis.v2/infercnv_sample.snatac.hmm/multiome.malignant.All/analysis.All/motif_targetPeaks.PLAG_family/motif_targets_peaks.ST2/MA1548.1.targetPeaks.cutoff10percentage.peaks'

cutoff <- 0.5
outdir <- '.'

GetoptLong(
    "f_peak=s",   "peak",
    "cutoff=f",   "overlap percentage",
    "outdir=s",    "output path"
)

dir.create(outdir)

# tumor diff peaks
d_tar_data <- read.table(f_peak, sep='\t', header=T)
# Cancer    peak

# format to GRanges
tar_peaks <- as.data.frame(str_split_fixed(unique(d_tar_data$peak), "-", 3))
tar_peaks <- unique(tar_peaks)
# bed to GRanges
gr.tar_peaks <- GRanges(seqnames=tar_peaks$V1, 
                ranges=IRanges(start=as.numeric(tar_peaks$V2), end=as.numeric(tar_peaks$V3)))
length(gr.tar_peaks)




RunSingle <- function(f_ref, gr.tar_peaks)
{
    # plag peaks
    ref_peaks <- readLines(f_ref)
    #colnames(d_ref) <- c('peak', 'group', 'motif')
    # format to GRanges
    ref_peaks <- as.data.frame(str_split_fixed(unique(ref_peaks), "-", 3))
    ref_peaks <- unique(ref_peaks)
    # bed to GRanges
    gr.ref_peaks <- GRanges(seqnames=ref_peaks$V1, 
                    ranges=IRanges(start=as.numeric(ref_peaks$V2), end=as.numeric(ref_peaks$V3)))
    length(gr.ref_peaks)


    # correct way
    hits <- findOverlaps(gr.ref_peaks, gr.tar_peaks)

    # overlap ranges
    overlaps <- pintersect(gr.ref_peaks[queryHits(hits)], gr.tar_peaks[subjectHits(hits)])
    percentOverlap <- width(overlaps) / width(gr.tar_peaks[subjectHits(hits)])
    hits <- hits[ percentOverlap >= cutoff ]
    #     queryHits subjectHits
    #      <integer>   <integer>
    #  [1]         3        6849


    gr.ref_peaks<- gr.ref_peaks[queryHits(hits)]
    gr.tar_peaks <- gr.tar_peaks[subjectHits(hits)]



    # only queryHits start was automatically +1, not sure why
    bed_ref_peaks <- as.data.frame(gr.ref_peaks)[,1:3]
    bed_ref_peaks$start <- bed_ref_peaks$start - 1
    bed_ref_peaks$peak <- paste0(bed_ref_peaks$seqnames, '-', bed_ref_peaks$start, '-', bed_ref_peaks$end)

    # target
    bed_tar_peaks <- as.data.frame(gr.tar_peaks)[,1:3]
    bed_tar_peaks$peak <- paste0(bed_tar_peaks$seqnames, '-', bed_tar_peaks$start, '-', bed_tar_peaks$end)


    # add ref peak to target data
    bed_tar_peaks$ref_peak <- bed_ref_peaks$peak
    head(bed_tar_peaks)

    # merge to original target peak data
    d_tar_data_plus <- merge(d_tar_data, bed_tar_peaks[,c('peak', 'ref_peak')], by='peak', all.x = TRUE)




    #report <- as.data.frame(table(d_tar_data_plus$Cancer))
    #colnames(report) <- c('Cancer', 'n_peaks')
    with_ref <- as.data.frame(table(d_tar_data_plus[complete.cases(d_tar_data_plus), 'Cancer']))
    colnames(with_ref) <- c('Cancer', 'ref_peaks')

    #report <- merge(report, with_ref, by='Cancer', all.x=T)

    return(with_ref)
}




with_ref_plag1 <- RunSingle(f_plag1, gr.tar_peaks)
colnames(with_ref_plag1) <- c('Cancer', 'n_peaks')
with_ref_plag1$motif <- 'PLAG1'

with_ref_plagl1 <- RunSingle(f_plagl1, gr.tar_peaks)
colnames(with_ref_plagl1) <- c('Cancer', 'n_peaks')
with_ref_plagl1$motif <- 'PLAGL1'

with_ref_plagl2 <- RunSingle(f_plagl2, gr.tar_peaks)
colnames(with_ref_plagl2) <- c('Cancer', 'n_peaks')
with_ref_plagl2$motif <- 'PLAGL2'

report <- rbind(with_ref_plag1, with_ref_plagl1, with_ref_plagl2)



# PLAGL1 '#F37121', PLAGL2 '#387EC1', PLAG1 '#29B780'
col_set <- c('#29B780', '#387EC1', '#F37121')

col_group <- c(
        "PLAG1" = "#29B780",  
        "PLAGL1" = "#F37121",                
        "PLAGL2" = "#387EC1"
    )

# sorted
p <- ggplot(report, aes(x=forcats::fct_rev(fct_reorder(Cancer, n_peaks)), y=n_peaks, fill=motif)) +
        geom_bar(position="dodge", stat = "identity", width = 0.65) +
            #geom_bar(data=report, aes(x=Cancer, y=n_peaks), position="dodge", stat = "identity", width = 0.65) +
            theme_classic(base_line=0.1) +
            scale_fill_manual(values=col_group) +
            labs(y = "Number of peaks", x = "Cancer type") +
            theme(plot.title = element_text(hjust = 0.5, size=8),
                    text=element_text(hjust = 0.5, size=7)) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
            theme(legend.title = element_blank(),
                    legend.position=c(0.8, 0.9),
                    legend.text = element_text(size = 6),
                    legend.key.size=unit(2,"mm")) 


pdf(paste0(outdir, '/cancer_types.withPLAGFamilyPeaks_fromST2group.pdf'), width = 2.5, height = 1.5)
print(p)
dev.off()



