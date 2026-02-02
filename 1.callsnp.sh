
cellsnp=~/miniconda3/envs/cellsnp/bin/cellsnp-lite

path=~/common/Analysis/Project_PrimaryRelapse/snatac/out_cellranger_atac_v2.2
outdir=~/common/Analysis/Project_PrimaryRelapse/snatac/out_cnv

sample='95_0020_P'

${cellsnp} \
    -s $path/$sample/outs/possorted_bam.bam \
    -b $path/outs/filtered_peak_bc_matrix/barcodes.tsv \
    -O $outdir/$sample/out_pileup \
    -R $outdir/$sample/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz \
    -p 25 \
    --minMAF 0 \
    --minCOUNT 2

