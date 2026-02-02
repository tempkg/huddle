
cellsnp=~/software/miniconda3/envs/cellsnp/bin/cellsnp-lite
eagle=~/software/miniconda3/envs/cellsnp/bin/eagle

nc=8 # number of cores to use
sample="95_0020_P"
phase_panel="~/software/1000G/hg38/1000G_hg38"
vcf_genome1k="~/software/1000G/hg38/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
gma_gz="~/software/1000G/hg38/genetic_map_hg38_withX.txt.gz"
outdir="out_phase"

Rscript pileup_and_phaseX.R \
  --cellsnp ${cellsnp} \
  --eagle ${eagle} \
  --label ${sample} \
  --samples ${sample} \
  --bams possorted_bam.bam \
  --barcodes barcodes.tsv \
  --gmap ${gma_gz} \
  --snpvcf ${vcf_genome1k} \
  --paneldir ${phase_panel} \
  --ncores ${nc} \
  --cellTAG CB \
  --UMItag Auto,None \
  --outdir ${outdir}



