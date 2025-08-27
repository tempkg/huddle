path='out_cellranger_arc'
outdir='/path'

zcat $path/ZR0/outs/atac_fragments.tsv.gz | grep -v '#' | awk '{print $1 "\t" $2 "\t" $3 "\t" "ZR0_" $4 "\t" $5}' > $outdir/combined_fragments.tsv
zcat $path/ZR1/outs/atac_fragments.tsv.gz | grep -v '#' | awk '{print $1 "\t" $2 "\t" $3 "\t" "ZR1_" $4 "\t" $5}' >> $outdir/combined_fragments.tsv
zcat $path/ZR2/outs/atac_fragments.tsv.gz | grep -v '#' | awk '{print $1 "\t" $2 "\t" $3 "\t" "ZR2_" $4 "\t" $5}' >> $outdir/combined_fragments.tsv

sort -k1,1 -k2,2n $outdir/combined_fragments.tsv > $outdir/combined_fragments_sorted.tsv
bgzip $outdir/combined_fragments_sorted.tsv
tabix -p bed $outdir/combined_fragments_sorted.tsv.gz



