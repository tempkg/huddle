
file='sample.list'

outdir=./atac_fragments
mkdir -p $outdir

touch $outdir/combined_fragments.tsv

fpath=/research/groups/mackgrp/home/common/Analysis/Project_PrimaryRelapse/snatac/out_cellranger_atac.v2.2

cat $file | while read sample
do
	zcat $fpath/$sample/outs/fragments.tsv.gz | grep -v '#' | awk -v prefix="$sample" '{print $1 "\t" $2 "\t" $3 "\t" prefix "_" $4 "\t" $5}' >> $outdir/combined_fragments.tsv
done

sort -k1,1 -k2,2n $outdir/combined_fragments.tsv > $outdir/combined_fragments_sorted.tsv
bgzip $outdir/combined_fragments_sorted.tsv
tabix -p bed $outdir/combined_fragments_sorted.tsv.gz



