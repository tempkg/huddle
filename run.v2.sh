
dir='out_cellranger.v9'

ls $dir | while read sample
do
	echo $sample
	dir2=$dir/$sample/outs
	sh ~/space/hpc.bsub.sh 16 2 $name Rscript soupx.R --dir $dir2
done

