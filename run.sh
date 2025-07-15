

cat data.tsv | while read name form path
do
	echo $name
	sh ~/space/hpc.bsub.sh 16 2 $name ~/kit/rscript442 runqc.R --dir $path --form $form

done

