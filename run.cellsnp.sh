cellsnp=/research/groups/mackgrp/home/common/Software/miniconda3/envs/cellsnp/bin/cellsnp-lite

NAME=
BAM=
BARCODE=
DIR_PATH="/research/groups/mackgrp/home/common/Analysis/Project_PrimaryRelapse/snatac/out_cellsnps"

OUT_DIR=$DIR_PATH/$NAME.pileup
mkdir -p $OUT_DIR

#${cellsnp} -s $BAM -b $BARCODE -O $OUT_DIR -p 25 --minMAF 0.1 --minCOUNT 20 --cellTAG None --UMItag UB --gzip 
${cellsnp} -s $BAM -b $BARCODE -O $OUT_DIR -p 25 --minMAF 0.05 --minCOUNT 10 --cellTAG CB --UMItag None 

