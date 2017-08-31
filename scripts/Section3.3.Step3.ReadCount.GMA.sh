#!/bin/sh
### 2017-07-18 by Jiyoung Lee 
### To count reads on features (gene models) from BAM.
### FeatureCounts handles single/paired-end data seperately.

PROJ_PATH=$(pwd)

### Setting input files ==================================================================
GTF=$PROJ_PATH/raw_data/Gmax_275_Wm82.a2.v1.gene_exons.geneid.gff3

### Readcounting using FeatureCounts =====================================================
cd $PROJ_PATH/raw_data/fastq/GMA;
for f1 in *.fastq.gz;
do
	f2=${f1%.fastq.gz};
	echo -e "`date +"%b%d|%T"`\t$f2";

	mkdir $PROJ_PATH/processed_data/rc/$f2

	bamsuffix=Aligned.sortedByCoord.out.bam
	BAM=$f2$bamsuffix

	### Processing Paired-end fastq fiels ================================================
  featureCounts -t exon \
			-g gene_id \
			-p \
			-a $GTF \
			-o $PROJ_PATH/processed_data/rc/$f2/$f2.readcount.txt \
			$PROJ_PATH/processed_data/bam/$f2/$BAM >$PROJ_PATH/logs/FtCnt_$f2.log 2>$PROJ_PATH/logs/FtCnt_$f2.err

