#!/bin/sh
### 2017-07-18 by Jiyoung Lee 
### To count reads on features (gene models) from BAM.
### FeatureCounts handles single/paired-end data seperately.

PROJ_PATH=$(pwd)

### Setting input files ==================================================================
GTF=$PROJ_PATH/raw_data/Araport11_GFF3_genes_transposons.201606.gtf

### Readcounting using FeatureCounts =====================================================
cd $PROJ_PATH/raw_data/fastq/ATH;
for f1 in *.fastq.gz;
do
	f2=${f1%.fastq.gz};

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

