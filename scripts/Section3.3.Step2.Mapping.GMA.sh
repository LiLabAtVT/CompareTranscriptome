#!/bin/sh
### RNA-seq pipeline: Mapping fastq files

### Setting input files ==================================================================
GNM=$PROJ_PATH/raw_data/Gmax_275_v2.0.fa
GTF=$PROJ_PATH/raw_data/Gmax_275_Wm82.a2.v1.gene_exons.gff3
IDX=$PROJ_PATH/raw_data/GMX_STAR-2.5.2b_index

### RNA-seq pipeline =====================================================================
mkdir -p $PROJ_PATH/processed_data/bam;
cd $PROJ_PATH/raw_data/fastq/GMA
for f1 in *.fastq.gz;
do
	f2=${f1%.fastq.gz};
	
	### Processing Paired-end fastq fiels ============================================
		f3=$f2"_1.fastq.gz"
		f4=$f2"_2.fastq.gz"

		cd $PROJ_PATH/processed_data/bam
		mkdir $f2
		STAR --genomeDir $IDX \
			--readFilesIn $PROJ_PATH/raw_data/fastq/GMA/$f3 $PROJ_PATH/raw_data/fastq/GMA/$f4 \
			--outFileNamePrefix $PROJ_PATH/processed_data/bam/$f2/$f2 \
			--outSAMtype BAM SortedByCoordinate >$PROJ_PATH/logs/STAR_$f2.log;

done
