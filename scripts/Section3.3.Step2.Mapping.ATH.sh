#!/bin/sh
### RNA-seq pipeline: Mapping fastq files

PROJ_PATH=$(pwd)
### Setting input files ==================================================================
GNM=$PROJ_PATH/raw_data/TAIR10_Chr.all.fasta;
GTF=$PROJ_PATH/raw_data/Araport11_GFF3_genes_transposons.201606.gtf
IDX=$PROJ_PATH/raw_data/ATH_STAR-2.5.2b_index

### RNA-seq pipeline =====================================================================
mkdir -p $PROJ_PATH/processed_data/ATH/bam;
cd $PROJ_PATH/raw_data/ATH
for f1 in *.fastq.gz;
do
	f2=${f1%.fastq.gz};
	
	### Processing Paired-end fastq fiels ============================================
		f3=$f2"_1.fastq.gz"
		f4=$f2"_2.fastq.gz"

		cd $PROJ_PATH/processed_data/bam
		mkdir $f2
		STAR --genomeDir $IDX \
			--readFilesIn $PROJ_PATH/raw_data/fastq/ATH/$f3 $PROJ_PATH/raw_data/fastq/ATH/$f4 \
			--outFileNamePrefix $PROJ_PATH/processed_data/bam/$f2/$f2 \
			--outSAMtype BAM SortedByCoordinate >$PROJ_PATH/logs/STAR_$f2.log;

done
