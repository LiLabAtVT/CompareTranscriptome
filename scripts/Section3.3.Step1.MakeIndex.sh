#!/bin/sh
### RNA-seq pipeline: make index for read mapping.

### The default directory for running this script is ATH_GMA, which is created in seciton 2.
WORKDIR=$(pwd)

IDX=$WORKDIR/raw_data/ATH_STAR-2.5.2b_index
GNM=$WORKDIR/raw_data/TAIR10_Chr.all.fasta
GTF=$WORKDIR/raw_data/Araport11_GFF3_genes_transposons.201606.gtf
STAR --runMode genomeGenerate \
     --genomeDir $IDX \
     --genomeFastaFiles $GNM \
     --sjdbGTFfile $GTF

IDX=$WORKDIR/raw_data/GMX_STAR-2.5.2b_index
GNM=$WORKDIR/raw_data/Gmax_275_v2.0.fa
GTF=$WORKDIR/raw_data/Gmax_275_Wm82.a2.v1.gene_exons.gff3
STAR --runMode genomeGenerate \
     --genomeDir $IDX \
     --genomeFastaFiles $GNM \
     --sjdbGTFfile $GTF \
     --sjdbGTFtagExonParentTranscript Parent

