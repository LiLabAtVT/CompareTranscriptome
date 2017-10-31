#!/bin/sh
### 2017-04-11 by Jiyoung Lee
### To run fastq-dump to download SRR fils from SRA

WORKDIR=$(pwd)
OUTDIR=$2

mkdir -p $WORKDIR/raw_data/fastq/$OUTDIR
mkdir -p $WORKDIR/logs/fastq

### Setup the PATH environmental variable.
### This step is not necessary if the previous step has already copied files to $PROJ_PATH/software/bin folder.
# export PATH=$WORKDIR/software/sratoolkit.2.8.2-1-centos_linux64/bin:$PATH

echo "---------[Download SRR.fastq]---------"
for SRR in `cat $1`
do
	echo -e "`date +"%b%d|%T"`\t\c"; echo "$SRR"
	fastq-dump --gzip --split-3 $SRR --outdir $WORKDIR/raw_data/fastq/$OUTDIR >$WORKDIR/logs/fastq/$SRR.log 2>$WORKDIR/logs/fastq/$SRR.err
done
