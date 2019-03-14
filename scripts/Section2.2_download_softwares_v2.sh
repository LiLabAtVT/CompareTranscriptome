#!/bin/sh
### 2017-08-18 by Jiyoung Lee

### Scripts for section 2.2 to section 2.4 
### Execute this script from current work directory (default: ATH_GMA)
### TO execute this script, run the following command:
### sh ./scripts/Section2.2_download_softwares.sh

### Setting up the PATH for installed softwares
workdir=$(pwd)
softwarepath=$workdir/software/bin/
echo $softwarepath
export PATH=$PATH:$softwarepath

### Installing BLAST  ===============================================================
cd software

echo -e "\n\n`date +"%b%d|%T"`\t\c";
echo "======== Installing BLAST ========"
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz
tar -xzf ncbi-blast-2.6.0+-x64-linux.tar.gz
./ncbi-blast-2.6.0+/bin/makeblastdb -h	
cp ./ncbi-blast-2.6.0+/bin/* $softwarepath

### Installing fastq-dump  ===============================================================
echo -e "\n\n`date +"%b%d|%T"`\t\c";
echo "======== Installing fastq-dump ========"
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar -xzf sratoolkit.current-centos_linux64.tar.gz
./sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump --version

find ./sratoolkit.2.8.2-1-centos_linux64/bin/ -type f -exec cp {}  $softwarepath \;
 
### Installing STAR ======================================================================
echo -e "\n\n`date +"%b%d|%T"`\t\c";
echo "======== Installing STAR ========"
wget https://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz
tar -xzf 2.5.2b.tar.gz
./STAR-2.5.2b/bin/Linux_x86_64_static/STAR --version
cp ./STAR-2.5.2b/bin/Linux_x86_64_static/* $softwarepath
 
### Installing featureCounts =============================================================
echo -e "\n\n`date +"%b%d|%T"`\t\c";
echo "======== Installing featureCounts ========"
wget https://sourceforge.net/projects/subread/files/subread-1.5.1/subread-1.5.1-Linux-x86_64.tar.gz/download 
tar -zxf download
./subread-1.5.1-Linux-x86_64/bin/featureCounts -v
find ./subread-1.5.1-Linux-x86_64/bin/ -type f -exec cp {}  $softwarepath \;

### Check R version  =============================================================
echo -e "\n\n`date +"%b%d|%T"`\t\c";
echo "======== Checking R, Rscript versions ========"
R --version
Rscript --version

### If R is not installed, use the following command to install R
# sudo yum install R

### Check Python version  =============================================================
echo -e "\n\n`date +"%b%d|%T"`\t\c";
echo "======== Checking Python version ========"
python -V

### If python is not installed, use the following command to install python
# sudo yum install python

