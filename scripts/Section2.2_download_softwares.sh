#!/bin/sh
### 2017-08-18 by Jiyoung Lee

### Scripts for section 2.2 to section 2.4 
### Execute this script from current work directory (default: ATH_GMA)
### TO execute this script, run the following command:
### sh ./scripts/Section2.2_download_softwares.sh

### Installing BLAST  ===============================================================
cd software

echo -e "`date +"%b%d|%T"`\t\c";
echo "Installing BLAST"
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz
tar -xzf ncbi-blast-2.6.0+-x64-linux.tar.gz
cp ./ncbi-blast-2.6.0+/bin/* $softwarepath

### Installing fastq-dump  ===============================================================
echo -e "`date +"%b%d|%T"`\t\c";
echo "Installing fastq-dump"
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
tar -xzf sratoolkit.current-centos_linux64.tar.gz
./sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump --version

cp ./sratoolkit.2.8.2-1-centos_linux64/bin/* $softwarepath
 
### Installing STAR ======================================================================
echo -e "`date +"%b%d|%T"`\t\c";
echo "Installing STAR"
wget https://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz
tar -xzf 2.5.2b.tar.gz
./STAR-2.5.2b/bin/Linux_x86_64_static/STAR –version
cp ./STAR-2.5.2b/bin/Linux_x86_64_static/* $softwarepath
 
### Installing featureCounts =============================================================
echo -e "`date +"%b%d|%T"`\t\c";
echo "Installing featureCounts"
wget https://sourceforge.net/projects/subread/files/subread-1.5.1/subread-1.5.1-Linux-x86_64.tar.gz/download 
tar -zxvf download
./subread-1.5.1-Linux-x86_64/bin/featureCounts -v
cp ./subread-1.5.1-Linux-x86_64/bin/* $softwarepath


### Check R version  =============================================================
R --version
Rscript --version

### If R is not installed, use the following command to install R
# sudo yum install R

### Check Python version  =============================================================
python -V

### If python is not installed, use the following command to install python
# sudo yum install python

