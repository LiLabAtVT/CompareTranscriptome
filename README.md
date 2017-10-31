# Comparative Transcriptome Analysis
## 1. Introduction
### 1.1 Repository of scripts for comparative transcriptome analysis

Expression analysis is commonly used to understand the tissue or stress specificity of genes in large gene families. The goal of comparative transcriptome analysis is to identify conserved co-expressed genes in two or more species. To compare transcriptomes between any two species, we took a three step approach. 

1. Establish homologous relationships between proteins in the two species.
2. Identify expression data obtained from experiments that are performed under similar conditions or tissue types. 
3. Compare the expression patterns between the two data sets. 

In this protocol, we will compare published time course seed embryo expression data from Arabidopsis with data from the same tissue in soybean as a demonstration of  how to apply computational tools to comparative transcriptome analysis.  

### 1.2 How to obtain example scripts and datasets 

To use scripts and files provided in this protocol, you can have two options. 

1. Intall `git` on your Linux or Mac OX and then clone this repository with HTTPS or with SSH on Terminal using `git clone [URL] [Folder_Name]`
```
git clone https://github.com/LiLabAtVT/CompareTranscriptomeMIMB.git ATH_GMA
```
2. Download a Zip file of this repository from the green menu, 'Clone or download' on the top-right side of this main page, and move files to where you want to locate them.

For a project folder, we named `ATH_GMA` as default name, and used it as the working directory. If you replace the project folder name with another name, please be careful to replace it for scripts or paths. 


## 2. Materials 
### 2.1 Set up folder structure for data analysis.

For effective data management, we suggest you to create a folder structure using the script `Section2.1_setup_directory.sh` on `ATH_GMA` folder.
```bash 
$ cd ATH_GMA	#If you are not in ATH_GMA folder, move to the folder. 
$ sh ./scripts/Section2.1_setup_directory.sh
```
After executing this script, you will get a folder structure like a figure below. 
![alt text](https://github.com/LiLabAtVT/CompareTranscriptomeMIMB/raw/master/docs/figures/Screen%20Shot%202017-10-30%20at%2016.45.14.png)


### 2.2 Software installation

This step is for installation of majority of software on Linux for this protocol. By running `Section2.2_download_softwares.sh`, you can automatically download and install softwars for this protocol. If softwares are successfully installed, you will get a tool version information for each tool. If you have an error for one software, you don't need to run the script again and just install the software with error manually. Also you can change versions of tool by editing the script or by providing different repository URLs. You can find more descriptions from the manuscript. 
```bash
$ cd ATH_GMA	#If you are already in ATH_GMA folder, do not type it. 
$ sh ./scripts/Section2.2_download_softwares.sh
$ cat ./scripts/Section2.2_download_softwares.sh
```

### 2.3 Install R, DESeq2, and edgeR packages for RNA-Seq data analysis.

In addition to softwares from the previous step, we will use R that is a programing language and environment for statistical data analysis.

If you don't have R in your system, you need to install R before moving to next step. If you are not sure whether your system has `R`, you can simply type `R` on your terminal. 
```bash
$R

R version 3.3.1 (2016-06-21) -- "Bug in Your Hair"
Copyright (C) 2016 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)
...
Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

>
```
To install R packages for this analysis, we recommend to type commends in R interactive console rather than running this Rscript. You can start `R` in a Linux terminal or Terminal for OSX by simply typing `R`. 
```bash
$ cd ATH_GMA	#If you are already in ATH_GMA folder, do not type it. 
$ cat ./scripts/Section2.3_install_r_packages.R
$ R
```
Commands preceded by `$` are executed under a Linux terminal, and commands preceded by `>` are executed under the R environments.

We will install R packages though packager installers like bioconductor or from local file. To install the OrthoClust package, the user should download the script for the OrthoClust package. If you clone or download the software directory of this GitHub repository, you can find `OrthoClust_1.0.tar.gz` for the OrthoClust package. 

```
> source('https://bioconductor.org/biocLite.R')
> biocLite('DESeq2')
> biocLite('edgeR')
> install.packages('igraph')
# During process, you may need to choose a CRAN mirror close to you, if you haven't done it.
> setwd("./software")
> install.packages("OrthoClust_1.0.tar.gz", repos=NULL, type="source")
```

This process could take for a while and return wordy messages. If you see `* DONE (PACKAGE NAME)`, it usually means the package installed successfully. By typing `library('PACKAGE NAME')`, you can verify that the installation process has been completed. `OrthoClust` requires `igraph` package as one of its dependencies.


### 2.4 Download protein and genome sequences for Arabidopsis and soybean.

The Araport and DOE phytozome database require free registration to access and download data.

1. Arabidopsis: Araport web site (www.araport.org) 
- genomic sequences: TAIR10_Chr.all.fasta.gz on https://www.araport.org/downloads/TAIR10_genome_release/assembly
- protein-coding sequences: Araport11_genes.201606.pep.fasta.gz on https://www.araport.org/downloads/Araport11_latest
- gene annotation: Araport11_GFF3_genes_transposons.201606.gtf.gz on https://www.araport.org/downloads/Araport11_latest

2. soybean: DOE phytozome database (https://phytozome.jgi.doe.gov/pz/portal.html#!bulk?org=Org_Gmax)
- genomic sequences: Gmax_275_v2.0.fa.gz
- protein-coding sequences: Gmax_275_Wm82.a2.v1.protein.fa.gz
- gene annotation: Gmax_275_Wm82.a2.v1.gene_exons.gff3.gz

After downloading files, you need to locate them into `ATH_GMA/raw_data` folder to let scripts utilize these files. Then, you can use gunzip command such as `gunzip [gz compressed file]` to de-compress '*.gz' files. 
```bash
$ cd ATH_GMA	#If you are already in ATH_GMA folder, do not type it. 
$ cd raw_data
$ gunzip Araport11_genes.201606.pep.fasta.gz
$ gunzip Gmax_275_Wm82.a2.v1.protein.fa.gz
$ cd ..			#To move back to ATH_GMA folder, you can type it.
```

### 2.5 Download raw data from published RNA-Seq experiments

To download raw sequencing data using SRR accession number, we can you fastq-dump. If you already installed fastq-dump from `2.2 Software installation`, you can simply type `fastq-dump` on a terminal. You can find fastq-dump on `ATH_GMA/software/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump`. 

An example below is to download raw data of SRR2927328 to `./raw_data` folder. 
```bash
$ cd ATH_GMA	#If you are already in ATH_GMA folder, do not type it. 
$ fastq-dump --split-3 SRR2927328 –outdir ./raw_data
```

Since we have multiple SRR files, we provide `ATH_GMA/scripts/Section2.5_download_fastq.sh` with text files with a list of SRR accession numbers. You can have a test run with `ATH_GMA/raw_data/PRJNAtest.txt` to test the execution time for downloading one file. Downloaded file will be located in `ATH_GMA/raw_data/fastq/[GivenName_by_User]`. 
```bash
$ cd ATH_GMA
$ sh ./scripts/Section2.5_download_fastq.sh ./raw_data/PRJNAtest.txt TEST
$ sh ./scripts/Section2.5_download_fastq.sh ./raw_data/PRJNA301162.txt ATH
$ sh ./scripts/Section2.5_download_fastq.sh ./raw_data/PRJNA197379.txt GMA
```


## 3. Methods
### 3.1 Identification of homologous pairs using BLAST.
### 3.2 Obtaining reciprocal best hit (RBH) genes
### 3.3 Gene expression data processing. 
### 3.4 OrthoClust analysis.
### 3.5 Visualization of OrthoClust results
## 4. Notes
### 4.2 Obtaining one way best hit genes from each species

As we explained in Notes 4.4 section, we provide scripts and sample datasets to identify k-best-hit genes from two species. 
```bash
$ cd ATH_GMA
$ sh ./scripts/Section4.4_k-BestHits.sh
```

A script, `Section4.4_k-BestHits.sh` is an example for k = 5. By modifying or adding commands, you can try different `k` values. For example, in this command line, `python ../scripts/OrthologousGenes_OneWayTopNBestHit.py 5 ARATH2GLYMA.BLAST_ARATH2GLYMA.subset.txt`, you can put positive integer number (such as 2, 3, 4) instead of 5.   

This shell script, `Section4.4_k-BestHits.sh` calls a python script `OrthologousGenes_OneWayTopNBestHit.py`. This python script requires input files located in processed_data folder. 

To generate these input files, 2e separated BLAST results into 1) ARATH2GLYMA.BLAST_ARATH2GLYMA.txt generated from Arabidopsis queries to soybean DB and 2) ARATH2GLYMA.BLAST_GLYMA2ARATH.txt generated from soybean queries to Arabidopsis DB. As sample files, we uploaded 2,000 lines as subsets of input files. 

