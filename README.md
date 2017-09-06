# Comparative Transcriptome Analysis
### 1. Introduction
#### 1.1 Repository of scripts for comparative transcriptome analysis

Expression analysis is commonly used to understand the tissue or stress specificity of genes in large gene families. The goal of comparative transcriptome analysis is to identify conserved co-expressed genes in two or more species. To compare transcriptomes between any two species, we took a three step approach. 

1. Establish homologous relationships between proteins in the two species.
2. Identify expression data obtained from experiments that are performed under similar conditions or tissue types. 
3. Compare the expression patterns between the two data sets. 

In this protocol, we will compare published time course seed embryo expression data from Arabidopsis with data from the same tissue in soybean as a demonstration of  how to apply computational tools to comparative transcriptome analysis.  

To use the scripts provided in this protocol, you need to install git, and clone the repository.

```
git clone https://github.com/LiLabAtVT/CompareTranscriptomeMIMB.git ATH_GMA
```

## 2. Materials 
### 2.1 Set up folder structure for data analysis.

```bash 
$ cd ATH_GMA
$ sh ./scripts/Section2.1_setup_directory.sh
```

### 2.2 Software installation

```bash
$ cd ATH_GMA
$ sh ./scripts/Section2.2_download_softwares.sh
```

### 2.3 Install R, DESeq2, and edgeR packages for RNA-Seq data analysis.

To install packages for this analysis, we recommend to type commends in R interactive console rather than running this Rscript. You can start `R` in a Linux terminal or Terminal for OSX by simply typing `R`. 
```bash
$ cd ATH_GMA
$ cat ./scripts/Section2.3_install_r_packages.R
$ R
```
Commands preceded by `$` are executed under a Linux terminal, and commands preceded by `>` are executed under the R environments.

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

### 2.5 Download raw data from published RNA-Seq experiments

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

