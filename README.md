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

### 2.4 Download protein and genome sequences for Arabidopsis and soybean.

### 2.5 Download raw data from published RNA-Seq experiments

## 3. Methods
### 3.1 Identification of homologous pairs using BLAST.
### 3.2 Obtaining reciprocal best hit (RBH) genes
### 3.3 Gene expression data processing. 
### 3.4 OrthoClust analysis.
### 3.5 Visualization of OrthoClust results
