# Comparative Transcriptome Analysis
## 1. Introduction
### 1.1 Repository of scripts for comparative transcriptome analysis

Expression analysis is commonly used to understand the tissue or stress specificity of genes in large gene families. The goal of comparative transcriptome analysis is to identify conserved co-expressed genes in two or more species. In order to compare transcriptomes between any two species, we took a three-step approach: 

1. Establish homologous relationships between proteins in the two species.
2. Identify expression data obtained from experiments that are performed under similar conditions or tissue types. 
3. Compare the expression patterns between the two data sets. 

In this protocol, we will compare published time course seed embryo expression data from Arabidopsis with data from the same tissue in soybean as a demonstration of how to apply computational tools to comparative transcriptome analysis.  

### 1.2 How to obtain example scripts and datasets 

To use scripts and files provided in this protocol, there are two options. 

1. Intall `git` in Linux or Mac OS and then clone this repository with HTTPS or with SSH in a Terminal using `git clone [URL] [Folder_Name]`
```
git clone https://github.com/LiLabAtVT/CompareTranscriptomeMIMB.git ATH_GMA
```
2. Download a Zip file of this repository from the green menu, 'Clone or download', on the top-right side of this main page, and move files to where the user want to locate them.

For a project folder, we named `ATH_GMA` as default name, and used it as the working directory. When the project folder name is replaced with another name, please replace it for scripts and paths in the github files. 

![alt text](https://github.com/LiLabAtVT/CompareTranscriptomeMIMB/raw/master/docs/figures/Clone_Download_GitRep.png)

### 1.3 LINUX commands
If the user is a beginner at Linux or OS X Terminal, we strongly suggest reviewing basic Linux commands before or during further practices.

## 2. Materials 
### 2.1 Set up folder structure for data analysis

For effective data management, we suggest creating a folder structure using the script `Section2.1_setup_directory.sh` on `ATH_GMA` folder.
```bash 
$ cd ATH_GMA	# Please move to this folder to run the script below. 
$ sh ./scripts/Section2.1_setup_directory.sh
```

### The initial folder structure
After executing this script, the folder structure like a figure below will be created. 
![alt text](https://github.com/LiLabAtVT/CompareTranscriptomeMIMB/raw/master/docs/figures/Inital_FolderStructure.png)


### 2.2 Software installation

This step covers installation of majority of software applying within this protocol. `Section2.2_download_software.sh` allows users automatically download and install software. 

```bash
$ cd ATH_GMA	# Please move to this folder to run the script below. 
$ sh ./scripts/Section2.2_download_software.sh
$ cat ./scripts/Section2.2_download_software.sh
```
If software is successfully installed, the version information for each tool will be displayed. If errors appear rather than the version information for a tool, the tool needs to be installed manually. The versions of software can be changed with different repository URLs for the software. 

### 2.3 Install R, DESeq2, and edgeR packages for RNA-Seq data analysis

In addition to the software from the previous step, we will use `**R**`, which is a programing language and environment for statistical data analysis. If `R` is not installed on a system, it will be had to be installed before moving to next step. In order to determine whether the system has `R`, users can simply type R into the Terminal. 

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
To install R packages for this protocol, we recommend typing commands in a R interactive console rather than running this Rscript in a Terminal. R can be started in a Linux or OSX Terminal by typing `R`. 

```bash
$ cd ATH_GMA	# Please move to this folder to run the script below. 
$ cat ./scripts/Section2.3_install_r_packages.R
$ R
```
In this protocol, commands preceded by `$` are executed under a Linux Terminal, and commands preceded by `>` are executed under a R console.

All R packages for this protocol will be able to be installed though package installers like Bioconductor except `OrthoClust`. `OrthoClust is required to be installed from its package file, and `OrthoClust_1.0.tar.gz` is included in `software` directory of this GitHub repository. 
Some packages require other packages as their dependencies such as `OrthoClust` requires `igraph`. We recommend allowing R Package Installer to install dependences together.  

```
> source('https://bioconductor.org/biocLite.R')
> biocLite('DESeq2')
> biocLite('edgeR')
> install.packages('igraph')
# During this installing process, R may request to choose a CRAN mirror from the user.
> setwd("./software")
> install.packages("OrthoClust_1.0.tar.gz", repos=NULL, type="source")
```

This process could take a while and return wordy messages. When a package has been installed successfully, `* DONE (PACKAGE NAME)` usually appears on a console. The installation process can be verified by typing `library('PACKAGE NAME')`. 


### 2.4 Download protein and genome sequences for Arabidopsis and soybean

The Araport and DOE phytozome database both require free registration to access and download their data.

1. Arabidopsis: Araport web site (www.araport.org) 
   - genomic sequences: TAIR10_Chr.all.fasta.gz on https://www.araport.org/downloads/TAIR10_genome_release/assembly
   - protein-coding sequences: Araport11_genes.201606.pep.fasta.gz on https://www.araport.org/downloads/Araport11_latest
   - gene annotation: Araport11_GFF3_genes_transposons.201606.gtf.gz on https://www.araport.org/downloads/Araport11_latest

2. soybean: DOE phytozome database (https://phytozome.jgi.doe.gov/pz/portal.html#!bulk?org=Org_Gmax)
   - genomic sequences: Gmax_275_v2.0.fa.gz
   - protein-coding sequences: Gmax_275_Wm82.a2.v1.protein.fa.gz
   - gene annotation: Gmax_275_Wm82.a2.v1.gene_exons.gff3.gz

After downloading sequence files and gene annotation files, these files are needed to be moved into the `ATH_GMA/raw_data` directory to let future scripts to work on these files. Then, compressed '*.gz' files can be de-compressed with gunzip command such as `gunzip [gz compressed file]`. 

```bash
$ cd ATH_GMA 
$ cd raw_data
$ gunzip Araport11_genes.201606.pep.fasta.gz
$ gunzip Gmax_275_Wm82.a2.v1.protein.fa.gz
$ cd ..			# Move back to the ATH_GMA directory.
```

### 2.5 Download raw data from published RNA-Seq experiments

`fastq-dump` is a tool from the Sequence Read Archive (SRA) toolkit, and allows the user to download raw sequencing data with SRR accession ID from SRA. If fastq-dump is already installed from the `2.2 Software installation` section, then typing `fastq-dump` in a Terminal shows its usage and version information, and it is located in the `ATH_GMA/software/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump` directory. 

The example below is to download a sequencing file (accession ID: SRR2927328) to the `./raw_data` directory. 
```bash
$ cd ATH_GMA	# Please move to this folder to run the script below. 
$ fastq-dump --split-3 SRR2927328 –outdir ./raw_data
```

We provide the `Section2.5_download_fastq.sh` script and two text files with a list of SRR accession IDs (PRJNA301162.txt for Arabidopsis and PRJNA197379.txt for soybean) for allowing easy download with multiple SRR files. Since downloading sped would be various depending on network environments, the approximated execution time can be estimated by running the script with a test file, `PRJNAtest.txt`. Obtained files will be located in the `ATH_GMA/raw_data/fastq/[TEST|ATH|GMA]` directory.
 
```bash
$ cd ATH_GMA
$ sh ./scripts/Section2.5_download_fastq.sh ./raw_data/PRJNAtest.txt TEST
$ sh ./scripts/Section2.5_download_fastq.sh ./raw_data/PRJNA301162.txt ATH
$ sh ./scripts/Section2.5_download_fastq.sh ./raw_data/PRJNA197379.txt GMA
```

### The second folder structure
Successfully completed steps in **`2. Materials`** section result in a folder structure shown in the figure below. 

![alt text](https://github.com/LiLabAtVT/CompareTranscriptomeMIMB/raw/master/docs/figures/Materials_FolderStructure.png)

## 3. Methods

### A workflow of comparative transcriptome analysis between soybean and Arabidopsis 
A workflow (See Figure 2 in the book chapter for this protocol) is composed of three major parts: 1) identification of orthologous pairs between two species using BLAST, 2) RNA-seq analysis to obtain co-expression networks, and 3) running OrthoClust to cluster genes with orthologous relations. 

### Setting up the PATH for installed software.
After installing the software from the `2.2 Software installation` section, we recommend the user either to add paths of the installed software in the `bashrc` or `bash_profile` file in the system for continues uses or to export their paths on every use.
 
```bash
# command lines below from line 9~13 in Section2.2_download_software_v2.sh. 
### Setting up the PATH for installed software
$ cd ATH_GMA
$ workdir=$(pwd)
$ softwarepath=$workdir/software/bin/
$ echo $softwarepath
$ export PATH=$PATH:$softwarepath
```

### 3.1 Identification of homologous pairs using BLAST

This section aims to identify homologous protein pairs between two species by using BLAST, which provide sequence alignments results. To ease this step, we provide the `Section3.2.1_BLAST.sh` with the following usage below.
 
```bash
$ cd ATH_GMA 
$ sh ./scripts/Section3.2.1_BLAST.sh
```

The `Section3.2.1_BLAST.sh` script is composed of three steps, and the same results can be achieved by following three instructions below. Please pay attention to locate the right working directory for each step.

1. Prepare BLAST input data by merge two protein sequence files into one file.
```bash
$ cd ATH_GMA 
# cat [first_file] [second_file] > [merged_file]
$ cat ./raw_data/Araport11.pep.fasta ./raw_data/GLYMA2.pep.fasta > ./processed_data/ATHGMA.pep.fasta
```

2. Build BLAST database with `makeblastdb`. 
```bash
$ cd ATH_GMA/processed_data	# This step is performed in "processed_data" directory.
$ makeblastdb -in ATHGMA.pep.fasta \
              -out ATHGMA.pep.blastdb \
              -dbtype prot \
              -logfile makeblastdb.log
```

3. Perform BLAST analysis for protein sequences with `blastp`.
```bash
$ cd ATH_GMA/processed_data	# This step is performed in "processed_data" directory.
$ blastp -evalue 0.00001 \
       -outfmt 6 -db ATHGMAX.pep.blastdb \
       -query ATHGMA.pep.fasta > ATHGMA.pep.blastout 
```

For more information about BLAST options, please refer to `3.2.1 Identification of homologous pairs using BLAST` in the book chapter or BLAST user manuals from NCBI.

### 3.2 identification of the reciprocal best hit (RBH) genes

The RBH genes are pairs of genes between two species with the best BLAST hit mutually among all BLAST results. The shell script ` Section3.2.2_RBH.sh` calls the python script` ReciprocalBlastHit.py` inside to extract RBH pairs from BLAST results from the previous section.

```bash
$ cd ATH_GMA 
$ sh ./scripts/Section3.2.2_RBH.sh
```
OR 
```bash
$ cd ATH_GMA/processed_data	# This step is performed in "processed_data" directory.
$ python ../scripts/ReciprocalBlastHit.py ATHGMA.pep.blastout ARATH GLYMA ARATH2GLYMA.RBH.txt 
```

A subset of RBH genes is provided in  [ARATH2GLYMA.RBH.subset.txt](https://raw.githubusercontent.com/LiLabAtVT/CompareTranscriptomeMIMB/master/processed_data/ARATH2GLYMA.RBH.subset.txt). 

The table below shows the first 10 lines of the `ARATH2GLYMA.RBH.subset.txt` file as an example. E-value of each gene is extracted from BLAST results.

| Soybean gene | Arabidopsis gene | E-value of Soybean gene | E-value of Arabidopsis gene |
| :---: | :---: | :---: | :---: |
| Glyma.01G001300	| AT2G07050 | 0.0 | 0.0 |
| Glyma.01G005800 | AT4G29310 | 0.0 | 0.0 |
| Glyma.01G006100 | AT4G26300 | 0.0 | 0.0 |
| Glyma.01G010100 | AT1G32090 | 0.0 | 0.0 |
| Glyma.01G015400 | AT2G35470 | 6.26e-16 | 4.31e-18 |
| Glyma.01G019400 | AT5G65670 | 1.9e-149 | 5.47e-145 |
| Glyma.01G019700 | AT5G65640 | 5.48e-87 | 3.47e-87 |
| Glyma.01G021300 | AT4G38040 | 0.0 | 0.0 |
| Glyma.01G021400 | AT2G22600 | 3.69e-115 | 1.55e-117 |
| Glyma.01G022500 | AT5G10510 | 0.0 | 0.0 |

RBH genes are widely used in comparative genomic analysis, and other methods can be used to identify homologous genes for downstream analysis. Please see [Note 4.2](https://github.com/LiLabAtVT/CompareTranscriptomeMIMB#42-obtaining-one-way-best-hit-genes-from-each-species) for a description of one way best hit method.

### 3.3 Gene expression data processing

This section covers RNA-seq data analysis to generate co-expression networks from gene expression data. It can be composed of three parts with following five steps: A) read mapping (step 1 and step 2), B) read counting (step 3) and C) FPKM calculation (step 4 and step 5). 

Read mapping is a process to map RNA-seq reads to the respective reference genome. For more details of STAR options, please refer to `3.3 Gene expression data processing` in the book chapter of this protocol.

**1. Create genome index by STAR before read mapping**

   Read mapping tools require creating genome index files to increase mapping speed. Indexing genome only needs to be performed once for the first time to use STAR with the respective references sequences for each species. 

   ```bash
   $ cd ATH_GMA
   $ sh ./scripts/Section3.3.Step1.MakeIndex.sh
   ```
   `Section3.3.Step1.MakeIndex.sh` contains two command lines to build index for Arabidopsis and soybean. Below is an example for Arabidopsis genome. 
   ```bash
   $ WORKDIR=$(pwd)
   $ IDX=$WORKDIR/raw_data/ATH_STAR-2.5.2b_index
   $ GNM=$WORKDIR/raw_data/TAIR10_Chr.all.fasta
   $ GTF=$WORKDIR/raw_data/Araport11_GFF3_genes_transposons.201606.gtf
   $ STAR	--runMode genomeGenerate \
			--genomeDir $IDX \
			--genomeFastaFiles $GNM \
			--sjdbGTFfile $GTF
   ```
   When indexing is successfully finished, index files can be found in the `ATH_GMA/raw_data/[ATH|GMA]_STAR-2.5.2b_index` folder. 
   ```
   # These are example files from indexing by STAR with Arabidopsis genome.
   $ cd ATH_GMA/raw_data/ATH_STAR-2.5.2b_index
   $ ls -sh	# This command prints out file names and their approximate human-readable file sizes. (zero size here does not mean empty.) 
   total 2.8G
      0 chrLength.txt      9.8M exonGeTrInfo.tab     0 genomeParameters.txt  3.3M sjdbList.fromGTF.out.tab
      0 chrNameLength.txt  4.0M exonInfo.tab      1.2G SA                    3.3M sjdbList.out.tab
      0 chrName.txt        512K geneInfo.tab      1.5G SAindex               3.0M transcriptInfo.tab
      0 chrStart.txt       142M Genome            3.5M sjdbInfo.txt
   ```

**2. Read mapping by STAR**

   The `Section3.3.Step2.Mapping.[ATH|GMA].sh` script covers read mapping for all sequencing reads files in `ATH_GMA/raw_data/fastq/[ATH|GMA]/` directories. 

   ```bash
   $ cd ATH_GMA
   $ sh ./scripts/Section3.3.Step2.Mapping.ATH.sh
   $ sh ./scripts/Section3.3.Step2.Mapping.GMA.sh
   ```
   The following command below is an example for SRR2927328 (*_1 and *_2, paired-end reads files have two inputs).
   ```bash
   $ STAR	--genomeDir $IDX \
			--readFilesIn $WORKDIR/raw_data/SRR2927328_1.fastq.gz   $WORKDIR/raw_data/SRR2927328_2.fastq.gz \
			--outFileNamePrefix $WORKDIR/processed_data/bam/SRR2927328/SRR2927328 \
			--outSAMtype BAM SortedByCoordinate
   ```

   STAR mapping results are generated in the `ATH_GMA/processed_data/bam/[SRR ID]` directory. Mapping statistics can be checked from `[SRR ID]Log.final.out` files. 
   ```
   $ cd ATH_GMA/processed_data/bam/SRR2927328
   $ ls -sh
   total 2.1G
   2.1G SRR2927328Aligned.sortedByCoord.out.bam  256K SRR2927328Log.out           3.8M SRR2927328SJ.out.tab
      0 SRR2927328Log.final.out                     0 SRR2927328Log.progress.out
   $ cat SRR2927328Log.final.out
   ```

**3. Read counting by featureCounts**

   With mapping results from the previous step, FeatureCounts will calculate the number of mapped reads on the corresponding genome. The `Section3.3.Step3.ReadCount.[ATH|GMA].sh` scripts are designed to perform featureCounts for all mapping results under `ATH_GMA/processed_data/bam/[SRR ID]` folders. 

   ```bash
   $ cd ATH_GMA
   $ sh ./scripts/Section3.3.Step3.ReadCount.ATH.sh
   $ sh ./scripts/Section3.3.Step3.ReadCount.GMA.sh
   ```
   The following command lines below are requirements to process read mapping results of SRR2927328 as an example.
   ```bash
   $ WORKDIR=$(pwd)
   $ GTF=$WORKDIR/raw_data/Araport11_GFF3_genes_transposons.201606.gtf
   $ BAM=$WORKDIR/processed_data/bam
   $ RC=$WORKDIR/processed_data/rc
   $ featureCounts	-t exon \
				-g gene_id \
				-p \
				-a $GTF \
				-o $RC/SRR2927328.readcount.txt \
				$BAM/SRR2927328/SRR2927328Aligned.sortedByCoord.out.bam
   ```

**4. FPKM calculation using DESeq2 and edgeR**

   We provide the unified Rscript, `ATH_GMA/scripts/Section3.3.Step4.FPKM.R` and csv files for replicate structure of the samples (PRJNA301162.csv for Arabidopsis and PRJNA197379.csv for soybean). The Rscript is composed of multiple functional blocks with detailed descriptions of input requirements, contents of outputs, and functions of each step. It generates several out files such as a result of differentially expressed genes analysis, normalized FPKM values for each SRR ID, and averaged FPKM values for each condition group. For more details, please look at the script, `ATH_GMA/scripts/Section3.3.Step4.FPKM.R`.

   ```bash
   $ cd ATH_GMA
   $ Rscript ./scripts/Section3.3.Step4.FPKM.R ./processed_data/fpkm/ATH	# for Arabidopsis data
   $ Rscript ./scripts/Section3.3.Step4.FPKM.R ./processed_data/fpkm/GMA	# for soybean data
   ```

**5. Generation of co-expression networks from gene expression profile**

     A co-expression network can be generated from correlation values between gene expression levels of averaged FPKM values by calculating Pearson correlation coefficient values and their p-values. Gene pairs with p-value < 0.001 and Pearson Correlation Coefficient > 0.99 are selected to generate a co-expression network. Pairs of genes are used as node, and correlation values of gene pairs are used for their edge value in a co-expression network. 

     A Rscript we provide, `ATH_GMA/scripts/Section3.3.Step5_FPKM2NETWORK.R`, converts a gene expression matrix to a co-expression network. This R script is written for an input file, `GMX_FPKM.csv`. To generate a co-expression network for Arabidopsis data, `GMX_FPKM.csv` can be replaced with `ATH_FPKM.csv` in the R script. For more details, please read descriptions in the script, `ATH_GMA/scripts/Section3.3.Step5_FPKM2NETWORK.R`. 
 
     ```bash 
     $ cd ATH_GMA
     $ Rscript ./scripts/Section3.3.Step5_FPKM2NETWORK.R
     ```

### 3.4 OrthoClust analysis

OrthoClust is a published integrated clustering method by incorporating co-association data of multiple species with orthologous relation between species. In order to cluster genes from Arabidopsis and soybean based on their co-expression data and RBH pairs, OrthoClust analysis will be performed in this section. 

With output files from the previous section, we can run `OrthoClust` and save its results into file for further analysis using a Rscript `ATH_GMA/scripts/Section3.4.Step1_OrthoClust.R`.
To perform OthoClust analysis, we require three input data files: 1) the gene co-expression network from soybean; 2) the gene co-expression network from Arabidopsis; and 3) the orthologous gene pairs between two species. **Table 2** in the book chapter shows examples of these input data files. 

```bash
$ cd ATH_GMA
$ Rscript ./scripts/Section3.4.Step1_OrthoClust.R
```

### 3.5 Visualization of OrthoClust results by Cytoscape

Cytoscape is a platform to visualize networks and to incorporate open source applications for network analysis. Here we use Cytoscape to visualize OrthoClust results and there are three steps. 

**1. Preparation of Cytoscape input files** 
     To prepare visualization of OrthoClust results, we can choose one module among several hundred modules. `Section3.4.Step2_CytoscapeInput.R` will show demonstrations for extracting genes in an interesting module and making plots with these genes.

     We suggest running commands of `Section3.4.Step2_CytoscapeInput.R` in the R environments such as R console rather than executing the Rscript in a Terminal. This is because OrthoClust randomly assigns genes into modules, generated results for each run from the user are different from what the book chapter reports. As an example, we chose module 8 (`ModuleName=8`) on the 39th line of `Section3.4.Step2_CytoscapeInput.R` from OrthoClust result. However, gene assignment will be different for every run, so the user needs to change this module number according to results from `print(ModulesOfInterest)`, the 38th line of the Rscript. Also, running R commands in R console allow the user to explore different module numbers. `Orthoclust_Results_Summary.csv` file from the previous section provide module IDs and the number of genes for each module. We suggest checking the file to choose interesting modules to explore in this section.

     ```
     line 37: ModulesOfInterest= ortho_module_list_sumOrdered[c(1:10),1:6]
     line 38: print(ModulesOfInterest)
     line 39: ModuleName=8	#This integer number can be changed according to the user’s interest. 
     ```
     R command lines in this section can be run in the `ATH_GMA`, the working directory by pasting commands from `Section3.4.Step2_CytoscapeInput.R` on R environments.
     ```bash
     $ cd ATH_GMA
     $ R
     > 
     ```
     For quick process, running R script in a Terminal also works. 
     ```bash
     $ cd ATH_GMA
     $ Rscript ./scripts/Section3.4.Step2_CytoscapeInput.R
     ```
**2. Gene expression patterns for modules**
     Among results from `Section3.4.Step2_CytoscapeInput.R`, there will be a PDF named `ExpressionProfile_Module8.pdf` or `ExpressionProfile_Module[Chosen Integer].pdf` for a selected module from the previous step. This figure shows expression levels of genes from the module with averaged expression value along different time points, and it would be similar with
**Figure 4. Expression plots of genes from Arabidopsis and soybea** in the book chapter. This Rscrip can be modified and improved according to the need and interests of the user. 
     
**3. Network visualization using Cytoscape.**
     Please refer to the book chapter subtitled, **3.4.3 Visualization of OrthoClust results as a network** for this visualization step using Cytoscape.


## 4. Notes
### 4.2 Obtaining one way best hit genes from each species

As we explained in Notes 4 in the book chapter, we provide a script and sample datasets to identify k-best-hit genes from two species. 

```bash
$ cd ATH_GMA
$ sh ./scripts/Section4.4_k-BestHits.sh
```

This shell script, `Section4.4_k-BestHits.sh` calls a python script `OrthologousGenes_OneWayTopNBestHit.py`. The shell script provides an example to identify five best hit genes (k=5). By modifying or adding commands in the shell script, different `k` values can be used. For example, in the command line, `python ../scripts/OrthologousGenes_OneWayTopNBestHit.py 5 ARATH2GLYMA.BLAST_ARATH2GLYMA.subset.txt`, other positive integer numbers (such 2, 3, 4) can be replaced with `5`.  

The python script requires input files located in processed_data folder. To generate these input files for the python script, separated BLAST results are needed 1) `ARATH2GLYMA.BLAST_ARATH2GLYMA.txt` generated from Arabidopsis queries to soybean DB and 2) `ARATH2GLYMA.BLAST_GLYMA2ARATH.txt` generated from soybean queries to Arabidopsis DB. As sample files, we uploaded 2,000 lines as subsets of input files. 


