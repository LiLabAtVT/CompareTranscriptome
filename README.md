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

![alt text](https://github.com/LiLabAtVT/CompareTranscriptomeMIMB/raw/master/docs/figures/Clone_Download_GitRep.png)

### 1.3 LINUX commands
If you are a beginner of linux or OS X terminal, we strongly suggest you to review basic linux commands before or during further practices.

## 2. Materials 
### 2.1 Set up folder structure for data analysis.

For effective data management, we suggest you to create a folder structure using the script `Section2.1_setup_directory.sh` on `ATH_GMA` folder.
```bash 
$ cd ATH_GMA	#If you are not in ATH_GMA folder, move to the folder. 
$ sh ./scripts/Section2.1_setup_directory.sh
```

### The initial folder structure
After executing this script, you will get a folder structure like a figure below. 
![alt text](https://github.com/LiLabAtVT/CompareTranscriptomeMIMB/raw/master/docs/figures/Inital_FolderStructure.png)


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

### The second folder structure
After successfully completing steps in `2. Materials` section, you will get a folder structure like a figure below. We only put related scripts among all scripts for the protocol. 
![alt text](https://github.com/LiLabAtVT/CompareTranscriptomeMIMB/raw/master/docs/figures/Materials_FolderStructure.png)

## 3. Methods

### A workflow of comparative transcriptome analysis between soybean and Arabidopsis. 
A workflow (See Figure 2 in a book chapter) is composed of three major parts: 1) identification of ortholous pairs between two species using BLAST, 2) RNA-seq analysis to get co-expression networks, and 3) running OrthoClust to cluster genes with orthologous relations. Blue fonts indicate scripts and bold fonts refer softwares used in this workflow.

### c.f. Setting up the PATH for installed softwares.
After installing softwares, if you newly connect terminal but haven't add your PATH for installed softwares, we recommend you either to `export` PATH for the softwares again or to add the path in `bashrc` or `bash_profile` file in your system for permanent uses. 
```bash
# command lines below from line 9~13 in Section2.2_download_softwares_v2.sh. 
### Setting up the PATH for installed softwares
$ cd ATH_GMA
$ workdir=$(pwd)
$ softwarepath=$workdir/software/bin/
$ echo $softwarepath
$ export PATH=$PATH:$softwarepath
```

### 3.1 Identification of homologous pairs using BLAST.

Analysis in this section is composed of three steps. You can run `Section3.2.1_BLAST.sh` with command lines below. 
```bash
$ cd ATH_GMA	#If you are already in ATH_GMA folder, do not type it. 
$ sh ./scripts/Section3.2.1_BLAST.sh
```

Or you follow each step one by one. In this section, please pay attention to your location or working directory. For more details of BLAST options, please refer to `3.2.1 Identification of homologous pairs using BLAST.` in [our article in publisher name](https://github.com/LiLabAtVT/CompareTranscriptomeMIMB). 

1. Prepare data for BLAST analysis by merging two protein sequences from section 2.4 into one file `cat [first_file] [second_file] > [merged_file]`.  
```bash
$ cd ATH_GMA	#If you are already in ATH_GMA folder, do not type it. 
$ cat ./raw_data/Araport11.pep.fasta ./raw_data/GLYMA2.pep.fasta > ./processed_data/ATHGMA.pep.fasta
```

2. Build BLAST database with `makeblastdb`. 
```bash
$ cd processed_data	# this step is performed in "processed_data" folder
# Two commandlines below are doing the exactly same function. 
$ makeblastdb -in ATHGMA.pep.fasta -out ATHGMA.pep.blastdb -dbtype prot -logfile makeblastdb.log
$ makeblastdb -in ATHGMA.pep.fasta \
              -out ATHGMA.pep.blastdb \
              -dbtype prot \
              -logfile makeblastdb.log
```

3. Perform BLAST analysis for protein sequences with `blastp`.
```bash
$ cd processed_data	# this step is performed in "processed_data" folder
blastp -evalue 0.00001 \
       -outfmt 6 -db ATHGMAX.pep.blastdb \
       -query ATHGMA.pep.fasta > ATHGMA.pep.blastout 
```

### 3.2 Obtaining reciprocal best hit (RBH) genes

To identify homologous genes in two species, the BLAST results from protein sequence alignment were first parsed to identify the best BLAST hit for each soybean protein in the Arabidopsis protein lists. 

```bash
$ cd ATH_GMA	#If you are already in ATH_GMA folder, do not type it. 
$ sh ./scripts/Section3.2.2_RBH.sh
```
OR 
```bash
$ cd processed_data	# this step is performed in "processed_data" folder
$ python ../scripts/ReciprocalBlastHit.py ATHGMA.pep.blastout ARATH GLYMA ARATH2GLYMA.RBH.txt 
```

An example file [ARATH2GLYMA.RBH.subset.txt](https://raw.githubusercontent.com/LiLabAtVT/CompareTranscriptomeMIMB/master/processed_data/ARATH2GLYMA.RBH.subset.txt) of RBH genes is provided. The user can use this file to perform the following analysis without running the RBH script. 

Although RBH genes are widely used in comparative genomic analysis, other methods can be used to identify homologous genes for downstream analysis (see [Note 4.2](https://github.com/LiLabAtVT/CompareTranscriptomeMIMB#42-obtaining-one-way-best-hit-genes-from-each-species)). 

### 3.3 Gene expression data processing. 

This section is mainly composed of three steps: A) read mapping (step 1 and step 2); B) read counting (step 3) and C) FPKM calculation. Three steps will be performed with following sub-steps. 

Read mapping is a processto align RNA-seq reads to the respective reference genome. For more details of STAR options, please refer to `3.3 Gene expression data processing` in [our article in publisher name](https://github.com/LiLabAtVT/CompareTranscriptomeMIMB).

1. Create genome index by STAR before read mapping. 
This step is only required to perform once with the respective references sequences for each species. 
```bash
$ cd ATH_GMA
$ sh ./scripts/Section3.3.Step1.MakeIndex.sh
```
`Section3.3.Step1.MakeIndex.sh` contains two command lines to build index for Arabidopsis and soybean. Below is an example for Arabidopsis. 
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
When indexing is successfully finished, you can see newly generated files in `ATH_GMA/raw_data/[ATH|GMA]_STAR-2.5.2b_index` folder. 
```
# These are example files from indexing by STAR with Arabidopsis genome.
$ cd ATH_GMA/raw_data/ATH_STAR-2.5.2b_index
$ls -sh		# ls -sh shows file names and their human-readable file sizes. (zero size here does not mean empty.) 
total 2.8G
   0 chrLength.txt      9.8M exonGeTrInfo.tab     0 genomeParameters.txt  3.3M sjdbList.fromGTF.out.tab
   0 chrNameLength.txt  4.0M exonInfo.tab      1.2G SA                    3.3M sjdbList.out.tab
   0 chrName.txt        512K geneInfo.tab      1.5G SAindex               3.0M transcriptInfo.tab
   0 chrStart.txt       142M Genome            3.5M sjdbInfo.txt
```

2. Read mapping by STAR.
`Section3.3.Step2.Mapping.[ATH|GMA].sh` is for mapping all sequencing reads files in `ATH_GMA/raw_data/fastq/[ATH|GMA]/` folders. 

```bash
$ cd ATH_GMA
$ sh ./scripts/Section3.3.Step2.Mapping.ATH.sh
$ sh ./scripts/Section3.3.Step2.Mapping.GMA.sh
```
The following command below shows one example of such SRR ids (*_1 and *_2, paired-end reads files have two inputs).
```bash
$ STAR	--genomeDir $IDX \
		--readFilesIn $WORKDIR/raw_data/SRR2927328_1.fastq.gz   $WORKDIR/raw_data/SRR2927328_2.fastq.gz \
		--outFileNamePrefix $WORKDIR/processed_data/bam/SRR2927328/SRR2927328 \
		--outSAMtype BAM SortedByCoordinate
```

These are results from STAR mapping step. You can see mapping statistics from `*Log.final.out` files. To print out the file you can type `cat SRR2927328Log.final.out`.
```
$ cd ATH_GMA/processed_data/bam/SRR2927328/SRR2927328
$ ls -sh
total 2.1G
2.1G SRR2927328Aligned.sortedByCoord.out.bam  256K SRR2927328Log.out           3.8M SRR2927328SJ.out.tab
   0 SRR2927328Log.final.out                     0 SRR2927328Log.progress.out
$ cat SRR2927328Log.final.out
```

3. Read counting with featureCounts.
With mapping results from the previous step, FeatureCounts will calculate how many reads map to each gene region. 
`Section3.3.Step3.ReadCount.[ATH|GMA]` is designed to handle all mapping results (BAM files) under `ATH_GMA/processed_data/bam/[SRRid]` folrders. 

```bash
$ cd ATH_GMA
$ sh ./scripts/Section3.3.Step3.ReadCount.ATH.sh
$ sh ./scripts/Section3.3.Step3.ReadCount.GMA.sh
```
The following command below shows one example of such SRR ids (*_1 and *_2, paired-end reads files have two inputs).
```bash
$ WORKDIR=$(pwd)
$ GTF=$WORKDIR/raw_data/Araport11_GFF3_genes_transposons.201606.gtf
$ BAM=$WORKDIR/processed_data/bam
$ RC=$WORKDIR/processed_data/rc
$ featureCounts -t exon \
		-g gene_id \
		-p \
		-a $GTF \
		-o $RC/SRR2927328.readcount.txt \
		$BAM/SRR2927328/SRR2927328Aligned.sortedByCoord.out.bam
```

4. FPKM calculation using DESeq2 and edgeR.

We provide the unified R script, `ATH_GMA/scripts/Section3.3.Step4.FPKM.R` with table fils for replicate structure of the samples (PRJNA301162.csv for Arabidopsis and PRJNA197379.csv for soybean). 
The R script is composed of many blocks and will generate multiple out files. Hence, we paid attention to explain input requirements, contents of outputs, and functions of each step. For more details, please look at the script, `ATH_GMA/scripts/Section3.3.Step4.FPKM.R`. 

```bash
$ cd ATH_GMA
$ Rscript ./scripts/Section3.3.Step4.FPKM.R ./processed_data/fpkm/ATH	# for Arabidopsis data
$ Rscript ./scripts/Section3.3.Step4.FPKM.R ./processed_data/fpkm/GMA	# for soybean data
```

5. Expression data will be summarized and converted to gene co-expression networks. 

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

