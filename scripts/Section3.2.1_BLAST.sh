#!/bin/sh
### Scripts for section 3.2.1 

### Prepare data for BLAST analysis  ===============================================================
### Before execute this step, make sure the raw data of protein sequences were downloaded in raw_data folder.
### Araport11.pep.fasta is the fasta file for Arabidopsis proteins
### GLYMA2.pep.fasta is the fasta file for soybean
### The actual names for these files may differ in different version of the genome annotation.
cat ./raw_data/Araport11.pep.fasta ./raw_data/GLYMA2.pep.fasta > ./processed_data/ATHGMA.pep.fasta


### Build BLAST database  ===============================================================
cd processed_data
makeblastdb -in ATHGMA.pep.fasta \
            -out ATHGMA.pep.blastdb \
            -dbtype prot \
            -logfile makeblastdb.log

### Perform BLAST analysis  ===============================================================
# this step is performed in "processed_data" folder
blastp -evalue 0.00001 \
  -outfmt 6 -db ATHGMAX.pep.blastdb \
  -query ATHGMA.pep.fasta > ATHGMA.pep.blastout 




