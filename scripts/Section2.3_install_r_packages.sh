#!/bin/sh
### Scripts for section 2.3 install R packages 
### To install R packages for this analysis using a bash command

### Check R version  =============================================================
echo -e "\n\n`date +"%b%d|%T"`\t\c";
echo "======== Checking R, Rscript versions ========"
R --version
Rscript --version

echo "======== Installing R packages ==============="
Rscript ./scripts/Section2.3_install_r_packages.R 