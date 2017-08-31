#!/bin/sh
### 2017-08-18 by Jiyoung Lee

### Scripts for section 2.1 Set up folder structure for data analysis

### Setting a working directory structure for tools ======================================
workdir=$(pwd)
mkdir raw_data processed_data scripts results software
mkdir processed_data/fpkm
mkdir software/bin

softwarepath=$workdir/software/bin
echo $softwarepath
cd software

### setting up the PATH for installed softwares
export PATH=$PATH:$softwarepath

