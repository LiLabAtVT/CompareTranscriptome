#!/bin/sh

### Scripts for section 4.4 (Note 4.4) 
### USAGE: sh ./scripts/Section4.4_k-BestHits.sh
### Identify k-best-hits genes between two species using BLAST results

cd processed_data

### ex) k = 5, blast results with Arabidopsis query - soybean DB
python ../scripts/OrthologousGenes_OneWayTopNBestHit.py 5 ARATH2GLYMA.BLAST_ARATH2GLYMA.subset.txt

### ex) k = 5, blast results with soybean query - Arabidopsis DB
python ../scripts/OrthologousGenes_OneWayTopNBestHit.py 5 ARATH2GLYMA.BLAST_GLYMA2ARATH.subset.txt

### ex) k = 3, blast results with Arabidopsis query - soybean DB
#python ../scripts/OrthologousGenes_OneWayTopNBestHit.py 3 ARATH2GLYMA.BLAST_ARATH2GLYMA.subset.txt

### ex) k = 3, blast results with soybean query - Arabidopsis DB
#python ../scripts/OrthologousGenes_OneWayTopNBestHit.py 3 ARATH2GLYMA.BLAST_GLYMA2ARATH.subset.txt
