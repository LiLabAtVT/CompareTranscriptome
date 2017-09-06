#!/usr/bin/python
### 2017-08-07 by Jiyoung Lee 
### A script for Note 4.4 and section 3.2.2. 
### USAGE: python XXX.py [INT] [INPUT]
### Identify k-best-hits genes between two species using BLAST results

#==> ex) ARATH2GLYMA.BLAST_ARATH2GLYMA.subset.txt <==
#AT1G02980.1	Glyma.17G025200.1	65.995	744	249	4	1	742	3	744	0.0	1046
#AT1G02980.1	Glyma.07G249100.1	66.129	744	248	4	1	742	3	744	0.0	1046

#==> ex) ARATH2GLYMA.BLAST_GLYMA2ARATH.subset.txt <==
#Glyma.10G183400.1	AT3G18460.1	48.438	128	64	1	11	136	54	181	1.32e-37	129
#Glyma.10G183400.1	AT3G18470.1	53.125	128	57	2	10	136	4	129	1.59e-35	122

#==> How we decided best hits among blast results with 4 criteria.
#http://homepages.ulb.ac.be/~dgonze/TEACHING/stat_scores.pdf
# 1.	 qseqid	   query (e.g., gene) sequence id
# 2.	 sseqid	   subject (e.g., reference genome) sequence id
# 3.	 pident	   percentage of identical matches	<==4th
# 4.	 length	   alignment length			<==3rd
# 11.	 evalue	   expect value				<==1st
# 12.	 bitscore  bit score				<==2nd

import sys

TOPN = int(sys.argv[1])
print "=== Top "+str(TOPN)+" One Way Best Hit results ==="

import numpy as np
inFile = sys.argv[2]
outFile = inFile[0:-4]+"_Top-"+str(TOPN)+"_OneWayBestHit.txt"
print "- Output file    : ", 
print outFile

with open(inFile, "r") as f:
	blastList=[]

	for line in f:
		bList=[]
		line=line.strip("\n\r").split('\t')
		qr= '.'.join(line[0].split(".")[0:-1])
		db= '.'.join(line[1].split(".")[0:-1])
		evl=float(line[10])
		bit=float(line[11])
		lth=int(line[3])
		idn=float(line[2])
		bList=[qr, db, evl, bit, lth, idn]
		blastList.append(bList)

print "- blastList All  : ", 		
print len(blastList)

import itertools
blastList.sort()
blastListUnq=list(blastList for blastList,_ in itertools.groupby(blastList))

print "- blastList Uniq : ", 
print len(blastListUnq)

blastListUnqArray=np.array(blastListUnq, dtype=object)

print "- length of Query: ", 
queries=list( set( np.array(blastListUnq)[:,0] ) )
print len(queries)

headerList=[["Query", "DB", "E-val", "Bits", "length", "PerctIdentity"]]
blastArrayTop5s = np.array(headerList)

#print blastArrayTop5s
for q in queries:
	blastQarray = blastListUnqArray[np.where(blastListUnqArray[:,0] == q)].tolist()
	blastArraySort= np.array(sorted(sorted(sorted(sorted(sorted(blastQarray,key=lambda e:e[5],reverse=True),key=lambda e:e[4],reverse=True),key=lambda e:e[3],reverse=True),key=lambda e:e[2],reverse=False),key=lambda e:e[0],reverse=False))[0:TOPN]
	blastArrayTop5s=np.concatenate((blastArrayTop5s, blastArraySort),axis=0)

print "- Num of results : ", 
print len(blastArrayTop5s)-1   #exclude the header line for counting lines.

import pandas as pd 
df = pd.DataFrame(blastArrayTop5s)
df.to_csv(outFile, header=None, index=False, sep="\t")

