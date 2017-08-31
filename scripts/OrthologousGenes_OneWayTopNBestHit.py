#!/usr/bin/python

#==> ARATH2GLYMA.BLAST_ARATH2GLYMA.txt <==
#AT1G02980.1	Glyma.17G025200.1	65.995	744	249	4	1	742	3	744	0.0	1046
#AT1G02980.1	Glyma.07G249100.1	66.129	744	248	4	1	742	3	744	0.0	1046
#==> ARATH2GLYMA.BLAST_GLYMA2ARATH.txt <==
#Glyma.10G183400.1	AT3G18460.1	48.438	128	64	1	11	136	54	181	1.32e-37	129
#Glyma.10G183400.1	AT3G18470.1	53.125	128	57	2	10	136	4	129	1.59e-35	122

#http://homepages.ulb.ac.be/~dgonze/TEACHING/stat_scores.pdf
# 1.	 qseqid	   query (e.g., gene) sequence id
# 2.	 sseqid	   subject (e.g., reference genome) sequence id
# 3.	 pident	   percentage of identical matches	<==4
# 4.	 length	   alignment length			<==3
# 11.	 evalue	   expect value				<==1
# 12.	 bitscore  bit score				<==2

### One way top N best hit. You can set your N number. 
N=5

import numpy as np
#inFile = "ARATH2GLYMA.BLAST_ARATH2GLYMA.txt"
inFile = "ARATH2GLYMA.BLAST_GLYMA2ARATH.txt"
outFile = inFile[0:-4]+"_Top5BestHit.txt"
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
print len(blastList), len(blastList[0])

import itertools
blastList.sort()
blastListUnq=list(blastList for blastList,_ in itertools.groupby(blastList))

print len(blastListUnq), len(blastListUnq[0])

blastListUnqArray=np.array(blastListUnq, dtype=object)

queries=list( set( np.array(blastListUnq)[:,0] ) )
print len(queries)

headerList=[["Query", "DB", "E-val", "Bits", "length", "PerctIdentity"]]
blastArrayTop5s = np.array(headerList)

print blastArrayTop5s
for q in queries:
	blastQarray = blastListUnqArray[np.where(blastListUnqArray[:,0] == q)].tolist()
	blastArraySort= np.array(sorted(sorted(sorted(sorted(sorted(blastQarray,key=lambda e:e[5],reverse=True),key=lambda e:e[4],reverse=True),key=lambda e:e[3],reverse=True),key=lambda e:e[2],reverse=False),key=lambda e:e[0],reverse=False))[0:N]
	blastArrayTop5s=np.concatenate((blastArrayTop5s, blastArraySort),axis=0)

import pandas as pd 
df = pd.DataFrame(blastArrayTop5s)
df.to_csv(outFile, header=None, index=False, sep="\t")

