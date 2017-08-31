#!/usr/bin/python

'''
This script will load BLAST results in the tabular format and
identify reciprocal best hit genes.

The BLAST output file is created using "-outfmt 6" paramter which
created tabular format output. Useful fields in this tabular format
include:
column 1: Query ID
column 2: Subject ID
column 11: e value.

For Arabidopsis to soybean blast, the original fasta files were formatted
such that each ATH gene will have an id start with ARATH and each soybean
gene will have an id start with GLYMA.

For example, the following line is from GLYMA fasta files blast to ARATH gene

GLYMA|Glyma.10G183400.1.p ARATH|AT3G18460.1 48.438  128 64  1 11  136 54  181 1.32e-37  129


Algorithm details:
Step 1. sort the data into two list, one for GLYMA blast to ARATH, one for ARATH
blast to GLYMA.

Step 2. for each GLYMA (ARATH) genes, find the one with lowest E value. create a
dictionary to save this result for each species.

Step 3. check each genes and find RBH hits.

Usage:

python ReciprocalBlastHit.py [inputfn] [species1id] [species2id] [outputfn]

Future changes:
make this more flexible by adding parameters to the script such that it can be
used in any pairwise comparisons.

'''

import sys
from operator import itemgetter, attrgetter
from pprint import pprint

# load two species into two lists
spec1 = []
spec2 = []
spec1n = sys.argv[2]
spec2n = sys.argv[3]
outfn = sys.argv[4]
for line in open(sys.argv[1],'r'):
    if line.find(spec1n) == 0:
        spec1.append(line.strip())
    if line.find(spec2n) == 0:
        spec2.append(line.strip())

def getgenename(input):
    # this is a parser that convert species specific isoform names
    # to gene names.
    #
    if input.find('ARATH')==0:
        return 'ARATH|'+input.split('|')[1].split('.')[0]
    elif input.find('GLYMA')==0:
        return 'GLYMA|'+'.'.join(input.split('|')[1].split('.')[0:2])
    else:
        return 'UnKnownSpecies'

def getbesthit(speclist):
    allhits = {} # this holds all BLAST hits for each gene
    for eg in speclist:
        tmp = eg.split('\t')
        queryname = getgenename(tmp[0])
        targetname = getgenename(tmp[1])
        eval = float(tmp[10])
        if not queryname in allhits:
            allhits[queryname]=[]
        allhits[queryname].append([targetname,eval])
    # sort allhits and get best ones.
    besthits= {}
    for eg in allhits:
        tmp = allhits[eg]
        besthita=sorted(tmp, key=itemgetter(1))
        besthitg=besthita[0]
        besthits[eg] = besthitg
        #pprint(besthita)
        #pprint(besthitg)
        #print eg,besthitg

    return besthits

spec1bh = getbesthit(spec1)
spec2bh = getbesthit(spec2)

#pprint(spec1bh)
#pprint(spec2bh)

outf = open(outfn,'w')
for eg in spec1bh:
    bhg = spec1bh[eg][0] # name of best hit
    bhe = spec1bh[eg][1] # e value for best hit
    if not bhg in spec2bh:
        continue # skip this if best hit is not even in other species list
    else:
        rbhg = spec2bh[bhg][0] # find best hit gene in species 2
        rbhe = spec2bh[bhg][1] # and associated e value.
        if rbhg == eg:
            outf.write(bhg+','+rbhg+','+`bhe`+','+`rbhe`+'\n')






