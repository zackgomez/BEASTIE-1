#!/usr/bin/env python
#=========================================================================
# Copyright (C)William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import ProgramName
import random
import copy
import gzip
from GffTranscriptReader import GffTranscriptReader
from Pipe import Pipe
from ConfigFile import ConfigFile
from Translation import Translation

#============================================
# DEPENDENCIES:
#   module load htslib
#   python version 3.x
#   SAM file must be bzipped and indexed with tabix
#============================================

class SamRecord:
    def __init__(self,fields):
        self.pos1=int(fields[3])
        self.tlen=int(fields[8])
        #self.end=self.pos1+int(fields[8])

def nextSamRec(fh,filename):
    while(True):
        line=fh.readline()
        if(line is None):
            fh.close()
            fh=open(filename,"rt")
            line=fh.readline()
        fields=line.rstrip().split()
        if(len(fields)<11): continue
        if(fields[0][0]=="@"): continue
        if(int(fields[8])==0): continue
        break
    return SamRecord(fields)

def tabixSAM(filename,Chr,begin,end):
    records=[]
    cmd="tabix "+filename+" "+Chr+":"+str(begin)+"-"+str(end)
    lines=Pipe.run(cmd)
    for line in lines.split("\n"):
        if(len(line)==0 or line[0]=="#"): continue
        fields=line.rstrip().split()
        rec=SamRecord(fields)
        records.append(rec)
    return records

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <indexed-sam.gz> <gff>\n")
(samFile,gffFile)=sys.argv[1:]

# Load GFF
gffReader=GffTranscriptReader()
genes=gffReader.loadGenes(gffFile)

# Process SAM records gene-by-gene
for gene in genes:
    transcript=gene.longestTranscript()
    transcript.exons=transcript.getRawExons()
    transcript.recomputeBoundaries()
    length=transcript.getLength()
    recs=tabixSAM(samFile,gene.getSubstrate(),transcript.getBegin(),
                  transcript.getEnd())
    for rec in recs:
        if(rec.tlen<0): continue
        begin=transcript.mapToTranscript(rec.pos1)
        if(begin<0): continue
        end=transcript.mapToTranscript(rec.pos1+rec.tlen)
        if(end<0): continue
        fragLen=abs(end-begin)
        if(fragLen<100): continue
        print(fragLen,length,sep="\t")
    
