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
from GffTranscriptReader import GffTranscriptReader
from Pipe import Pipe

class Variant:
    def __init__(self,fields):
        self.Chr=fields[0]
        self.pos=int(fields[1])
        self.ref=fields[3]
        self.alt=fields[4]

class SamRecord:
    def __init__(self,fields):
        self.fragLen=int(fields[8])
        if(self.fragLen<0): self.fragLen=-self.fragLen
        self.readLen=len(fields[9])
        self.qual=fields[10]

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

def loadVCF(filename):
    records=[]
    with open(filename,"rt") as IN:
        for line in IN:
            if(len(line)==0 or line[0]=="#"): continue
            fields=line.rstrip().split()
            variant=Variant(fields)
            records.append(variant)
    return records

def simRead(L,genome,path,Chr,start,end):
    cmd=path+"/twoBitToFa -seq="+Chr+" -start="+str(start)+\
        " -end="+str(end)+" "+genome+" stdout"
    output=Pipe.run(cmd)
    fields=output.rstrip().split("\n")
    seq=fields[1] ### ?
    return seq
    
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=7):
    exit(ProgramName.get()+" <#reads> <rna.sam> <vcf> <gff> <genome.2bit> <path-to-twoBitToFa>\n")
(NUM_READS,samFile,vcfFile,gffFile,genome2bit,twoBitDir)=sys.argv[1:]
NUM_READS=int(NUM_READS)

# Load VCF
#print("reading VCF...",flush=True)
loadVCF(vcfFile)

# Load GFF
gffReader=GffTranscriptReader()
#print("reading GFF...",flush=True)
genes=gffReader.loadGenes(gffFile)

# Simulate
print("simulating...",flush=True)
nextReadID=1
IN=open(samFile,"rt")
for i in range(NUM_READS):
    rec1=nextSamRec(IN,samFile)
    rec2=nextSamRec(IN,samFile)
    Chr="chr1" ###
    start1=100000; end1=100100; start2=200000; end2=200100; ###
    seq1=simRead(rec1.readLen,genome2bit,twoBitDir,Chr,start1,end1)
    seq2=simRead(rec2.readLen,genome2bit,twoBitDir,Chr,start2,end2)
    print("@"+str(nextReadID),"SIM")
    #print("ACGT....")
    print(seq1)
    print("+")
    print(rec1.qual)
    print("@"+str(nextReadID),"SIM")
    #print("ACGT....")
    print(seq2)
    print("+")
    print(rec2.qual)
    nextReadID+=1
