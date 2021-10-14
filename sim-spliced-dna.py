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
from GffTranscriptReader import GffTranscriptReader
from Pipe import Pipe
from ConfigFile import ConfigFile
from Translation import Translation

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

#def simRead(L,genome,path,Chr,start,end):
#    cmd=path+"/twoBitToFa -seq="+Chr+" -start="+str(start)+\
#        " -end="+str(end)+" "+genome+" stdout"
#    output=Pipe.run(cmd)
#    fields=output.rstrip().split("\n")
#    seq=fields[1] ### ?
#    return seq

def simRead(transcript,rec1,rec2,genome,path):
    L=len(transcript.sequence)
    if(L<rec1.fragLen or L<rec1.readLen or L<rec2.readLen): return None
    lastStart=L-rec1.fragLen
    #print("L=",L,"fraglen=",rec1.fragLen)
    start1=random.randrange(lastStart+1)
    end1=start1+rec1.readLen
    end2=start1+rec1.fragLen
    start2=end2-rec2.readLen
    seq1=readFrom2bit(genome,path,transcript.getSubstrate(),start1,end1)
    seq2=readFrom2bit(genome,path,transcript.getSubstrate(),start2,end2)
    return (seq1,rec1.qual,seq2,rec2.qual)
    
def readFrom2bit(genome,path,Chr,start,end):
    cmd=path+"/twoBitToFa -seq="+Chr+" -start="+str(start)+\
        " -end="+str(end)+" "+genome+" stdout"
    output=Pipe.run(cmd)
    fields=output.rstrip().split("\n")
    seq=fields[1]
    return seq

def loadTranscriptSeqs(gene,genome,path):
    Chr=gene.getSubstrate()
    for transcript in gene.transcripts:
        transcript.sequence=""
        transcript.getRawExons()
        for exon in transcript.rawExons:
            exon.sequence=readFrom2bit(genome,path,Chr,exon.begin,exon.end)
            transcript.sequence+=exon.sequence
    if(gene.getStrand()=="-"):
        transcript.sequence=Translation.reverseComplement(transcript.sequence)
        
def pickTranscript(gene):
    n=gene.getNumTranscripts()
    i=random.randrange(n)
    return gene.getIthTranscript(i)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <config-file> <#reads-per-gene>\n")
(configFile,NUM_READS)=sys.argv[1:]
NUM_READS=int(NUM_READS)

# Load config file
configFile=ConfigFile(configFile)
twoBitDir=configFile.lookupOrDie("util-dir")
genome2bit=configFile.lookupOrDie("genome")
vcfFile=configFile.lookupOrDie("vcf")
samFile=configFile.lookupOrDie("aligned-rna")
gffFile=configFile.lookupOrDie("gff")

# Load VCF
print("reading VCF...",file=sys.stderr,flush=True)
loadVCF(vcfFile)

# Load GFF
gffReader=GffTranscriptReader()
print("reading GFF...",file=sys.stderr,flush=True)
genes=gffReader.loadGenes(gffFile)

# Simulate
print("simulating...",file=sys.stderr,flush=True)
nextReadID=1
IN=open(samFile,"rt")
for gene in genes:
    loadTranscriptSeqs(gene,genome2bit,twoBitDir)
    for i in range(NUM_READS):
        transcript=pickTranscript(gene)
        rec1=nextSamRec(IN,samFile)
        rec2=nextSamRec(IN,samFile)
        Chr="chr1" ### DEBUGGING
        start1=100000; end1=100100; start2=200000; end2=200100; ### DEBUGGING
        #seq1=simRead(rec1.readLen,genome2bit,twoBitDir,Chr,start1,end1)
        #seq2=simRead(rec2.readLen,genome2bit,twoBitDir,Chr,start2,end2)
        sim=simRead(transcript,rec1,rec2,genome2bit,twoBitDir)
        if(sim is None): break # gene is shorter than fragment length
        (seq1,qual1,seq2,qual2)=sim
        print("@READ"+str(nextReadID),"SIM")
        print(seq1)
        print("+")
        print(qual1)
        print("@READ"+str(nextReadID),"SIM")
        print(seq2)
        print("+")
        print(qual2)
        nextReadID+=1


        
