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
#   VCF file must be bzipped and indexed with tabix
#
# EXAMPLE CONFIG FILE:
#   util-dir = /hpc/home/bmajoros/twobit
#   genome = hg19.2bit
#   aligned-rna = HG00096.sam
#   fragment-lengths = fragmentlengths.txt
#   vcf = HG00096.vcf.gz
#   gff = gencode1000.gff
#   original-read-len = 75
#============================================

# GLOBAL VARIABLES:
matches=0 # number of sites containing alt or ref allele in transcript
mismatches=0 # number of sites having neither ref nor alt allele in transcript

# CLASSES:

class Variant:
    def __init__(self,fields):
        self.Chr=fields[0]
        self.genomicPos=int(fields[1])-1 # convert to 0-based
        self.ref=fields[3]
        self.alt=fields[4]

class SamRecord:
    def __init__(self,fields):
        #self.fragLen=int(fields[8])
        #if(self.fragLen<0): self.fragLen=-self.fragLen
        self.readLen=len(fields[9])
        self.qual=fields[10]

# FUNCTIONS:
        
def makeAltGene(gene,variants):
    global matches; global mismatches
    altGene=copy.deepcopy(gene)
    changes=False
    for transcript in gene.transcripts:
        array=list(transcript.sequence)
        for variant in variants:
            pos=transcript.mapToTranscript(variant.genomicPos)
            if(pos<0): continue
            c=array[pos]
            if(gene.getStrand()=="-"):
                c=Translation.reverseComplement(c)
            if(c==variant.ref or c==variant.alt): matches+=1
            else: mismatches+=1
            array[pos]=variant.alt
            changes=True
        transcript.sequence="".join(array)
    return altGene if changes else None
        
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

def tabixVCF(vcfFile,Chr,begin,end):
    records=[]
    cmd="tabix "+vcfFile+" "+Chr+":"+str(begin)+"-"+str(end)
    lines=Pipe.run(cmd)
    for line in lines.split("\n"):
        if(len(line)==0 or line[0]=="#"): continue
        fields=line.rstrip().split()
        variant=Variant(fields)
        records.append(variant)
    return records

def simRead(refTranscript,altTranscript,rec1,rec2,genome,path,fragLen):
    L=len(refTranscript.sequence)
    if(L<fragLen or L<readLen or L<readLen): return None
    lastStart=L-fragLen
    start1=random.randrange(lastStart+1)
    end1=start1+rec1.readLen
    end2=start1+fragLen
    start2=end2-rec2.readLen
    refSeq1=refTranscript.sequence[start1:end1]
    refSeq2=refTranscript.sequence[start2:end2]
    altSeq1=altTranscript.sequence[start1:end1]
    altSeq2=altTranscript.sequence[start2:end2]
    return (refSeq1,altSeq1,rec1.qual,refSeq2,altSeq2,rec2.qual)
    
def readFrom2bit(genome,path,Chr,start,end):
    cmd=path+"/twoBitToFa -seq="+Chr+" -start="+str(start)+\
        " -end="+str(end)+" "+genome+" stdout"
    output=Pipe.run(cmd)
    fields=output.rstrip().split("\n")
    seq=""
    for i in range(1,len(fields)): seq+=fields[i]
    return seq

def loadTranscriptSeqs(gene,genome,path):
    Chr=gene.getSubstrate()
    strand=gene.getStrand()
    for transcript in gene.transcripts:
        transcript.sequence=""
        transcript.exons=transcript.getRawExons()
        for exon in transcript.rawExons:
            exon.sequence=readFrom2bit(genome,path,Chr,exon.begin,exon.end)
            exon.sequence=exon.sequence.upper()
            if(strand=="-"):
                exon.sequence=Translation.reverseComplement(exon.sequence)
            if("N" in exon.sequence):
                raise Exception("N FOUND: Chr="+Chr+"begin="+str(exon.begin)+\
                                "end="+str(exon.end))
            transcript.sequence+=exon.sequence
        
def pickTranscript(refGene,altGene):
    n=refGene.getNumTranscripts()
    i=random.randrange(n)
    return (refGene.getIthTranscript(i),altGene.getIthTranscript(i))

def getGeneVariants(gene,vcfFile):
    Chr=gene.getSubstrate()
    return tabixVCF(vcfFile,gene.getSubstrate(),gene.getBegin(),
                    gene.getEnd())

def printRead(header,seq,qual,FH):
    print(header+"\n"+seq+"\n+\n"+qual,file=FH)

def loadFragLens(filename):
    array=[]
    with open(filename,"rt") as IN:
        for line in IN:
            L=int(line)
            array.append(L)
    return array
    
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <config-file> <per-base-read-depth> <out-read1.gz> <out-read2.gz>\n")
(configFile,DEPTH,outFile1,outFile2)=sys.argv[1:]
DEPTH=int(DEPTH)

# Load config file
configFile=ConfigFile(configFile)
twoBitDir=configFile.lookupOrDie("util-dir")
genome2bit=configFile.lookupOrDie("genome")
vcfFile=configFile.lookupOrDie("vcf")
samFile=configFile.lookupOrDie("aligned-rna")
gffFile=configFile.lookupOrDie("gff")
readLen=int(configFile.lookupOrDie("original-read-len"))
fragLenFile=configFile.lookupOrDie("fragment-lengths")

# Load GFF and fragment lengths
gffReader=GffTranscriptReader()
print("reading GFF...",file=sys.stderr,flush=True)
genes=gffReader.loadGenes(gffFile)
fragLens=loadFragLens(fragLenFile)

# Create output files
OUT1=gzip.open(outFile1,"wt")
OUT2=gzip.open(outFile2,"wt")

# Simulate
print("simulating...",file=sys.stderr,flush=True)
nextReadID=1
IN=open(samFile,"rt")
for gene in genes:
    loadTranscriptSeqs(gene,genome2bit,twoBitDir)
    length=gene.longestTranscript().getLength()
    numReads=int(DEPTH*length/readLen)
    variants=getGeneVariants(gene,vcfFile)
    if(len(variants)==0): continue
    refGene=gene
    altGene=makeAltGene(gene,variants)
    if(altGene is None): continue # no variants in exons
    #print(len(variants),"variants in this gene")
    for i in range(numReads):
        (refTranscript,altTranscript)=pickTranscript(refGene,altGene)
        rec1=nextSamRec(IN,samFile)
        rec2=nextSamRec(IN,samFile)
        fragLen=fragLens[random.randrange(len(fragLens))]
        sim=simRead(refTranscript,altTranscript,rec1,rec2,genome2bit,
                    twoBitDir,fragLen)
        if(sim is None): break # gene is shorter than fragment length
        (refSeq1,altSeq1,qual1,refSeq2,altSeq2,qual2)=sim
        printRead("@READ"+str(nextReadID)+" SIM",refSeq1,qual1,OUT1)
        printRead("@READ"+str(nextReadID)+" SIM",refSeq2,qual2,OUT2)
        nextReadID+=1
        printRead("@READ"+str(nextReadID)+" SIM",altSeq1,qual1,OUT1)
        printRead("@READ"+str(nextReadID)+" SIM",altSeq2,qual2,OUT2)
        nextReadID+=1
matchRate=float(matches)/float(matches+mismatches)
print(matchRate*100,"% exonic sites had ref or alt: ",matches,"/",
      matches+mismatches,sep="",file=sys.stderr)


        
