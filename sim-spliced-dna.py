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
from GffTranscriptReader import GffTranscriptReader
from Pipe import Pipe
from ConfigFile import ConfigFile
from Translation import Translation

# DEPENDENCIES:
#   module load htslib

class Variant:
    def __init__(self,fields):
        self.Chr=fields[0]
        self.genomicPos=int(fields[1])
        self.ref=fields[3]
        self.alt=fields[4]

class SamRecord:
    def __init__(self,fields):
        self.fragLen=int(fields[8])
        if(self.fragLen<0): self.fragLen=-self.fragLen
        self.readLen=len(fields[9])
        self.qual=fields[10]

def makeAltGene(gene,variants):
    altGene=copy.deepcopy(gene)
    changes=False
    for transcript in gene.transcripts:
        array=list(transcript.sequence)
        for variant in variants:
            pos=transcript.mapToTranscript(variant.genomicPos)
            if(pos<0): continue
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

def simRead(refTranscript,altTranscript,rec1,rec2,genome,path):
    L=len(refTranscript.sequence)
    if(L<rec1.fragLen or L<rec1.readLen or L<rec2.readLen): return None
    lastStart=L-rec1.fragLen
    start1=random.randrange(lastStart+1)
    end1=start1+rec1.readLen
    end2=start1+rec1.fragLen
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

def printRead(header,seq,qual):
    print(header+"\n"+seq+"\n+\n"+qual)

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
#print("reading VCF...",file=sys.stderr,flush=True)
#loadVCF(vcfFile)

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
    variants=getGeneVariants(gene,vcfFile)
    if(len(variants)==0): continue
    refGene=gene
    altGene=makeAltGene(gene,variants)
    if(altGene is None): continue # no variants in exons
    #print(len(variants),"variants in this gene")
    for i in range(NUM_READS):
        (refTranscript,altTranscript)=pickTranscript(refGene,altGene)
        rec1=nextSamRec(IN,samFile)
        rec2=nextSamRec(IN,samFile)
        sim=simRead(refTranscript,altTranscript,rec1,rec2,genome2bit,
                    twoBitDir)
        if(sim is None): break # gene is shorter than fragment length
        (refSeq1,altSeq1,qual1,refSeq2,altSeq2,qual2)=sim
        printRead("@READ"+str(nextReadID)+" SIM",refSeq1,qual1)
        printRead("@READ"+str(nextReadID)+" SIM",refSeq2,qual2)
        nextReadID+=1
        printRead("@READ"+str(nextReadID)+" SIM",altSeq1,qual1)
        printRead("@READ"+str(nextReadID)+" SIM",altSeq2,qual2)
        nextReadID+=1


        
