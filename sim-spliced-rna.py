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
from Rex import Rex

#============================================
# DEPENDENCIES:
#   module load htslib
#   python version 3.x
#   VCF file must be bzipped and indexed with tabix.
#   VCF must contain only one sample (individual).
#   VCF must include all sites, including homozygous and heterozygous sites.
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
rex=Rex()

# CLASSES:

class Variant:
    def __init__(self,fields):
        self.Chr=fields[0]
        self.genomicPos=int(fields[1])-1 # convert to 0-based
        self.ref=fields[3]
        self.alt=fields[4]
        genotype=fields[9]
        self.genotype=None
        if(rex.find("(\d)\|(\d)",genotype)):
            self.genotype=(int(rex[1]),int(rex[2]))
            if(self.genotype[0] not in {0,1} or
               self.genotype[1] not in {0,1}):
                self.genotype=None
    def isOK(self):
        return self.genotype is not None
    def isHet(self):
        return self.genotype[0]!=self.genotype[1]
    def isHomozygousAlt(self):
        return self.genotype[0]>0 and self.genotype[1]>0

class SamRecord:
    def __init__(self,fields):
        #self.fragLen=int(fields[8])
        #if(self.fragLen<0): self.fragLen=-self.fragLen
        self.readLen=len(fields[9])
        self.qual=fields[10]

# FUNCTIONS:
        
def makeAltCopy(gene,haplotype,variants):
    #print("make")
    # maternal=makeAltCopy(gene,0,variants)
    global matches; global mismatches
    altGene=copy.deepcopy(gene)
    for transcript in gene.transcripts:
        array=list(transcript.sequence)
        for idx, variant in enumerate(variants):
            pos=transcript.mapToTranscript(variant.genomicPos)
            #pos=variant.genomicPos
            if(pos<0): continue
            print("    %sth variant, pos is %s"%(idx+1,pos))
            if haplotype == 0:
                print("    Assigned haplotype is maternal: %s"%(haplotype))
            else:
                print("    Assigned haplotype is paternal: %s"%(haplotype))
            print("    variant genotype: %s"%(str(variant.genotype)))
            print("    array size is %s, index is %s, pos is %s"%(len(array),idx,pos))
            #print(array[idx])
            #print(array[pos])
            #break
            if(int(haplotype)==0): # count mismatches (sanity check)
                allele=array[pos]
                print("    allele : %s ; variant.ref : %s ; variant.alt : %s"%(allele,variant.ref,variant.alt))
                if(gene.getStrand()=="-"): allele=Translation.reverseComplement(allele)
                if(allele==variant.ref or allele==variant.alt): 
                    matches+=1
                    if(allele==variant.ref):
                        print("        match to REF")
                    else:
                        print("        match to ALT")
                else: 
                    print("        mismatch")
                    mismatches+=1
            if(variant.genotype[int(haplotype)]>0):
                array[pos]=variant.alt
        transcript.sequence="".join(array)
    return altGene

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
        if(len(fields)!=10): raise Exception("Expecting 10 fields in VCF file")
        variant=Variant(fields)
        if(not variant.isOK()): continue
        if(not variant.isHet() and not variant.isHomozygousAlt()): continue
        records.append(variant)
    return records

def simRead(refTranscript,altTranscript,rec1,rec2,genome,path,fragLen):
    L=len(refTranscript.sequence)
    if(L<fragLen or L<readLen): return None
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
    exit(ProgramName.get()+" <config-file> <per-base-read-depth> <out-read1.gz> <out-read2.gz> <prefix>\n")
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
file_prefix=configFile.lookupOrDie("prefix")

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
nextReadID=0
IN=open(samFile,"rt")
num_gene=0

# output list for debugging purpose
#
with open(file_prefix+"_geneID-numReads.txt", 'w') as file_handler:
    for idx,gene in enumerate(genes):
        num_gene+=1
        geneid = gene.getId()
        loadTranscriptSeqs(gene,genome2bit,twoBitDir)
        length=gene.longestTranscript().getLength()
        numReads=int(DEPTH*length/readLen)
        #
        file_handler.write("%s,%s\n"%(geneid,numReads))
        #
        variants=getGeneVariants(gene,vcfFile)
        print("\t%sth gene: %s, %s reads, %s variants in this gene"%(idx+1,geneid,numReads,len(variants)))
        #if(len(variants)==0): continue
        maternal=makeAltCopy(gene,0,variants)
        #print(maternal)
        paternal=makeAltCopy(gene,1,variants)
        #print(paternal)
        #if(paternal is None): continue # no variants in exons
        for i in range(numReads):
            #print(">>>> num reads : %s"%(i))
            (matTranscript,patTranscript)=pickTranscript(maternal,paternal)
            rec1=nextSamRec(IN,samFile)
            rec2=nextSamRec(IN,samFile)
            fragLen=fragLens[random.randrange(len(fragLens))]
            while length < fragLen:
                fragLen=fragLens[random.randrange(len(fragLens))]
                #print("%s >= %s"%(length,fragLen))
            sim=simRead(matTranscript,patTranscript,rec1,rec2,genome2bit,
                            twoBitDir,fragLen)
            if(sim is None): 
                #print(">>>>>>> break")
                break # gene is shorter than fragment length
            (refSeq1,altSeq1,qual1,refSeq2,altSeq2,qual2)=sim
            printRead("@READ"+str(nextReadID)+" SIM",refSeq1,qual1,OUT1)
            printRead("@READ"+str(nextReadID)+" SIM",refSeq2,qual2,OUT2)
            nextReadID+=1
            printRead("@READ"+str(nextReadID)+" SIM",altSeq1,qual1,OUT1)
            printRead("@READ"+str(nextReadID)+" SIM",altSeq2,qual2,OUT2)
            nextReadID+=1

geneId_file=file_prefix+"_geneID_debug.txt"
numReads_file=file_prefix+"_numReads_debug.txt"
print("saved debugging lists!")
print(">> total geneID processed: %d"%(num_gene))
print(">> total ReadID processed: %d"%(nextReadID))
#matches=0 # number of sites containing alt or ref allele in transcript
#mismatches=0 # number of sites having neither ref nor alt allele in transcript
print("matches: %s , mismatches : %d"%(matches,mismatches))
matchRate=float(matches)/float(matches+mismatches)
print("%s percent of exonic sites had ref or alt "%(matchRate*100))

# MY SAM FILE:
# HWI-ST661:130:C037KACXX:6:2307:8186:172941      99      chr1    13170   109     75M     =       13278   183     CTGTTGGGGAGGCAGCTGTAACTCAAAGCCTTAGCCTCTGTTCCCACGAAGGCAGGGCCATCAGGCACCAAAGGG     @@@DDFFFHHFFFJIIIIHECCCGHIICEHJCEBFHEGIGIIDGHBDGGHDEGG?HFFFFFEEED>ADDBDBDDB     RG:Z:0  NM:i:3  XT:A:U  md:Z:75 XA:Z:chrY,-59359485,75M,1;chrX,-155256479,75M,1;chr9,+13283,75M,1;chr2,-114357772,75M,1;chr15,-102517926,75M,2;chr12,-92369,75M,2;

# SCARLETT'S SAM FILE:
# ERR188231.12567592      163     chr10   93307   255     75M     =       93442   210     GCAAAGTAACTGCTGTTCTTATCTTGAATGTTGAACATTTGTTCATCCACCTCCCTCATGGGCATGCGACCCCTG     CCCFFFDFHHHHHJJJJIJJJJJJJJJIJIIIJJJJIJJJJIJJIJJIJJGIIJIIJIJIJJIJJJJIIIJHHFF     MD:Z:75 PG:Z:MarkDuplicates     NH:i:1  HI:i:1  jI:B:i,-1       NM:i:0  jM:B:c,-1       nM:i:1  AS:i:146


        
