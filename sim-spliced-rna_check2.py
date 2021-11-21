#!/usr/bin/env python
# =========================================================================
# Copyright (C)William H. Majoros (bmajoros@alumni.duke.edu)
# =========================================================================
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
    generators,
    nested_scopes,
    with_statement,
)
from builtins import (
    bytes,
    dict,
    int,
    list,
    object,
    range,
    str,
    ascii,
    chr,
    hex,
    input,
    next,
    oct,
    open,
    pow,
    round,
    super,
    filter,
    map,
    zip,
)

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
from Bio.Seq import Seq

# ============================================
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
# ============================================

# GLOBAL VARIABLES:
matches = 0  # number of sites containing alt or ref allele in transcript
mismatches = 0  # number of sites having neither ref nor alt allele in transcript
rex = Rex()

# CLASSES:


class Variant:
    def __init__(self, fields):
        self.Chr = fields[0]
        self.genomicPos = int(fields[1]) - 1  # convert to 0-based
        self.ref = fields[3]
        self.alt = fields[4]
        genotype = fields[9]
        self.genotype = None
        if rex.find("(\d)\|(\d)", genotype):
            self.genotype = (int(rex[1]), int(rex[2]))
            if self.genotype[0] not in {0, 1} or self.genotype[1] not in {0, 1}:
                self.genotype = None

    def isOK(self):
        return self.genotype is not None

    def isHet(self):
        return self.genotype[0] != self.genotype[1]

    def isHomozygousAlt(self):
        return self.genotype[0] > 0 and self.genotype[1] > 0


class SamRecord:
    def __init__(self, fields):
        # self.fragLen=int(fields[8])
        # if(self.fragLen<0): self.fragLen=-self.fragLen
        self.readLen = len(fields[9])
        self.qual = fields[10]


# FUNCTIONS:


def makeAltCopy(gene, haplotype, variants):  # maternal=makeAltCopy(gene,0,variants)
    global matches
    global mismatches
    altGene = copy.deepcopy(gene)
    for transcript in gene.transcripts:
        array = list(transcript.sequence)
        for idx, variant in enumerate(variants):
            # transcript pos
            pos = transcript.mapToTranscript(variant.genomicPos)
            # pos=variant.genomicPos
            if pos < 0:
                continue
            if int(haplotype) == 0:  # count mismatches (sanity check)
                allele = array[pos]
                if gene.getStrand() == "-":
                    allele = Translation.reverseComplement(allele)
                if allele == variant.ref or allele == variant.alt:
                    matches += 1
                else:  #     print("        mismatch")
                    mismatches += 1
            if variant.genotype[int(haplotype)] > 0:
                array[pos] = variant.alt
        transcript.sequence = "".join(array)
    return altGene


#################################### randomly choose a quality score line


def nextSamRec(fh, filename):
    while True:
        line = fh.readline()
        if line is None:
            fh.close()
            fh = open(filename, "rt")
            line = fh.readline()
        fields = line.rstrip().split()
        if len(fields) < 11:
            continue
        if fields[0][0] == "@":
            continue
        if int(fields[8]) == 0:
            continue
        break
    return SamRecord(fields)


def tabixVCF(vcfFile, Chr, begin, end):
    records = []
    cmd = "tabix " + vcfFile + " " + Chr + ":" + str(begin) + "-" + str(end)
    lines = Pipe.run(cmd)
    for line in lines.split("\n"):
        if len(line) == 0 or line[0] == "#":
            continue
        fields = line.rstrip().split()
        if len(fields) != 10:
            raise Exception("Expecting 10 fields in VCF file")
        variant = Variant(fields)
        if not variant.isOK():
            continue
        if not variant.isHet() and not variant.isHomozygousAlt():
            continue
        records.append(variant)
    return records


def simRead(refTranscript, altTranscript, rec1, rec2, genome, path, fragLen):
    L = len(refTranscript.sequence)
    if L < fragLen or L < readLen:
        return None
    # L = 100
    # fragLen = 80
    lastStart = L - fragLen  # 20
    start1 = random.randrange(lastStart + 1)  # 10
    end1 = start1 + rec1.readLen  # 10+75 = 85
    end2 = start1 + fragLen  # 10+80 = 90
    start2 = end2 - rec2.readLen  # 90-75 = 15
    refSeq = refTranscript.sequence[start1:end1]
    refSeq_rev = Seq(refSeq).reverse_complement()
    # print("REF RWD: %s"%(refSeq))
    # print("REF REV: %s"%(refSeq))
    altSeq = altTranscript.sequence[start2:end2]
    altSeq_rev = Seq(altSeq).reverse_complement()
    # print("ALT RWD: %s"%(altSeq))
    # print("ALT REV: %s"%(altSeq))
    return (refSeq, refSeq_rev, altSeq, altSeq_rev, rec1.qual, rec2.qual)


def readFrom2bit(genome, path, Chr, start, end):
    cmd = (
        path
        + "/twoBitToFa -seq="
        + Chr
        + " -start="
        + str(start)
        + " -end="
        + str(end)
        + " "
        + genome
        + " stdout"
    )
    output = Pipe.run(cmd)
    fields = output.rstrip().split("\n")
    seq = ""
    for i in range(1, len(fields)):
        seq += fields[i]
    return seq


def loadTranscriptSeqs(gene, genome, path):
    Chr = gene.getSubstrate()
    strand = gene.getStrand()
    for transcript in gene.transcripts:
        transcript.sequence = ""
        transcript.exons = transcript.getRawExons()
        for exon in transcript.rawExons:
            exon.sequence = readFrom2bit(genome, path, Chr, exon.begin, exon.end)
            exon.sequence = exon.sequence.upper()
            if strand == "-":
                exon.sequence = Translation.reverseComplement(exon.sequence)
            if "N" in exon.sequence:
                raise Exception(
                    "N FOUND: Chr="
                    + Chr
                    + "begin="
                    + str(exon.begin)
                    + "end="
                    + str(exon.end)
                )
            transcript.sequence += exon.sequence


def pickTranscript(refGene, altGene):
    n = refGene.getNumTranscripts()
    i = random.randrange(n)
    return (refGene.getIthTranscript(i), altGene.getIthTranscript(i))


def getGeneVariants(gene, vcfFile):
    Chr = gene.getSubstrate()
    return tabixVCF(vcfFile, gene.getSubstrate(), gene.getBegin(), gene.getEnd())


def printRead(header, seq, qual, FH):
    print(header + "\n" + seq + "\n+\n" + qual, file=FH)


def loadFragLens(filename):
    array = []
    with open(filename, "rt") as IN:
        for line in IN:
            L = int(line)
            array.append(L)
    return array


# =========================================================================
# main()
# =========================================================================
if len(sys.argv) != 5:
    exit(
        ProgramName.get()
        + " <config-file> <per-base-read-depth> <out-read1.gz> <out-read2.gz> <prefix>\n"
    )
(configFile, DEPTH, outFile1, outFile2) = sys.argv[1:]
DEPTH = int(DEPTH)

# Load config file
configFile = ConfigFile(configFile)
twoBitDir = configFile.lookupOrDie("util-dir")
genome2bit = configFile.lookupOrDie("genome")
vcfFile = configFile.lookupOrDie("vcf")
samFile = configFile.lookupOrDie("aligned-rna")
gffFile = configFile.lookupOrDie("gff")
readLen = int(configFile.lookupOrDie("original-read-len"))
fragLenFile = configFile.lookupOrDie("fragment-lengths")
file_prefix = configFile.lookupOrDie("prefix")

# Load GFF and fragment lengths
gffReader = GffTranscriptReader()
print("reading GFF...", file=sys.stderr, flush=True)
genes = gffReader.loadGenes(gffFile)
fragLens = loadFragLens(fragLenFile)

# Create output files
OUT1 = gzip.open(outFile1, "wt")
OUT2 = gzip.open(outFile2, "wt")

# Simulate
print("simulating...", file=sys.stderr, flush=True)
nextReadID = 0
IN = open(samFile, "rt")
num_gene = 0
n_break = 0
# output list for debugging purpose
#
counter = 0
with open(file_prefix + "_geneID-numReads.txt", "w") as file_handler:
    for idx, gene in enumerate(genes):
        counter += 1
        num_gene += 1
        geneid = gene.getId()
        chrN = gene.getSubstrate()
        loadTranscriptSeqs(gene, genome2bit, twoBitDir)
        length = gene.longestTranscript().getLength()
        numReads = int(DEPTH * length / readLen)
        # all variants in gene region from VCF
        variants = getGeneVariants(gene, vcfFile)
        # if(len(variants)==0): continue
        maternal = makeAltCopy(gene, 0, variants)
        paternal = makeAltCopy(gene, 1, variants)
        print(
            "%s,gene: %s,#reads: %s,#variants: %d"
            % (chrN, geneid, numReads, len(variants))
        )
        file_handler.write("%s,%s,%s,%d\n" % (chrN, geneid, numReads, len(variants)))
        for i in range(numReads):
            (matTranscript, patTranscript) = pickTranscript(maternal, paternal)
            length = matTranscript.getLength()
            rec1 = nextSamRec(IN, samFile)
            rec2 = nextSamRec(IN, samFile)
            fragLen = fragLens[random.randrange(len(fragLens))]
            n_while_loop = 0
            if length < fragLen:
                fragLen = length
            sim = simRead(
                matTranscript, patTranscript, rec1, rec2, genome2bit, twoBitDir, fragLen
            )
            if sim is None:
                n_break += 1
                break  # gene is shorter than fragment length
            (refSeq, refSeq_rev, altSeq, altSeq_rev, qual1, qual2) = sim
            printRead("@READ" + str(nextReadID) + " SIM", refSeq, qual1, OUT1)
            printRead("@READ" + str(nextReadID) + " SIM", refSeq_rev, qual1, OUT2)
            nextReadID += 1
            printRead("@READ" + str(nextReadID) + " SIM", altSeq, qual2, OUT1)
            printRead("@READ" + str(nextReadID) + " SIM", altSeq_rev, qual2, OUT2)
            nextReadID += 1

print("saved debugging lists!")
print(">> total geneID processed: %d" % (num_gene))
print(">> total ReadID processed: %d" % (nextReadID))
print(">> total num breaks /gene is shorter than fragment length : %d" % (n_break))
# matches=0 # number of sites containing alt or ref allele in transcript
# mismatches=0 # number of sites having neither ref nor alt allele in transcript
print("matches: %s , mismatches : %d" % (matches, mismatches))
matchRate = float(matches) / float(matches + mismatches)
print("%s percent of exonic sites had ref or alt " % (matchRate * 100))

# MY SAM FILE:
# HWI-ST661:130:C037KACXX:6:2307:8186:172941      99      chr1    13170   109     75M     =       13278   183     CTGTTGGGGAGGCAGCTGTAACTCAAAGCCTTAGCCTCTGTTCCCACGAAGGCAGGGCCATCAGGCACCAAAGGG     @@@DDFFFHHFFFJIIIIHECCCGHIICEHJCEBFHEGIGIIDGHBDGGHDEGG?HFFFFFEEED>ADDBDBDDB     RG:Z:0  NM:i:3  XT:A:U  md:Z:75 XA:Z:chrY,-59359485,75M,1;chrX,-155256479,75M,1;chr9,+13283,75M,1;chr2,-114357772,75M,1;chr15,-102517926,75M,2;chr12,-92369,75M,2;

# SCARLETT'S SAM FILE:
# ERR188231.12567592      163     chr10   93307   255     75M     =       93442   210     GCAAAGTAACTGCTGTTCTTATCTTGAATGTTGAACATTTGTTCATCCACCTCCCTCATGGGCATGCGACCCCTG     CCCFFFDFHHHHHJJJJIJJJJJJJJJIJIIIJJJJIJJJJIJJIJJIJJGIIJIIJIJIJJIJJJJIIIJHHFF     MD:Z:75 PG:Z:MarkDuplicates     NH:i:1  HI:i:1  jI:B:i,-1       NM:i:0  jM:B:c,-1       nM:i:1  AS:i:146
