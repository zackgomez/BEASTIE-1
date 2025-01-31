#!/usr/bin/env python
# =========================================================================
# Copyright (C)William H. Majoros (bmajoros@alumni.duke.edu)
# =========================================================================
# This version has latest modification by Scarlett
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
from collections import defaultdict
from GffTranscriptReader import GffTranscriptReader
from Pipe import Pipe
from ConfigFile import ConfigFile
from Translation import Translation
from Rex import Rex
from Bio.Seq import Seq
from pathlib import Path

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
    transcriptIdToBiSNPcount = {}
    transcriptIdToBiSNPpos = defaultdict(set)
    if haplotype == 0:
        print("    ref/paternal")
    else:
        print("    alt/maternal")
    transcript_num = 0
    # loop through each transcript
    if_bi = False
    for transcript in altGene.transcripts:
        array = list(transcript.sequence)
        transcript_num += 1
        num = 0
        # loop through each bi-allelic SNP
        for variant in variants:
            if_rev = False
            trans_pos = transcript.mapToTranscript(variant.genomicPos)
            allele_in_ref = array[trans_pos]
            if len(variant.ref) == 1 and len(variant.alt) == 1 and trans_pos > 0:
                # if_bi = True
                print(
                    " genomic pos is %d, trans pos is %d"
                    % (variant.genomicPos, trans_pos)
                )
                if variant.genotype[0] != variant.genotype[1]:
                    if_bi = True
                    transcriptIdToBiSNPpos[transcript.getID()].add(trans_pos)
                    num += 1
                ########################## could be commented out
                # if int(haplotype) == 0:  # count mismatches (sanity check)
                #     allele = array[trans_pos]
                #     if gene.getStrand() == "-":
                #         allele = Translation.reverseComplement(allele)
                #     if allele == variant.ref or allele == variant.alt:
                #         matches += 1
                #         mismatches += 1
                ########################## could be commented out
                # 0|0 case and (1|0 or 0|1 case)
                if (
                    variant.genotype[0]
                    == variant.genotype[1]  # (maternal or paternal) with 0|0
                    and variant.genotype[0] == 0
                ) or (
                    variant.genotype[0] != variant.genotype[1]
                    and int(haplotype) == 0  # paternal copy with het SNP
                ):
                    allele = variant.ref
                # 1|1 case and (1|0 or 0|1 case)
                elif (
                    variant.genotype[0]
                    == variant.genotype[1]  # (maternal or paternal) with 1|1
                    and variant.genotype[0] == 1
                    or (
                        variant.genotype[0]
                        != variant.genotype[1]  # maternal copy with het SNP
                        and int(haplotype) == 1
                    )
                ):
                    allele = variant.alt
                # if it is reverse strand, then take care of it:
                if gene.getStrand() == "-":
                    if_rev = True
                    array[trans_pos] = Translation.reverseComplement(allele)
                else:
                    array[trans_pos] = allele
                print(
                    "        %dth transcript, %dth bi-allelic SNP"
                    % (transcript_num, num)
                )
                print(
                    "            > if reverse strand: %s, variant trans pos: %s, haplotype: %s, in ref genome: %s,ref: %s, alt: %s, write_in_sequence: %s"
                    % (
                        str(if_rev),
                        trans_pos,
                        variant.genotype,
                        allele_in_ref,
                        variant.ref,
                        variant.alt,
                        array[trans_pos],
                    )
                )
        transcript.sequence = "".join(array)
        transcriptIdToBiSNPcount[transcript.getID()] = num
        # print(transcriptIdToBiSNPpos)
        print(
            "        >>>>>> %dth transcript, in total %d bi-allelic SNP"
            % (transcript_num, num)
        )
        print(" ")
    if if_bi == True:
        return (
            altGene,
            transcriptIdToBiSNPcount,
            transcriptIdToBiSNPpos,
        )  # , transcriptIdToBiSNPpos
    else:
        return None, None, None


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
    # print(cmd)
    lines = Pipe.run(cmd)
    # print(lines)
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


def simRead(
    n_while,
    geneid,
    nextReadID,
    refTranscript,
    altTranscript,
    rec1,
    rec2,
    fragLen,
    trans_pos,
):
    L = len(refTranscript.sequence)
    if L < fragLen or L < readLen:
        return None
    # L = 100
    # fragLen = 80
    lower_bound = trans_pos - rec1.readLen + 1
    upper_bound = trans_pos + rec1.readLen
    if lower_bound < 1:
        lower_bound = 1
    if upper_bound > L:
        upper_bound = L
    if upper_bound - rec1.readLen - 20 < lower_bound:
        return (None, None, None, None, None, None)

    ref_start = random.randrange(lower_bound, upper_bound - rec1.readLen - 20)  # 10
    ref_end = ref_start + rec1.readLen  # 10+75 = 85

    ref_rev_end = ref_start + fragLen  # 10+80 = 90
    ref_rev_start = ref_rev_end - rec2.readLen  # 90-75 = 15

    alt_start = random.randint(lower_bound + 20, upper_bound - rec1.readLen)  # 10
    alt_end = alt_start + rec1.readLen  # 10+75 = 85

    alt_rev_end = alt_start + fragLen  # 10+80 = 90
    alt_rev_start = alt_rev_end - rec2.readLen  # 90-75 = 15

    if_overlap = False
    if (
        (alt_start > ref_start and alt_end > ref_end)
        and (ref_start < trans_pos and ref_end > trans_pos)
        and (alt_start < trans_pos and alt_end > trans_pos)
    ):
        if_overlap = True

    refSeq = refTranscript.sequence[ref_start:ref_end]
    refSeq_rev = Seq(
        refTranscript.sequence[ref_rev_start:ref_rev_end]
    ).reverse_complement()
    altSeq = altTranscript.sequence[alt_start:alt_end]
    altSeq_rev = Seq(
        altTranscript.sequence[alt_rev_start:alt_rev_end]
    ).reverse_complement()

    if str(refSeq) == str(altSeq):
        return (None, None, None, None, None, None)
    else:
        print(
            "while_loop:%d ,geneid: %s - readID: %s, pos %d, FWD seq length %d,REV seq length %d, REF FWD start-end : %d-%d, REF REV start-end : %d-%d,ALT FWD start-end : %d-%d, ALT REV start-end : %d-%d, if overlap? %s"
            % (
                n_while,
                geneid,
                nextReadID,
                trans_pos,
                rec1.readLen,
                rec2.readLen,
                ref_start,
                ref_end,
                ref_rev_start,
                ref_rev_end,
                alt_start,
                alt_end,
                alt_rev_start,
                alt_rev_end,
                if_overlap,
            )
        )
        print("REF FWD: %s" % (refSeq))
        print("REF REV: %s" % (refSeq_rev))
        print("ALT FWD: %s" % (altSeq))
        print("ALT REV: %s" % (altSeq_rev))
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


def pickTranscript(refGene, altGene, transcriptIdToBiSNPcount, transcriptIdToBiSNPpos):
    # print("pick transcript")
    # print(transcriptIdToBiSNPcount)
    n1 = refGene.getNumTranscripts()
    i1 = random.randrange(n1)

    num = transcriptIdToBiSNPcount[refGene.getIthTranscript(i1).getID()]
    while num == 0:
        n1 = refGene.getNumTranscripts()
        i1 = random.randrange(n1)
        num = transcriptIdToBiSNPcount[refGene.getIthTranscript(i1).getID()]
    print(
        "DEBUG: pick transcript step: chosen %d-th transcript %s, number of biSNP count %d"
        % (i1 + 1, refGene.getIthTranscript(i1).getID(), num)
    )
    print(transcriptIdToBiSNPpos)
    set1 = transcriptIdToBiSNPpos[refGene.getIthTranscript(i1).getID()]
    list1 = [int(x) for x in set1]
    # make sure to pick a transcript with biSNP
    list1_sorted = sorted(list1)
    print(list1_sorted)
    transpos_chosen = random.choice(list1)
    # print(list1_sorted)
    # transpos_chosen = list1_sorted[0]
    print("random choose pos %d" % (transpos_chosen))
    return (
        refGene.getIthTranscript(i1),
        altGene.getIthTranscript(i1),
        i1,
        refGene.getIthTranscript(i1).getID(),
        num,
        transpos_chosen,
    )


def getGeneVariants(gene, vcfFile):
    # Chr = gene.getSubstrate()
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
if len(sys.argv) != 6:
    exit(
        ProgramName.get()
        + " <config-file> <per-base-read-depth> <prefix> <out-read1.gz> <out-read2.gz>\n"
    )
(configFile, DEPTH, prefix, outFile1, outFile2) = sys.argv[1:]
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
out_path = configFile.lookupOrDie("out_path")

output_path = out_path + prefix
outpath_refalt = output_path + "/refalt"
Path(output_path).mkdir(parents=True, exist_ok=True)
Path(outpath_refalt).mkdir(parents=True, exist_ok=True)


# Load GFF and fragment lengths
gffReader = GffTranscriptReader()
print("reading GFF...", file=sys.stderr, flush=True)
genes = gffReader.loadGenes(gffFile)
genes.sort(key=lambda gene: (gene.getSubstrate(), gene.getBegin()))
fragLens = loadFragLens(fragLenFile)

# Create output files
OUT1 = gzip.open(outpath_refalt + "/" + outFile1, "wt")
OUT2 = gzip.open(outpath_refalt + "/" + outFile2, "wt")

# Simulate
print("simulating...", file=sys.stderr, flush=True)
nextReadID = 0
IN = open(samFile, "rt")
num_gene = 0
n_break = 0
# output list for debugging purpose
#
counter = 0
with open(output_path + "/simulation_info.txt", "w") as file_handler:
    file_handler.write(
        "chrN,geneid,n_while,num_variants,nth_reads,nth_transcript_picked,transcriptID,num_biSNP,REF_FWD,ALT_FWD,REF-quality,REF_REV,ALT_FWD,ALT_REV,ALT-quality\n"
    )
    for idx, gene in enumerate(genes):
        counter += 1
        num_gene += 1
        geneid = gene.getId()
        chrN = gene.getSubstrate()
        loadTranscriptSeqs(gene, genome2bit, twoBitDir)

        length = gene.longestTranscript().getLength()
        numReads = int(float(DEPTH / readLen) * length)
        # DEPTH : 10
        # length: gene length from longest transcript
        # readLen: 75
        # numReads = 10/75*(longest transcript length) = 0.133 * longest transcript length
        print("longest transcript length of this gene is %d" % (length))

        # all variants in gene region from VCF
        variants = getGeneVariants(gene, vcfFile)
        if len(variants) == 0:
            print("no variants")
            continue
        print(
            "%s,gene: %s,#reads: %s,total #variants in VCF: %d"
            % (chrN, geneid, numReads, len(variants))
        )
        maternal, transcriptIdToBiSNPcount, transcriptIdToBiSNPpos = makeAltCopy(
            gene, 1, variants
        )
        paternal, _, _ = makeAltCopy(gene, 0, variants)
        if isinstance(maternal, type(None)):
            print(
                "!!!!! DEBUG: skipping gene %s because it does not contain biSNP"
                % (geneid)
            )
            continue
        # initialization
        refSeq = None
        for i in range(numReads):
            # while loop incase
            n_while = 0
            while refSeq is None:
                n_while += 1
                (
                    matTranscript,
                    patTranscript,
                    th_transcript,
                    transcriptID,
                    bivariant,
                    transpos_chosen,
                ) = pickTranscript(
                    maternal, paternal, transcriptIdToBiSNPcount, transcriptIdToBiSNPpos
                )
                length = matTranscript.getLength()
                rec1 = nextSamRec(IN, samFile)
                rec2 = nextSamRec(IN, samFile)
                fragLen = fragLens[random.randrange(len(fragLens))]
                n_while_loop = 0
                if length < fragLen:
                    fragLen = length
                (refSeq, refSeq_rev, altSeq, altSeq_rev, qual1, qual2) = simRead(
                    n_while,
                    geneid,
                    nextReadID,
                    patTranscript,
                    matTranscript,
                    rec1,
                    rec2,
                    fragLen,
                    transpos_chosen,
                )

            file_handler.write(
                "%s,%s,%d,%d,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"
                % (
                    chrN,
                    geneid,
                    n_while,
                    len(variants),
                    i + 1,
                    th_transcript + 1,
                    transcriptID,
                    bivariant,
                    refSeq,
                    altSeq,
                    qual1,
                    refSeq_rev,
                    altSeq_rev,
                    qual2,
                )
            )
            # ref+alt
            printRead(
                "@SIM-" + str(nextReadID) + "-" + str(geneid) + ":REF:FWD" + " /1",
                refSeq,
                qual1,
                OUT1,
            )
            printRead(
                "@SIM-" + str(nextReadID) + "-" + str(geneid) + ":REF:REV" + " /2",
                refSeq_rev,
                qual2,
                OUT2,
            )
            nextReadID += 1
            # ref+alt
            printRead(
                "@SIM-" + str(nextReadID) + "-" + str(geneid) + ":ALT:FWD" + " /1",
                altSeq,
                qual1,
                OUT1,
            )
            printRead(
                "@SIM-" + str(nextReadID) + "-" + str(geneid) + ":ALT:REV" + " /2",
                altSeq_rev,
                qual2,
                OUT2,
            )
            nextReadID += 1

print("")
print("saved debugging lists!")
print(">> total geneID processed: %d" % (num_gene))
print(">> total ReadID processed: %d" % (nextReadID))
# print(">> total num breaks /gene is shorter than fragment length : %d" % (n_break))
# # matches=0 # number of sites containing alt or ref allele in transcript
# # mismatches=0 # number of sites having neither ref nor alt allele in transcript
# print("matches: %s , mismatches : %d" % (matches, mismatches))
# if mismatches == 0 and matches == 0:
#     print("No match or mismatch!")
# else:
#     matchRate = float(matches) / float(matches + mismatches)
#     print("%s percent of exonic sites had ref or alt " % (matchRate * 100))

# MY SAM FILE:
# HWI-ST661:130:C037KACXX:6:2307:8186:172941      99      chr1    13170   109     75M     =       13278   183     CTGTTGGGGAGGCAGCTGTAACTCAAAGCCTTAGCCTCTGTTCCCACGAAGGCAGGGCCATCAGGCACCAAAGGG     @@@DDFFFHHFFFJIIIIHECCCGHIICEHJCEBFHEGIGIIDGHBDGGHDEGG?HFFFFFEEED>ADDBDBDDB     RG:Z:0  NM:i:3  XT:A:U  md:Z:75 XA:Z:chrY,-59359485,75M,1;chrX,-155256479,75M,1;chr9,+13283,75M,1;chr2,-114357772,75M,1;chr15,-102517926,75M,2;chr12,-92369,75M,2;

# SCARLETT'S SAM FILE:
# ERR188231.12567592      163     chr10   93307   255     75M     =       93442   210     GCAAAGTAACTGCTGTTCTTATCTTGAATGTTGAACATTTGTTCATCCACCTCCCTCATGGGCATGCGACCCCTG     CCCFFFDFHHHHHJJJJIJJJJJJJJJIJIIIJJJJIJJJJIJJIJJIJJGIIJIIJIJIJJIJJJJIIIJHHFF     MD:Z:75 PG:Z:MarkDuplicates     NH:i:1  HI:i:1  jI:B:i,-1       NM:i:0  jM:B:c,-1       nM:i:1  AS:i:146
