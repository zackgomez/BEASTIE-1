#!/usr/bin/env python
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

from splice_reads import (
    constructTwoBitInput,
    run2bitBatch,
    annotate_transcripts,
    makeAltCopy,
    pickTranscript,
    nextSamRec,
    loadFragLens,
    printRead,
    loadTranscriptSeqs,
    readFrom2bit,
    Variant,
    SamRecord,
    chunk_iter,
)

# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import ProgramName
import random
import copy
import gzip
import os
from collections import defaultdict
from GffTranscriptReader import GffTranscriptReader
from Pipe import Pipe
from ConfigFile import ConfigFile
from Translation import Translation
from Rex import Rex
from Bio.Seq import Seq
from pathlib import Path
from datetime import datetime

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

def simRead_patmat(
    refTranscript, altTranscript, rec1, rec2, fragLen, nextReadID, geneid,
):
    #####################################################################
    L = len(refTranscript.sequence)
    if L < fragLen or L < readLen or L < len(rec1.qual) or L < len(rec2.qual):
        return (None, None, None, None, None,None,None,None)
    # L = 100
    # fragLen = 80
    
    lastStart = min(L - fragLen, L - len(rec1.qual),L - len(rec2.qual))  # 20

    start1 = random.randrange(lastStart + 1)  # 10
    end1 = start1 + len(rec1.qual)  #rec1.readLen  # 10+75 = 85
    LEN1 = end1-start1
    end2 = start1 + fragLen  # 10+80 = 90
    start2 = end2 - len(rec2.qual) #rec2.readLen  # 90-75 = 15
    assert start1 >= 0
    assert end1 <= L
    assert start2 >= 0
    assert end2 <= L
    LEN2 = end2-start2

    ######## forward strand, same sequence pos for mat/aptf fov
    refSeq = refTranscript.sequence[start1:end1]
    altSeq = altTranscript.sequence[start1:end1]
    ######## reverse strand, same sequence pos for mat/apt rev
    refSeq_rev = Seq(refTranscript.sequence[start2:end2]).reverse_complement()
    altSeq_rev = Seq(altTranscript.sequence[start2:end2]).reverse_complement()

    return (refSeq, refSeq_rev, altSeq, altSeq_rev, rec1.qual, rec2.qual,LEN1,LEN2)
    # print(
    #     "@SIM-%s-%s, if maternal: %s, FWD seq length %d,REV seq length %d, FWD start-end : %d-%d, REV start-end : %d-%d"
    #     % (
    #         nextReadID,
    #         geneid,
    #         if_mat,
    #         rec1.readLen,
    #         rec2.readLen,
    #         start1,
    #         end1,
    #         start2,
    #         end2,
    #     )
    # )

# =========================================================================
# main()
# =========================================================================
if len(sys.argv) != 7:
    exit(
        ProgramName.get()
        + " <config-file> <per-base-read-depth> <if_random> <prefix> <out-read1.gz> <out-read2.gz>\n"
    )
(configFile, DEPTH, if_random, prefix, outFile1, outFile2) = sys.argv[1:]

if_random = if_random == "True"
DEPTH = int(DEPTH)
# if if_random:
#     DEPTH = int(DEPTH)
# else:
#     DEPTH = int(DEPTH) / 2

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
outpath_refalt = output_path + "/simulated_fastq"
Path(output_path).mkdir(parents=True, exist_ok=True)
Path(outpath_refalt).mkdir(parents=True, exist_ok=True)

# Load GFF and fragment lengths
gffReader = GffTranscriptReader()
print(f"{datetime.now()} reading GFF...", file=sys.stderr, flush=True)
genes = gffReader.loadGenes(gffFile)
genes.sort(key=lambda gene: (gene.getSubstrate(), gene.getBegin()))
print(f"{datetime.now()} done reading GFF...", file=sys.stderr, flush=True)
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
# with open(output_path + "/simulation_info.txt", "w") as file_handler:
#     file_handler.write(
#         "identifier_random,chrN,geneid,num_variants,nth_reads,nth_transcript_picked,transcriptID,num_biSNP,random_FWD,FWD-quality,random_REV,REV-quality\n"
#     )

print(
    f"{datetime.now()} processing {len(genes)} genes for transcripts",
    file=sys.stderr,
    flush=True,
)

twoBitInputFile = "/tmp/twobitinput"
constructTwoBitInput(genes, twoBitInputFile)

print(
    f"{datetime.now()} running 2bit",
    file=sys.stderr,
    flush=True,
)
twobitId_to_seq = run2bitBatch(twoBitDir, twoBitInputFile, genome2bit)

print(
    f"{datetime.now()} done running 2bit {len(twobitId_to_seq)} ids",
    file=sys.stderr,
    flush=True,
)

annotate_transcripts(genes, twobitId_to_seq)

print(
    f"{datetime.now()} done annotating transcripts",
    file=sys.stderr,
    flush=True,
)
region_str_to_genes = {}
for idx, gene in enumerate(genes):
    counter += 1
    num_gene += 1

    region_str = f"{gene.getSubstrate()}:{gene.getBegin()}-{gene.getEnd()}"
    if region_str in region_str_to_genes:
        region_str_to_genes[region_str].append(gene)
    else:
        region_str_to_genes[region_str] = [gene]

print(
    f"{datetime.now()} Got {len(region_str_to_genes)} regions to tabix",
    file=sys.stderr,
    flush=True,
)

CHUNK_SIZE = 1000
regions_processed = 0
gene_to_variants = {}
for x in chunk_iter(iter(region_str_to_genes.keys()), CHUNK_SIZE):
    regions = " ".join(x)
    cmd = f"tabix --separate-regions {vcfFile} {regions}"
    output = Pipe.run(cmd)
    if len(output) == 0:
        continue

    lines = output.split("\n")
    for line in lines:
        if len(line) == 0:
            continue
        if line.startswith("#"):
            region_str = line[1:]
            if region_str not in region_str_to_genes:
                print(
                    f"{datetime.now()} bad region_str: '{region_str}'",
                    file=sys.stderr,
                    flush=True,
                )
            assert region_str in region_str_to_genes
            genes = region_str_to_genes[region_str]
            records = []
            for gene in genes:
                gene_to_variants[gene] = records
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
    regions_processed += len(x)
    print(
        f"{datetime.now()} finished {regions_processed} / {len(region_str_to_genes)} tabix regions",
        file=sys.stderr,
        flush=True,
    )

print(
    f"{datetime.now()} Got {len(gene_to_variants)} genes with variants",
    file=sys.stderr,
    flush=True,
)

processed_genes = 0
for (gene, variants) in gene_to_variants.items():
    chrN = gene.getSubstrate()
    geneid = gene.getId()
    length = gene.longestTranscript().getLength()
    numReads = int(float(DEPTH / readLen) * length)

    processed_genes += 1

    if len(variants) == 0:
        print(f"{chrN},gene: {geneid}, no variants, skip")
        continue
    print(
        "%s,gene: %s,#reads: %s,total #variants in VCF: %d"
        % (chrN, geneid, numReads, len(variants))
    )
    maternal, transcriptIdToBiSNPcount, transcriptIdToBiSNPpos = makeAltCopy(
        gene, 1, variants
    )
    paternal, _, _ = makeAltCopy(gene, 0, variants)
    # if isinstance(maternal, type(None)):
    #     print(
    #         "!!!!! DEBUG: skipping gene %s because it does not contain biSNP"
    #         % (geneid)
    #     )
    #     continue
    # # file_handler.write("%s,%s,%s,%d\n" % (chrN, geneid, numReads, len(variants)))
    for i in range(numReads):
        (
            matTranscript,
            patTranscript,
            th_transcript,
            transcriptID,
            bivariant,
        ) = pickTranscript(maternal, paternal, transcriptIdToBiSNPpos)
        length = matTranscript.getLength()
        rec1 = nextSamRec(IN, samFile)
        rec2 = nextSamRec(IN, samFile)
        fragLen = fragLens[random.randrange(len(fragLens))]
        if length < fragLen:
            fragLen = length

        (patSeq, patSeq_rev, matSeq, matSeq_rev, fwd_qual, rev_qual,fwd_LEN,rev_LEN) = simRead_patmat(
            patTranscript,
            matTranscript,
            rec1,
            rec2,
            fragLen,
            nextReadID,
            geneid,
            )

        if patSeq is None or matSeq is None:
            n_break += 1
            continue  # gene is shorter than fragment length
        
        ################## random or equal

        if if_random:
            random_prob = random.random()
            if_mat = False
            if random_prob >= 0.5:
                if_mat = True
            identifier_random = "@SIM-" + str(nextReadID) + "-" + str(geneid)
            print("%s,if maternal: %s,rec1 quality string length %s, forward strand length %s, rec2 quality strand length %s, reverse strand length %s!! "
                    %(identifier_random,
                    str(if_mat),
                    len(fwd_qual),
                    fwd_LEN,
                    len(rev_qual),
                    rev_LEN,
                )
            ) 
            ######## random decision on maternal / paternal
            if if_mat is True:
                randomSeq = matSeq
                randomSeq_rev = matSeq_rev
                print("%s,if mat:%s,MAT FWD: %s" % (identifier_random, str(if_mat),matSeq))
                print("%s,if mat:%s,MAT REV: %s" % (identifier_random, str(if_mat),matSeq_rev))
            else:
                randomSeq = patSeq
                randomSeq_rev = patSeq_rev
                print("%s,if mat:%s,PAT FWD: %s" % (identifier_random, str(if_mat),patSeq))
                print("%s,if mat:%s,PAT REV: %s" % (identifier_random, str(if_mat),patSeq_rev))
            print("%s,qual FWD: %s" % (identifier_random, fwd_qual))
            print("%s,qual REV: %s" % (identifier_random, rev_qual))
            print("")

            printRead(
                str(identifier_random)
                + ":if_maternal"
                + str(if_mat)
                + ":FWD"
                + " /1",
                randomSeq,
                fwd_qual,
                OUT1,
            )

            nextReadID += 1

            printRead(
                str(identifier_random)
                + ":if_maternal"
                + str(if_mat)
                + ":REV"
                + " /1",
                randomSeq_rev,
                rev_qual,
                OUT2,
            )
        else:
            identifier_PAT = "@SIM-" + str(nextReadID) + "-" + str(geneid)
            identifier_MAT = "@SIM-" + str(nextReadID + 1) + "-" + str(geneid)
            print("REF: %s,ALT: %s,rec1 quality string length %s, forward strand length %s, rec2 quality strand length %s, reverse strand length %s!! "
                    %(identifier_PAT,
                    identifier_MAT,
                    len(fwd_qual),
                    fwd_LEN,
                    len(rev_qual),
                    rev_LEN,
                )
            )
            print("%s,PAT FWD: %s" % (identifier_PAT, patSeq))
            print("%s,MAT FWD: %s" % (identifier_MAT, matSeq))
            print("%s,qual FWD: %s" % (identifier_PAT, fwd_qual))

            print("%s,PAT REV: %s" % (identifier_PAT, patSeq_rev))
            print("%s,MAT REV: %s" % (identifier_MAT, matSeq_rev))
            print("%s,qual REV: %s" % (identifier_MAT, rev_qual))

            print("")

            printRead(
                identifier_PAT + ":PAT:FWD" + " /1",
                patSeq,
                fwd_qual,
                OUT1,
            )
            printRead(
                identifier_PAT + ":PAT:REV" + " /2",
                patSeq_rev,
                rev_qual,
                OUT2,
            )

            nextReadID += 1

            printRead(
                identifier_MAT + ":MAT:FWD" + " /1",
                matSeq,
                fwd_qual,
                OUT1,
            )
            printRead(
                identifier_MAT + ":MAT:REV" + " /2",
                matSeq_rev,
                rev_qual,
                OUT2,
            )
        # file_handler.write(
        #     "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"
        #     % (
        #         identifier_random,
        #         str(chrN),
        #         geneid,
        #         str(len(variants)),
        #         str(i + 1),
        #         str(th_transcript + 1),
        #         transcriptID,
        #         bivariant,
        #         randomSeq,
        #         qual1,
        #         randomSeq_rev,
        #         qual2,
        #     )
        # )

        nextReadID += 1
    if processed_genes % 100 == 0:
        print(
            f"{datetime.now()} processed {processed_genes} / {len(gene_to_variants)} genes : {round(100*processed_genes/len(gene_to_variants),2)} %",
            file=sys.stderr,
            flush=True,
        )

    print(
        f"{datetime.now()} processed {processed_genes} / {len(gene_to_variants)} genes : {round(100*processed_genes/len(gene_to_variants),2)} %",
        file=sys.stderr,
        flush=True,
    )

print(
    f"{datetime.now()} DONE",
    file=sys.stderr,
    flush=True,
)

print("")
print("saved debugging lists!")
print(">> total geneID processed: %d" % (num_gene))
print(">> total ReadID processed: %d" % (nextReadID))
print(">> total num breaks /gene is shorter than fragment length : %d" % (n_break))
print("matches: %s , mismatches : %d" % (matches, mismatches))
if mismatches == 0 and matches == 0:
    print("No match or mismatch!")
else:
    matchRate = float(matches) / float(matches + mismatches)
    print("%s percent of exonic sites had ref or alt " % (matchRate * 100))
