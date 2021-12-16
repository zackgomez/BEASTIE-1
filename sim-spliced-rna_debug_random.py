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


def simRead_random_patmat(
    refTranscript, altTranscript, rec1, rec2, fragLen, nextReadID, geneid
):
    random_prob = random.random()
    if_mat = False
    if random_prob >= 0.5:
        if_mat = True

    L = len(refTranscript.sequence)
    if L < fragLen or L < readLen:
        return (None, None, None, None, None)
    # L = 100
    # fragLen = 80
    lastStart = L - fragLen  # 20
    start1 = random.randrange(lastStart + 1)  # 10
    end1 = start1 + rec1.readLen  # 10+75 = 85
    end2 = start1 + fragLen  # 10+80 = 90
    start2 = end2 - rec2.readLen  # 90-75 = 15

    print(
        "@SIM-%s-%s, if maternal: %s, FWD seq length %d,REV seq length %d, FWD start-end : %d-%d, REV start-end : %d-%d"
        % (
            nextReadID,
            geneid,
            if_mat,
            rec1.readLen,
            rec2.readLen,
            start1,
            end1,
            start2,
            end2,
        )
    )
    refSeq = refTranscript.sequence[start1:end1]
    refSeq_rev = Seq(refTranscript.sequence[start2:end2]).reverse_complement()
    altSeq = altTranscript.sequence[start1:end1]
    altSeq_rev = Seq(altTranscript.sequence[start2:end2]).reverse_complement()
    if if_mat is True:
        randomseq = altSeq
        randomseq_rev = altSeq_rev
        print("@SIM-%s-%s,ALT FWD: %s" % (nextReadID, geneid, altSeq))
        print("@SIM-%s-%s,ALT REV: %s" % (nextReadID, geneid, altSeq_rev))
    else:
        randomseq = refSeq
        randomseq_rev = refSeq_rev
        print("@SIM-%s-%s,REF FWD: %s" % (nextReadID, geneid, refSeq))
        print("@SIM-%s-%s,REF REV: %s" % (nextReadID, geneid, refSeq_rev))
    print("")
    return (randomseq, randomseq_rev, rec1.qual, rec2.qual, if_mat)


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
with open(output_path + "/simulation_info.txt", "w") as file_handler:
    file_handler.write(
        "identifier_random,chrN,geneid,num_variants,nth_reads,nth_transcript_picked,transcriptID,num_biSNP,random_FWD,FWD-quality,random_REV,REV-quality\n"
    )

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
            (randomSeq, randomSeq_rev, qual1, qual2, if_mat,) = simRead_random_patmat(
                patTranscript,
                matTranscript,
                rec1,
                rec2,
                fragLen,
                nextReadID,
                geneid,
            )
            if randomSeq is None:
                n_break += 1
                continue  # gene is shorter than fragment length

            identifier_random = "@SIM-" + str(nextReadID) + "-" + str(geneid)
            file_handler.write(
                "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"
                % (
                    identifier_random,
                    str(chrN),
                    geneid,
                    str(len(variants)),
                    str(i + 1),
                    str(th_transcript + 1),
                    transcriptID,
                    bivariant,
                    randomSeq,
                    qual1,
                    randomSeq_rev,
                    qual2,
                )
            )
            # ref+alt
            printRead(
                "@SIM-"
                + str(nextReadID)
                + "-"
                + str(geneid)
                + ":"
                + str(if_mat)
                + ":FWD"
                + " /1",
                randomSeq,
                qual1,
                OUT1,
            )

            nextReadID += 1

            # ref+alt
            printRead(
                "@SIM-"
                + str(nextReadID)
                + "-"
                + str(geneid)
                + ":"
                + str(if_mat)
                + ":FWD"
                + " /1",
                randomSeq_rev,
                qual2,
                OUT2,
            )

            nextReadID += 1
        if processed_genes % 100 == 0:
            print(
                f"{datetime.now()} processed {processed_genes} / {len(gene_to_variants)} genes : {round(processed_genes/len(gene_to_variants),2)} %",
                file=sys.stderr,
                flush=True,
            )

    print(
        f"{datetime.now()} processed {processed_genes} / {len(gene_to_variants)} genes : {round(processed_genes/len(gene_to_variants),2)} %",
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
