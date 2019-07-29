#!/usr/bin/env python3
import os
import sys
import glob
import argparse
import itertools

def decodeKvPairs(kv):
    keyValues = {}
    for part in kv.split(';'):
        kv = part.strip().split()
        if len(kv)==2:
            key = kv[0]
            value = kv[1].replace('"', '')
            keyValues[key] =  value
    return keyValues

import collections


def exonGTF_to_intronGTF(exon_path):
    geneToExonRanges = collections.defaultdict(list) #(gene,chrom,strand)-> [range(start,end), range(start,end)]
    prevChrom = None
    with open(exon_path) as f:
        for line in f:
            if line.startswith('#'):
                yield line
                continue
            chromosome, source, feature_type, feature_start, feature_end, at5, strand, at7, kv = line.split(None,8)

            if feature_type!='exon':
                continue

            geneToExonRanges[(
                decodeKvPairs(kv)['transcript_id'],chromosome,strand,source
                )
                ].append( range(int(feature_start),int(feature_end)+1))

            if prevChrom is not None and prevChrom!=chromosome:
                for intr in generate_introns(geneToExonRanges):
                    yield intr
                geneToExonRanges = collections.defaultdict(list)


            prevChrom = chromosome

    if prevChrom is not None:
        for intr in generate_introns(geneToExonRanges):
            yield intr
        geneToExonRanges = collections.defaultdict(list)

def generate_introns(geneToExonRanges):
    for g in geneToExonRanges:
        if len(geneToExonRanges[g])==0:
            continue
        geneStart = min([min(r) for r in geneToExonRanges[g]])
        geneEnd = max([max(r) for r in geneToExonRanges[g]])
        gene_id, chromosome, strand,source = g

        introns = []
        in_intron = None
        current_intron_start = None

        for i in sorted(

[geneStart,geneEnd]+
                        [i.start for i in geneToExonRanges[g]]+
                        [i.stop for i in geneToExonRanges[g]]):
            if any([ i in r for r in geneToExonRanges[g]]):
                if current_intron_start is not None and in_intron is True :
                    introns.append([current_intron_start, i])

                    yield '\t'.join([ str(x) for x in [
                        chromosome,
                        source,
                        'intron',
                        current_intron_start,
                        i-1,
                        '.',
                        strand,
                        '.', f'gene_id "{gene_id}";']])+'\n'

                in_intron = False
            else:
                if in_intron is False:
                    in_intron = True
                    current_intron_start=i

if __name__=='__main__':
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Convert exonic GTF file to intronic GTF file")
    argparser.add_argument('gffin')
    argparser.add_argument('-o', help="Output GTF file", type=str, required=True)
    args = argparser.parse_args()


    with open(args.o, 'w') as o:
        for intron in exonGTF_to_intronGTF(args.gffin):
            if intron is not None:
                o.write(intron)