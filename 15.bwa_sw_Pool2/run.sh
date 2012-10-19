#!/bin/sh

## Constructing 'bwa bwasw' heuristic Smith-Waterman-like alignment command for long reads
echo "first check our command line"
for N in {2..8}
do 
  echo "bwa bwasw -f Pool2_Pop$N.sam../05.reference/pool2.fasta ../00.raw/Pool2_BARCODE$N.fastq"
done

## Running 'bwa bwasw' heuristic Smith-Waterman-like alignment for long reads
echo "now running for real"
for N in {2..8}
  do bwa bwasw -f Pool2_Pop$N.sam ../05.reference/pool2.fasta  ../00.raw/Pool2_BARCODE$N.fastq
done

## Create an unfiltered  bam file just to compare sizes
samtools view -bS -o Pool2_Pop2_unfiltered.bam Pool2_Pop2.sam
