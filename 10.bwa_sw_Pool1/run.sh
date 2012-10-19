#!/usr/bin/env bash

## running 'bwa bwasw' heuristic Smith-Waterman-like alignment for long reads

echo "first check our command line"
for N in {2..8}; 
do 
  echo "bwa bwasw -f Pool1_Pop$N.sam ../05.reference/pool1.fasta ../00.raw/Pool1_BARCODE$N.fastq"
done

echo "now running for real"
for N in {2..8}; 
do 
  bwa bwasw -f Pool1_Pop$N.sam ../05.reference/pool1.fasta  ../00.raw/Pool1_BARCODE$N.fastq
done

#create an unfiltered  bam file just to compare sizes
samtools view -bS -o Pool1_Pop2_unfiltered.bam Pool1_Pop2.sam 
