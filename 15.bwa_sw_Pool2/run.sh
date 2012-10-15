#!/bin/sh
echo "first check our command line"
for N in {2..8}; do echo "bwa bwasw ../05.reference/pool2.fasta ../00.raw/Pool2_BARCODE$N.fastq -f  Pool2_Pop$N.sam";done
echo "now running for real"
for N in {2..8}; do bwa bwasw ../05.reference/pool2.fasta  ../00.raw/Pool2_BARCODE$N.fastq -f  Pool2_Pop$N.sam;done

#create an unfiltered  bam file just to compare sizes
samtools view -bS Pool2_Pop2.sam -o Pool2_Pop2_unfiltered.bam
