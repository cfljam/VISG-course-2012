#!/usr/bin/env bash

## Create a sorted bam file filtered with minimum quality of 20 (MAPQ column 5 in sam file) 
echo "using samtools to filter  sam files"
echo "creating a single filtered bam file"

samtools view -q 20 -Sb ../10.bwa_sw_Pool1/Pool1_Pop2.sam | samtools sort -  Pool1_Pop2.filtered

## Looping through all sam files
echo "now run on all populations"
for N in {2..8}
do  
  samtools view -q 20 -Sb ../10.bwa_sw_Pool1/Pool1_Pop$N.sam | samtools sort - Pool1_Pop$N.filtered
done
