#!/bin/sh
echo "using samtools to filter  sam files"
echo "creating a sinlgle filtered bam file"
 samtools view -q 20 -Sb ../10.bwa_sw_Pool1/Pool1_Pop2.sam | samtools sort -  Pool1_Pop2_filtered
echo "now run on all populations"
 for N in {2..8}; do  samtools view -q 20 -Sb ../10.bwa_sw_Pool1/Pool1_Pop$N.sam | samtools sort - Pool1_Pop$N.filtered ; done
