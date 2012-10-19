#!/usr/bin/env bash

## Call variants with mpileup
samtools mpileup  -d 10000 -f ../05.reference/pool1.fasta ../11.Samtools_Filter_Pool1/*.bam > Pool1.mpileup
samtools mpileup  -d 10000 -f ../05.reference/pool2.fasta ../16.Samtools_Filter_Pool2/*.bam > Pool2.mpileup
