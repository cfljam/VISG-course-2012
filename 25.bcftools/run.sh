#!/usr/bin/env bash

echo "Calling snps in Pool1"
samtools mpileup -ugf ../05.reference/pool1.fasta ../11.Samtools_Filter_Pool1/*.bam  | bcftools view -ve - > Pool1.vcf
echo "Calling snps in Pool2"
samtools mpileup -ugf ../05.reference/pool2.fasta ../16.Samtools_Filter_Pool2/*.bam  | bcftools view -ve - > Pool2.vcf

