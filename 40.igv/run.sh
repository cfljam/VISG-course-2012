#!/usr/bin/env sh

## Concatenate bams into Pool1.bam
samtools cat -o Pool1.bam ../11.Samtools_Filter_Pool1/*.bam

## Sort Pool1.bam
samtools sort Pool1.bam Pool1.sorted

## Index Pool1.bam
samtools index Pool1.sorted.bam Pool1.sorted.bai

## Concatenate bams into Pool2.bam
samtools cat -o Pool2.bam ../16.Samtools_Filter_Pool2/*.bam

## Sort Pool2.bam
samtools sort Pool2.bam Pool2.sorted

## Index Pool2.bam
samtools index Pool2.sorted.bam Pool2.sorted.bai
