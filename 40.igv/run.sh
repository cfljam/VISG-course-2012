#!/bin/sh
samtools cat -o Pool1.bam ../11.Samtools_Filter_Pool1/*.bam
samtools sort Pool1.bam Pool1.sorted
samtools index Pool1.sorted.bam Pool1.sorted.bai
samtools cat -o Pool2.bam ../16.Samtools_Filter_Pool2/*.bam
samtools sort Pool2.bam Pool2.sorted
samtools index Pool2.sorted.bam Pool2.sorted.bai
