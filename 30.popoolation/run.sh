#!/usr/bin/env sh

## Convert mpileup to sync with popoolation
perl ../popoolation2_1201/mpileup2sync.pl --fastq-type sanger  --input ../20.mpileup/Pool1.mpileup  --output Pool1.sync --min-qual 20
perl ../popoolation2_1201/mpileup2sync.pl --fastq-type sanger  --input ../20.mpileup/Pool2.mpileup  --output Pool2.sync --min-qual 20
cat Pool*.sync > combined.sync

## Run fishers exact test comparing snp frequencies among all populations
perl ../popoolation2_1201/fisher-test.pl --input combined.sync  --output combined.fet --min-count 6 --min-coverage 10 --max-coverage  10000 --suppress-noninformative

## Export test results to IGV format
perl ../popoolation2_1201/export/pwc2igv.pl  --input combined.fet  --output combined.fet.igv
