#!/usr/bin/env sh

## Flow control using a for loop to run 'bwa index'
for ref in `ls *.fasta`
do  
  bwa index $ref
done
