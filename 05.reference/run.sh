#!/bin/sh
for ref in `ls *.fasta`; do  bwa index $ref; done
