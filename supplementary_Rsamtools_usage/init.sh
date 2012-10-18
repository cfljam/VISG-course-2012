#!/usr/bin/env bash

## Setup, clear out any existing files
rm *.sam *.bam *.bai 2> /dev/null

# Default type
type="roche"


if [ $# -gt 0 ]
then
    echo "[ options are: illumina (default), roche ]"
    type=$1
fi

echo "[ NGS type:$type ]"

if [ $type == 'roche' ] 
then
## Create symbolic link to 454 example sam file
    echo "[ Using 454 data ]" >&2
    ln -s ../10.bwa_sw_Pool1/Pool1_Pop2.sam aln.sam 
else
    echo "[ Simulating illumina paired end data ]" >&2    
    echo "[ mkdir bwa ]" >&2
    mkdir bwa 2> /dev/null

    ## change to 'bwa' directory

    echo "[ change t bwa directory ]" >&2
    cd bwa
    ln -s ../../05.reference/pool1.fasta pool1.fasta
    
    ## Run wgsim to randomly generate an illumina paired end dataset
    nreads=2000
    len=70
    echo "[ wgsim call ]\n" >&2
    echo "wgsim -N $nreads -1 $len -2 $len pool1.fasta read1.fq read2.fq > snips.txt" >&2
    wgsim -N $nreads -S 42 -1 $len -2 $len pool1.fasta read1.fq read2.fq > snips.txt

    ## Align using bwa 
    echo "[ bwa index ]" >&2
    bwa index -a is pool1.fasta 

    echo "[ bwa gapped alignments ]" >&2
    bwa aln pool1.fasta read1.fq > aln_sa1.sai 
    bwa aln pool1.fasta read2.fq > aln_sa2.sai 
    echo "[ bwa generate alignment ]" >&2
    bwa sampe pool1.fasta aln_sa1.sai aln_sa2.sai read1.fq read2.fq > ../aln.sam 
    ## change to parent directory
    cd -    
fi

exit;
