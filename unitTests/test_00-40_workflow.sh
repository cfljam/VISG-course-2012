#!/usr/bin/env bash

dir=$HOME/VISG-course-2012

## Run each bash script
echo "[ 05 unit test ]" >&2
$dir/05.reference/run.sh

echo "[ 10 unit test ]" >&2
$dir/10.bwa_sw_Pool1/run.sh

echo "[ 11 unit test ]" >&2
$dir/11.Samtools_Filter_Pool1/run.sh

echo "[ 15 unit test ]" >&2
$dir/15.bwa_sw_Pool2/run.sh

echo "[ 16 unit test ]" >&2
$dir/16.Samtools_Filter_Pool2/run.sh

echo "[ 20 unit test ]" >&2
$dir/20.mpileup/run.sh

echo "[ 25 unit test ]" >&2
$dir/25.bcftools/run.sh

echo "[ 30 unit test ]" >&2
$dir/30.popoolation/run.sh

echo "[ 35 unit test ]" >&2
cd $dir/35.Visualise_FET
Rscript FET_plots.R
cd -

echo "[ 40 unit test ]" >&2
$dir/40.igv/run.sh