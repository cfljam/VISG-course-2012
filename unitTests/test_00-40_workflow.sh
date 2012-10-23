#!/usr/bin/env bash

dir=$HOME/VISG-course-2012

## Run each bash script
echo "[ 05 unit test ]" >&2
cd $dir/05.reference/
time ./run.sh
cd -

echo "[ 10 unit test ]" >&2
cd $dir/10.bwa_sw_Pool1/
./run.sh
cd -

echo "[ 11 unit test ]" >&2
cd $dir/11.Samtools_Filter_Pool1
./run.sh
cd -

echo "[ 15 unit test ]" >&2
cd $dir/15.bwa_sw_Pool2
./run.sh
cd -

echo "[ 16 unit test ]" >&2
cd $dir/16.Samtools_Filter_Pool2
./run.sh
cd -

echo "[ 20 unit test ]" >&2
cd $dir/20.mpileup
./run.sh

echo "[ 25 unit test ]" >&2
cd $dir/25.bcftools
./run.sh
cd -

echo "[ 30 unit test ]" >&2
cd $dir/30.popoolation
./run.sh
cd -

echo "[ 35 unit test ]" >&2
cd $dir/35.Visualise_FET
Rscript FET_plots.R
cd -

echo "[ 40 unit test ]" >&2
cd $dir/40.igv/
./run.sh
cd -
