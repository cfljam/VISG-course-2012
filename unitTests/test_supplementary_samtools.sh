#!/usr/bin/env bash

cd ../supplementary_samtools_usage

## Execute bash initialization script

./init.sh

## Execute bash script
./samtools.sh  2> log

echo "[ unit test finished ]\n" >&2


