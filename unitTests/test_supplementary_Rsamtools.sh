#!/usr/bin/env bash

cd ../supplementary_Rsamtools_usage

## Execute bash initialization script

./init.sh

## Execute bash script
Rscript Rsamtools.R  2> log

echo "[ unit test finished ]" >&2


