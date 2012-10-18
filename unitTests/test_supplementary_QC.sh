#!/usr/bin/env bash

cd ../supplementary_QC

## Execute R script

Rscript ./qc.R > log

echo "[ unit test finished ]\n" >&2


