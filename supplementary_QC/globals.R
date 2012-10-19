#!/usr/bin/env echo source() me in R
basedir   <- "~/VISG-course-2012"
path      <- file.path(basedir, "supplementary_QC")

# Set working directory
if(getwd() != path) {
  setwd(path)
}

indir     <- file.path(basedir, "00.raw")
fileNames <- dir(indir, pattern="fastq$")

