#!/usr/bin/env Rscript

# Clear working directory objects
rm(list=ls())

require(ShortRead)

# Getting vignette help
if(interactive()) vignette("Overview") ## Quality assessment, page 8


# Globals
cat("[ Code to source() ]\n")
cat(readLines("globals.R"), sep="\n")

# source globals code
source("globals.R")

cat("[ Directory and Fastq filenames ]\n")
print(indir)
print(fileNames)


## Path to first example Fastq file
firstFile <- file.path(indir, fileNames[1])


## Read in Fastq example
rfq       <- readFastq(firstFile, qualityType="Auto")


## Show method
print(rfq)


## Class of rfq is an instance of ShortReadQ S4 class
class(rfq)


## Print all S4 class slotNames
slotNames(rfq)


## sread accessor method
sread(rfq)


## id accessor method
id(rfq)


## quality accessor method
quality(rfq)


## Run Quality control script
qaSummary <- qa(rfq, lane="Roche 454")


## Create hmtml report
fname <- report(qaSummary, dest="/tmp")


## Open html report 
if(interactive()) browseURL(fname)


## Snippet to Investigate qaSummary structure
cat("[ qaSummary structure ]\n")
for(i in names(qaSummary)) {
  if(i == "perTile") next;
  if(i == "perCycle") {
    for(j in names(qaSummary[[i]])) {
      cat("[", i, "][", j, "]\n")
      print(head(qaSummary[[i]][[j]]))
      cat("...\n\n")
    }
  } else {
    cat("[", i, "]\n")
    print(head(qaSummary[[i]]))
    cat("...\n\n")
  }
}


## Generate qa reports for all Fastq files
for(f in fileNames) {
  # Truncate the ".fastq$" off names
  fdir      <- gsub("\\..+?$", "", f)
  # Generate a qc object
  qaSummary <- qa(indir, pattern=f, lane=fdir, type="fastq")
  # Generate an html report
  fname     <- report(qaSummary, dest=file.path(path, fdir))
  # Be verbose
  cat("[", fdir , "=>", fname, "]\n")
}

  
