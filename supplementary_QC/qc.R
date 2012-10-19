#!/usr/bin/env Rscript

require(ShortRead)

## Globals
basedir   <- "~/VISG-course-2012"
path      <- file.path(basedir, "supplementary_QC")
if(getwd() != path) {
  setwd(path)
}

indir     <- file.path(basedir, "00.raw")
outdir    <- file.path(basedir, "01.qc")
fileNames <- dir(path, pattern="fastq$")
suppressWarnings(dir.create(outdir))

## Read in Fastq example
rfq       <- readFastq(file.path(indir, fileNames[1]), qualityType="Auto")

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

## Investigating qaSummary content (include??)
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
  ## Truncate the ".fastq$" off names
  fdir      <- gsub("\\..+?$", "", f)
  ## Generate a qc object
  qaSummary <- qa(indir, pattern=f, lane=fdir, type="fastq")
  ## Generate an html report
  fname     <- report(qaSummary, dest=file.path(outdir, fdir))
  ## Be verbose
  cat("[", fdir , "=>", fname, "]\n")
}

  
