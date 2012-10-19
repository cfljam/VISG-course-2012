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


######################################################################
## 1. Path to first example Fastq file
######################################################################
firstFile <- file.path(indir, fileNames[1])


######################################################################
## 2. Read in Fastq example
######################################################################
rfq       <- readFastq(firstFile, qualityType="Auto")


######################################################################
## 3. Show method
######################################################################
print(rfq)


######################################################################
## 4. Class of rfq is an instance of ShortReadQ S4 class
######################################################################
class(rfq)


######################################################################
## 5. Print all S4 class slotNames
######################################################################
slotNames(rfq)


######################################################################
## 6. sread accessor method
######################################################################
sread(rfq)


######################################################################
## 7. id accessor method
######################################################################
id(rfq)


######################################################################
## 8. quality accessor method
######################################################################
quality(rfq)


######################################################################
## 9. Run Quality control script
######################################################################
qaSummary <- qa(rfq, lane="Roche 454")


######################################################################
## 10. Create hmtml report
######################################################################
fname <- report(qaSummary, dest="/tmp")


######################################################################
## 11. Open html report 
######################################################################
if(interactive()) browseURL(fname)


######################################################################
## 12. Snippet to Investigate qaSummary structure
######################################################################
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


######################################################################
## 13. Generate qa reports for all Fastq files
######################################################################
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

  
