#!/usr/bin/env Rscript

require(ShortRead)
require(Rsamtools)
require(chipseq)


## Getting vignette help
if(interactive()) vignette("Rsamtools-Overview")


## Setup
if(getwd() != "/VISG/03.Rsamtools_usage") {
  setwd("/VISG/03.Rsamtools_usage")
}


## Globals
PNG           <- TRUE
plotDir       <- "plots"

suppressWarnings(dir.create(plotDir))

basedir       <- "/VISG"
path <- file.path(basedir, "03.Rsamtools_usage")
samName       <- dir(path, pattern="\\.sam$")
fl            <- file.path(path, samName)
bamDest       <- "aln"


## sam to bam (also sorts and makes an .bai index file) 
bamName <- asBam(fl, bamDest, overwrite=TRUE)
## Sanity check:
dir(pattern="\\.ba[mi]$")


## Examining the bam file

## What columns are we interested in?
what <- scanBamWhat()

## Retrieving header information in a BAM file
ft    <- scanBamHeader(bamName)[[1]][["targets"]]
print(ft)

## Extract features 
## Which features are we interested in?
which <- GRanges(names(ft), IRanges(1, ft))
## which    <- RangesList(IRanges(1, ft[1]))
## names(which)<- ft[1]

##  Parameters for scanning BAM files, 'which' requires an indexed BAM file to exist
param <- ScanBamParam(which=which, what=what)

## Load sorted BAM file into R
bam <- scanBam(bamName, param = param)

object.size(bam)
str(bam)

## Details
class(bam)
length(bam)
sapply(bam, class)


## Each component is a list - see other tutorial for correct comments
for(i in 1:3) {
  cat("[", names(bam)[i], "]\n")
  print(class(bam[[i]]))
  print(length(bam[[i]]))
  print(sapply(bam[[i]], class))
}


## Examining bam[[1]] list elements
for(j in seq_len(length(bam[[1]]))) {
  cat("[", names(bam[[1]])[j], "list element ]\n")
  print(head(bam[[i]][[j]], n=10))
  cat("...\n\n")
}


## First bam list component
print(names(bam)[1])


## Referencing by name or index identical
identical(bam[["CO_Pool1_contig00004:1-1928"]], bam[[1]])


## Cigar string
head(bam[[1]][["cigar"]])


## grep on cigar string containing INDEL characters I or D
noINDELS <- grep("[ID]", bam[[1]][["cigar"]], invert=TRUE)


## Reads not containing INDELS in alignment
print(bam[[1]][["cigar"]][noINDELS])


## Post alignment quality
abc <- alphabetByCycle(bam[[1]][["seq"]])


## Alphabet plot
if(PNG) png(file.path(plotDir, "alphabetPlot_50cycles.png"))
matplot(t(abc[1:4,1:50]), type="l", lty=1, lwd=1, ylab="Nuclotide frequency", main=names(bam)[1])
if(PNG) dev.off()

if(PNG) png(file.path(plotDir, "alphabetPlot.png"))
matplot(t(abc[1:4, ]), type="l", lty=1, lwd=1, ylab="Nuclotide frequency",  main=names(bam)[1])
if(PNG) dev.off()


## Positive and negative strand indices
indPos <- which(bam[[1]][["strand"]] == "+")
indNeg <- which(bam[[1]][["strand"]] == "-")

## Alphabet frequencies
abc <- alphabetByCycle(bam[[1]][["seq"]][indPos])

## Visualize alphabet frequencies
if(PNG) png(file.path(plotDir, "alphabetPlot_pos.png"))
matplot(t(abc[1:4,]), type="l", lty=1, lwd=1, ylab="Nuclotide frequency",  main=paste(names(bam)[1], "+ strand"))
if(PNG) dev.off()


## Construct a ShortReadQ object
fq <- ShortReadQ(bam[[1]][["seq"]], bam[[1]][["qual"]], BStringSet(bam[[1]][["qname"]]))


## QA report on Fastq
qaSummary <- qa(fq, lane="Roche 454")
perCycle  <- qaSummary[["perCycle"]]


## Generate a Quality by cycle plot
if(PNG) png(file.path(plotDir, "cycleQuality.png"))
ShortRead:::.plotCycleQuality(perCycle$quality, main=names(bam)[1])
if(PNG) dev.off()


## Coverage plot - Both strands combined
IRanges   <- IRanges(start = bam[[1]][["pos"]], width=bam[[1]][["qwidth"]])
Cov       <- coverage(IRanges)
Peaks     <- slice(Cov, 0)

if(PNG) png(file.path(plotDir, "coverage.png"))
coverageplot(Peaks, main=names(bam)[1])
if(PNG) dev.off()

## Coverage plot - Positive and Negative strands
IRangesF   <- IRanges(start = bam[[1]][["pos"]][indPos], width=bam[[1]][["qwidth"]][indPos])
CovF       <- coverage(IRangesF)
PeaksF     <- slice(CovF, 0)

IRangesR   <- IRanges(start = bam[[1]][["pos"]][indNeg], width=bam[[1]][["qwidth"]][indNeg])
CovR       <- coverage(IRangesR)
PeaksR     <- slice(CovR, 0)

## Coverage plot - Separate strands
if(PNG) png(file.path(plotDir, "coverageStrands.png"))
coverageplot(PeaksF, PeaksR,  main=names(bam)[1])
if(PNG) dev.off()
