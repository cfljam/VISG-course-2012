#!/usr/bin/env Rscript

require(ShortRead)
require(Rsamtools)
require(chipseq)


# Getting vignette help
if(interactive()) vignette("Rsamtools-Overview")


# Globals
if(interactive()) {
  PNG         <- FALSE
} else {
  PNG         <- TRUE
}
  
plotDir       <- "plots"
suppressWarnings(dir.create(plotDir))
basedir       <- "~/VISG-course-2012"
path          <- file.path(basedir, "supplementary_Rsamtools_usage")

# Set working directory
if(getwd() != path) {
  setwd(path)
}
samName       <- dir(path, pattern="\\.sam$")
fl            <- file.path(path, samName)
bamDest       <- "aln"


## sam to bam (also sorts and makes an .bai index file) 
bamName <- asBam(fl, bamDest, overwrite=TRUE)


## Sanity check:
dir(pattern="\\.ba[mi]$")


## Examining the bam file

## What columns are we interested in? See help(BamInput)
what <- scanBamWhat()


## Retrieving header information in a BAM file
ft    <- scanBamHeader(bamName)[[1]][["targets"]]
print(ft)


## Which features are we interested in extracting? -requires an indexed BAM file to exist
which <- GRanges(names(ft), IRanges(1, ft))


##  Create a parameter object for scanning BAM files, 
param <- ScanBamParam(which=which, what=what)


## Load sorted BAM file into R
bam <- scanBam(bamName, param = param)
object.size(bam)


## bam file: class, length, element classes
class(bam)
length(bam)
sapply(bam, class)


## Each element of the list corresponds to a range specified by the which argumen
for(i in 1:3) {
  cat("[", names(bam)[i], "]\n")
  print(class(bam[[i]]))
  print(length(bam[[i]]))
  print(sapply(bam[[i]], class))
}


## First bam[[1]] list component
print(names(bam)[1])


## Each component is a list containing the elements specified by the 'what'
for(j in seq_len(length(bam[[1]]))) {
  cat("[", names(bam[[1]])[j], "list element ]\n")
  print(head(bam[[i]][[j]], n=10))
  cat("...\n\n")
}


## Referencing by name or index is the same
identical(bam[["CO_Pool1_contig00004:1-1928"]], bam[[1]])


## Cigar string see help(cigar-utils) for utility functions
head(bam[[1]][["cigar"]])


## grep on cigar string containing INDEL characters I or D
noINDELS <- grep("[ID]", bam[[1]][["cigar"]], invert=TRUE)


## Reads not containing INDELS in alignment
print(bam[[1]][["cigar"]][noINDELS])


## Alphabet by cycle
abc <- alphabetByCycle(bam[[1]][["seq"]])

## Printing the first four alphabet cycles
abc[,1:4]

## Post alignment quality: Alphabet plot
if(PNG) png(file.path(plotDir, "alphabetPlot.png"))
matplot(t(abc[1:4, ]), type="l", lty=1, lwd=1, ylab="Nuclotide frequency",  main=names(bam)[1])
if(PNG) dev.off()


## Post alignment quality: Alphabet plot (first 50 cycles only)
if(PNG) png(file.path(plotDir, "alphabetPlot_50cycles.png"))
matplot(t(abc[1:4,1:50]), type="l", lty=1, lwd=1, ylab="Nuclotide frequency", main=names(bam)[1])
if(PNG) dev.off()


## Extracting +/- strand indices
indPos <- which(bam[[1]][["strand"]] == "+")
indNeg <- which(bam[[1]][["strand"]] == "-")


## ## Alphabet by cycle +/- strands
abc.P <- alphabetByCycle(bam[[1]][["seq"]][indPos])
abc.N <- alphabetByCycle(bam[[1]][["seq"]][indNeg])

## Visualize alphabet frequencies +/- strands
if(PNG) png(file.path(plotDir, "alphabetPlot_pos.png"))
matplot(t(abc.P[1:4,]), type="l", lty=1, lwd=1, ylab="Nuclotide frequency",  main=paste(names(bam)[1], "+ strand"))
if(PNG) dev.off()

if(PNG) png(file.path(plotDir, "alphabetPlot_neg.png"))
matplot(t(abc.N[1:4,]), type="l", lty=1, lwd=1, ylab="Nuclotide frequency",  main=paste(names(bam)[1], "- strand"))
if(PNG) dev.off()


## Construct a ShortReadQ object for QC
fq <- ShortReadQ(bam[[1]][["seq"]], bam[[1]][["qual"]], BStringSet(bam[[1]][["qname"]]))


## qa() report from bam file sequences
qaSummary <- qa(fq, lane="Roche 454")
perCycle  <- qaSummary[["perCycle"]]


## Generate a quality by cycle plot
if(PNG) png(file.path(plotDir, "cycleQuality.png"))
ShortRead:::.plotCycleQuality(perCycle$quality, main=names(bam)[1])
if(PNG) dev.off()


## Coverage plot ignoring strand orientation
IRanges   <- IRanges(start = bam[[1]][["pos"]], width=bam[[1]][["qwidth"]])
Cov       <- coverage(IRanges)
Peaks     <- slice(Cov, 0)

if(PNG) png(file.path(plotDir, "coverage.png"))
coverageplot(Peaks, main=names(bam)[1])
if(PNG) dev.off()


## Coverage plot +/- strands
IRangesF   <- IRanges(start = bam[[1]][["pos"]][indPos], width=bam[[1]][["qwidth"]][indPos])
CovF       <- coverage(IRangesF)
PeaksF     <- slice(CovF, 0)

IRangesR   <- IRanges(start = bam[[1]][["pos"]][indNeg], width=bam[[1]][["qwidth"]][indNeg])
CovR       <- coverage(IRangesR)
PeaksR     <- slice(CovR, 0)

## Coverage plot +/- strands visualized separately
if(PNG) png(file.path(plotDir, "coverageStrands.png"))
coverageplot(PeaksF, PeaksR,  main=names(bam)[1])
if(PNG) dev.off()
