#!/usr/bin/env Rscript

## Clear working directory objects
rm(list=ls())

require(ShortRead)
require(Rsamtools)
require(chipseq)

# Getting vignette help
if(interactive()) vignette("Rsamtools-Overview")


# Globals
cat("[ Code to source() ]\n")
cat(readLines("globals.R"), sep="\n")

## source globals code
source("globals.R")

cat("[ path to bam:",fl, "]\n")
print(fl)


######################################################################
## 1. sam to bam (also sorts and makes an .bai index file) 
######################################################################
bamDest       <- "aln"
bamName <- asBam(fl, bamDest, overwrite=TRUE)


######################################################################
## Sanity check:
######################################################################
dir(pattern="\\.ba[mi]$")


######################################################################
## 2. Examining the bam file: w hat columns are we interested in? 
######################################################################
what <- scanBamWhat() ## See help(scanBam)


######################################################################
## 3. Retrieving header information in a BAM file
######################################################################
ft    <- scanBamHeader(bamName)[[1]][["targets"]]
print(ft)


######################################################################
## 4. Which features to extract? -requires an indexed BAM 
######################################################################
which <- GRanges(names(ft), IRanges(1, ft))


######################################################################
## 5. Create a parameter object for scanning BAM files, 
######################################################################
param <- ScanBamParam(which=which, what=what)


######################################################################
## 6. Load sorted BAM file into R
######################################################################
bam <- scanBam(bamName, param = param)
object.size(bam)


######################################################################
## 7. bam file: class, length, element classes
######################################################################
class(bam)
length(bam)
sapply(bam, class)


######################################################################
## 8. Each element of the list corresponds to a range specified by 'which'
######################################################################
names(bam)

######################################################################
## 9. First bam[[1]] list component
######################################################################
class(bam[[1]])
length(bam[[1]])
sapply(bam[[1]], class)


######################################################################
## 10. Each component is a list containing the elements specified by 'what'
######################################################################
for(j in seq_len(length(bam[[1]]))) {
  cat("[", names(bam[[1]])[j], "list element ]\n")
  print(head(bam[[1]][[j]], n=10))
  cat("...\n\n")
}


######################################################################
## 11. Referencing by name or index is the same
######################################################################
identical(bam[["CO_Pool1_contig00004:1-1928"]], bam[[1]])


######################################################################
## 12. Cigar string see help(cigar-utils) for utility functions
######################################################################
head(bam[[1]][["cigar"]])


######################################################################
## 13. grep on cigar string containing INDEL characters I or D
######################################################################
noINDELS <- grep("[ID]", bam[[1]][["cigar"]], invert=TRUE)


######################################################################
## 14. Reads not containing INDELS in alignment
######################################################################
print(bam[[1]][["cigar"]][noINDELS])


######################################################################
## 15. Alphabet by cycle
######################################################################
abc <- alphabetByCycle(bam[[1]][["seq"]])


######################################################################
## 16. Printing the first four alphabet cycles
######################################################################
print(abc[,1:4])


######################################################################
## 17. Post alignment quality: Alphabet plot
######################################################################
if(PNG) png(file.path(plotDir, "alphabetPlot.png"))
matplot(t(abc[1:4, ]), type="l", lty=1, lwd=1, ylab="Nuclotide frequency",  main=names(bam)[1])
if(PNG) dev.off()


######################################################################
## 18. Post alignment quality: Alphabet plot (first 50 cycles only)
######################################################################
if(PNG) png(file.path(plotDir, "alphabetPlot_50cycles.png"))
matplot(t(abc[1:4,1:50]), type="l", lty=1, lwd=1, ylab="Nuclotide frequency", main=names(bam)[1])
if(PNG) dev.off()


######################################################################
## 19. Extracting +/- strand indices
######################################################################
indPos <- which(bam[[1]][["strand"]] == "+")
indNeg <- which(bam[[1]][["strand"]] == "-")


######################################################################
## 20. Alphabet by cycle +/- strands
######################################################################
abc.P <- alphabetByCycle(bam[[1]][["seq"]][indPos])
abc.N <- alphabetByCycle(bam[[1]][["seq"]][indNeg])


######################################################################
## 21. Visualize alphabet frequencies +/- strands
######################################################################
if(PNG) png(file.path(plotDir, "alphabetPlot_pos.png"))
matplot(t(abc.P[1:4,]), type="l", lty=1, lwd=1, ylab="Nuclotide frequency",  main=paste(names(bam)[1], "+ strand"))
if(PNG) dev.off()

if(PNG) png(file.path(plotDir, "alphabetPlot_neg.png"))
matplot(t(abc.N[1:4,]), type="l", lty=1, lwd=1, ylab="Nuclotide frequency",  main=paste(names(bam)[1], "- strand"))
if(PNG) dev.off()


######################################################################
## 22. Construct a ShortReadQ object for QC
######################################################################
fq <- ShortReadQ(bam[[1]][["seq"]], bam[[1]][["qual"]], BStringSet(bam[[1]][["qname"]]))


######################################################################
## 23. qa() report from bam file sequences
######################################################################
qaSummary <- qa(fq, lane="Roche 454")
perCycle  <- qaSummary[["perCycle"]]


######################################################################
## 24. Generate a quality by cycle plot
######################################################################
if(PNG) png(file.path(plotDir, "cycleQuality.png"))
ShortRead:::.plotCycleQuality(perCycle$quality, main=names(bam)[1])
if(PNG) dev.off()


######################################################################
## 25. Coverage plot ignoring strand orientation
######################################################################
IRanges   <- IRanges(start = bam[[1]][["pos"]], width=bam[[1]][["qwidth"]])
Cov       <- coverage(IRanges)
Peaks     <- slice(Cov, 0)

if(PNG) png(file.path(plotDir, "coverage.png"))
coverageplot(Peaks, main=names(bam)[1])
if(PNG) dev.off()


######################################################################
## 26. Coverage plot +/- strands
######################################################################
IRangesF   <- IRanges(start = bam[[1]][["pos"]][indPos], width=bam[[1]][["qwidth"]][indPos])
CovF       <- coverage(IRangesF)
PeaksF     <- slice(CovF, 0)

IRangesR   <- IRanges(start = bam[[1]][["pos"]][indNeg], width=bam[[1]][["qwidth"]][indNeg])
CovR       <- coverage(IRangesR)
PeaksR     <- slice(CovR, 0)


######################################################################
## 27. Coverage plot +/- strands visualized separately
######################################################################
if(PNG) png(file.path(plotDir, "coverageStrands.png"))
coverageplot(PeaksF, PeaksR,  main=names(bam)[1])
if(PNG) dev.off()
