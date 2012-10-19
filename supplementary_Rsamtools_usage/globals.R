#!/usr/bin/echo source me dont execute me

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
