#!/usr/bin/env Rscript

##Visualise resuls from Popoolation2 Fishers Exact Test
require("reshape2")
require("ggplot2")

## Create FET names
x             <- t(outer(1:6, 2:7, paste, sep="."))
FETnames      <- paste("FET_", x[!upper.tri(x)], sep="")
my_col_names  <- c("reference_contig","position","N_SNPs","fraction","avg_minimum", FETnames)
measure_vars  <- my_col_names[6:18]
id_vars       <- my_col_names[1:2]

## Read dataset
mydata        <- read.table("../30.popoolation/combined.fet", header=FALSE, col.names=my_col_names)
mydata[6:26]  <- lapply(mydata[6:26], function(x) as.numeric(substr(x,5,9)))

meltdata      <- melt(mydata,id.vars=id_vars,measure.vars=measure_vars)
FET_gene_plot <- ggplot(meltdata,aes(position,y=value)) + geom_point()
FET_gene_plot + facet_wrap(~ reference_contig) + opts(title="Fishers Exact Test") + ylab("-log 10 P val")

ggsave("FET_plot.png")

write.csv(meltdata,file="FET_data.csv")
